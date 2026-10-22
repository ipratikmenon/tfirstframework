import Foundation
import Photos
import Vision
import UIKit
import RecallCore

/// Indexes the photo library: when and where each photo was taken, plus text found in screenshots
/// and labels for photos from the glasses album. Images are analysed on device and never copied.
@MainActor
final class PhotoService: ObservableObject {
    static let shared = PhotoService()

    @Published private(set) var authorization = PHPhotoLibrary.authorizationStatus(for: .readWrite)
    var isEnabled: Bool { UserDefaults.standard.bool(forKey: Prefs.enabled("photo")) }
    var isAuthorized: Bool { authorization == .authorized || authorization == .limited }

    func requestAuthorization() async {
        let status = await PHPhotoLibrary.requestAuthorization(for: .readWrite)
        authorization = status
    }

    struct Album: Identifiable, Hashable {
        let id: String
        let title: String
    }

    /// All albums the user could mark as the Ray-Ban Meta album.
    func albums() -> [Album] {
        guard isAuthorized else { return [] }
        var out: [Album] = []
        let result = PHAssetCollection.fetchAssetCollections(with: .album, subtype: .any, options: nil)
        result.enumerateObjects { c, _, _ in out.append(Album(id: c.localIdentifier, title: c.localizedTitle ?? "Untitled")) }
        return out.sorted { $0.title.localizedCaseInsensitiveCompare($1.title) == .orderedAscending }
    }

    func sync(maxNew: Int = 300) async {
        guard isEnabled, isAuthorized else { return }
        let defaults = UserDefaults.standard
        let last = defaults.object(forKey: Prefs.lastPhotoSync) as? Date
            ?? Calendar.current.date(byAdding: .day, value: -14, to: Date()) ?? Date()
        let options = PHFetchOptions()
        options.predicate = NSPredicate(format: "creationDate > %@", last as NSDate)
        options.sortDescriptors = [NSSortDescriptor(key: "creationDate", ascending: false)]
        options.fetchLimit = maxNew

        let glassesAssetIDs = glassesAssetIdentifiers()
        let assets = PHAsset.fetchAssets(with: options)
        var analysed = 0
        var newest = last
        var pending: [PHAsset] = []
        assets.enumerateObjects { a, _, _ in pending.append(a) }

        for asset in pending {
            guard let created = asset.creationDate else { continue }
            newest = max(newest, created)
            let isGlasses = glassesAssetIDs.contains(asset.localIdentifier)
            let isScreenshot = asset.mediaSubtypes.contains(.photoScreenshot)
            let source: MemorySource = isGlasses ? .glasses : .photo
            let title = isGlasses ? "Glasses photo" : isScreenshot ? "Screenshot" : asset.mediaType == .video ? "Video" : "Photo"
            var body = ""
            if (isScreenshot || isGlasses || analysed < 25) && asset.mediaType == .image && analysed < 40 {
                if let image = await loadImage(for: asset) {
                    let result = await Self.analyse(image, wantsLabels: !isScreenshot)
                    body = result
                    analysed += 1
                }
            }
            let coordinate = asset.location.map { Coordinate(latitude: $0.coordinate.latitude, longitude: $0.coordinate.longitude) }
            let entry = MemoryEntry(date: created, source: source, title: title, body: body, coordinate: coordinate,
                                    externalID: "photo-\(asset.localIdentifier)")
            MemoryStore.shared.upsert(entry)
        }
        defaults.set(newest, forKey: Prefs.lastPhotoSync)
        MemoryStore.shared.save()
    }

    private func glassesAssetIdentifiers() -> Set<String> {
        guard let albumID = UserDefaults.standard.string(forKey: Prefs.glassesAlbumID), !albumID.isEmpty,
              let collection = PHAssetCollection.fetchAssetCollections(withLocalIdentifiers: [albumID], options: nil).firstObject
        else { return [] }
        var ids = Set<String>()
        PHAsset.fetchAssets(in: collection, options: nil).enumerateObjects { a, _, _ in ids.insert(a.localIdentifier) }
        return ids
    }

    private func loadImage(for asset: PHAsset) async -> CGImage? {
        await withCheckedContinuation { continuation in
            let options = PHImageRequestOptions()
            options.deliveryMode = .highQualityFormat
            options.isNetworkAccessAllowed = false     // never pull from iCloud just to index
            options.isSynchronous = false
            var resumed = false
            PHImageManager.default().requestImage(for: asset, targetSize: CGSize(width: 1024, height: 1024),
                                                  contentMode: .aspectFit, options: options) { image, info in
                let degraded = (info?[PHImageResultIsDegradedKey] as? Bool) ?? false
                if degraded { return }
                guard !resumed else { return }
                resumed = true
                continuation.resume(returning: image?.cgImage)
            }
        }
    }

    /// Text recognised in the image plus (optionally) the strongest scene labels.
    nonisolated static func analyse(_ image: CGImage, wantsLabels: Bool) async -> String {
        await withCheckedContinuation { continuation in
            DispatchQueue.global(qos: .utility).async {
                let handler = VNImageRequestHandler(cgImage: image, options: [:])
                let text = VNRecognizeTextRequest()
                text.recognitionLevel = .accurate
                text.usesLanguageCorrection = true
                var requests: [VNRequest] = [text]
                let classify = VNClassifyImageRequest()
                if wantsLabels { requests.append(classify) }
                try? handler.perform(requests)

                var parts: [String] = []
                if wantsLabels {
                    let labels = (classify.results ?? []).filter { $0.confidence > 0.35 }.prefix(5)
                        .map { $0.identifier.replacingOccurrences(of: "_", with: " ") }
                    if !labels.isEmpty { parts.append("Looks like: " + labels.joined(separator: ", ")) }
                }
                let lines = (text.results ?? []).compactMap { $0.topCandidates(1).first?.string }
                if !lines.isEmpty { parts.append(lines.joined(separator: " ").prefix(600).description) }
                continuation.resume(returning: parts.joined(separator: ". "))
            }
        }
    }
}
