import Foundation
import SwiftData
import Vision
import UIKit
import RecallCore

/// Picks up items dropped into the App Group by the Share extension.
@MainActor
enum ShareInbox {
    static func ingest() async {
        guard let inbox = AppGroup.inboxURL,
              let files = try? FileManager.default.contentsOfDirectory(at: inbox, includingPropertiesForKeys: nil) else { return }
        let decoder = JSONDecoder()
        decoder.dateDecodingStrategy = .iso8601
        for file in files where file.pathExtension == "json" {
            defer { try? FileManager.default.removeItem(at: file) }
            guard let data = try? Data(contentsOf: file), let item = try? decoder.decode(InboxItem.self, from: data) else { continue }
            await handle(item, inbox: inbox)
        }
        MemoryStore.shared.save()
    }

    private static func handle(_ item: InboxItem, inbox: URL) async {
        var title = "Shared item"
        var body = item.text ?? ""
        var coordinate: Coordinate?
        var place: String?

        if let urlString = item.url {
            if let dest = MapsLinkParser.parse(urlString) {
                // A destination shared from Maps becomes a place memory and a suggestion for trip recall.
                title = dest.name ?? dest.address ?? "Shared place"
                coordinate = dest.coordinate
                place = dest.name ?? dest.address
                if let c = dest.coordinate { recordDestination(name: title, coordinate: c) }
            } else {
                title = URL(string: urlString)?.host ?? "Link"
                body = [item.text, urlString].compactMap { $0 }.joined(separator: "\n")
            }
        } else if let text = item.text {
            title = String(text.split(separator: "\n").first?.prefix(80) ?? "Note")
        }

        if let imageName = item.imageFile {
            let imageURL = inbox.appendingPathComponent(imageName)
            defer { try? FileManager.default.removeItem(at: imageURL) }
            title = "Shared image"
            if let data = try? Data(contentsOf: imageURL), let image = UIImage(data: data)?.cgImage {
                let text = await PhotoService.analyse(image, wantsLabels: true)
                if !text.isEmpty { body = text; title = String(text.prefix(60)) }
            }
        }

        let entry = MemoryEntry(date: item.date, source: .share, title: title, body: body, coordinate: coordinate,
                                placeName: place, externalID: "share-\(item.id.uuidString)")
        MemoryStore.shared.upsert(entry)

        // Spoken-style commitments in shared text (a forwarded message, for example) become drafts.
        if !body.isEmpty {
            for draft in TodoExtractor().extract(from: body) where draft.confidence >= 0.7 {
                MemoryStore.shared.addTodo(from: draft, status: .draft)
            }
        }
    }

    private static func recordDestination(name: String, coordinate: Coordinate) {
        let context = Persistence.shared.context
        let all = (try? context.fetch(FetchDescriptor<DestinationRecord>())) ?? []
        if let existing = all.first(where: { Geo.distance($0.coordinate, coordinate) < 100 }) {
            existing.lastUsed = Date()
            existing.useCount += 1
        } else {
            context.insert(DestinationRecord(name: name, coordinate: coordinate))
        }
    }
}
