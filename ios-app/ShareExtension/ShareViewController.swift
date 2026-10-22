import UIKit
import UniformTypeIdentifiers
import UserNotifications
import ImageIO
import RecallCore

/// Total Recall in the system share sheet. It saves what was shared into the App Group inbox and,
/// when the item is a place from Maps, immediately tells the user what they have noted near it.
final class ShareViewController: UIViewController {
    private let label = UILabel()
    private let icon = UIImageView()

    override func viewDidLoad() {
        super.viewDidLoad()
        view.backgroundColor = .systemBackground
        icon.image = UIImage(systemName: "checkmark.circle.fill")
        icon.tintColor = .systemGreen
        icon.contentMode = .scaleAspectFit
        label.text = "Saving\u{2026}"
        label.font = .preferredFont(forTextStyle: .headline)
        label.textAlignment = .center
        label.numberOfLines = 0
        label.adjustsFontForContentSizeCategory = true

        let stack = UIStackView(arrangedSubviews: [icon, label])
        stack.axis = .vertical
        stack.spacing = 12
        stack.alignment = .center
        stack.translatesAutoresizingMaskIntoConstraints = false
        view.addSubview(stack)
        NSLayoutConstraint.activate([
            icon.heightAnchor.constraint(equalToConstant: 56), icon.widthAnchor.constraint(equalToConstant: 56),
            stack.centerXAnchor.constraint(equalTo: view.centerXAnchor),
            stack.centerYAnchor.constraint(equalTo: view.centerYAnchor),
            stack.leadingAnchor.constraint(greaterThanOrEqualTo: view.leadingAnchor, constant: 24),
            stack.trailingAnchor.constraint(lessThanOrEqualTo: view.trailingAnchor, constant: -24),
        ])
    }

    override func viewDidAppear(_ animated: Bool) {
        super.viewDidAppear(animated)
        Task { await process() }
    }

    private func process() async {
        var item = InboxItem(id: UUID(), date: Date(), text: nil, url: nil, imageFile: nil)
        let providers = (extensionContext?.inputItems as? [NSExtensionItem])?.flatMap { $0.attachments ?? [] } ?? []

        for provider in providers {
            if provider.hasItemConformingToTypeIdentifier(UTType.url.identifier),
               let value = try? await provider.loadItem(forTypeIdentifier: UTType.url.identifier), let url = value as? URL {
                item.url = url.absoluteString
            } else if provider.hasItemConformingToTypeIdentifier(UTType.image.identifier),
                      let value = try? await provider.loadItem(forTypeIdentifier: UTType.image.identifier) {
                item.imageFile = saveImage(value, id: item.id)
            } else if provider.hasItemConformingToTypeIdentifier(UTType.plainText.identifier),
                      let value = try? await provider.loadItem(forTypeIdentifier: UTType.plainText.identifier), let text = value as? String {
                if text.hasPrefix("http"), item.url == nil { item.url = text } else { item.text = (item.text.map { $0 + "\n" } ?? "") + text }
            }
        }
        // Maps and Safari often attach the page title as text next to the URL.
        if let content = (extensionContext?.inputItems as? [NSExtensionItem])?.first?.attributedContentText?.string, !content.isEmpty, item.text == nil {
            item.text = content
        }

        write(item)
        var message = "Saved to Total Recall"
        if let urlString = item.url, let destination = MapsLinkParser.parse(urlString), let recall = await recall(for: destination) {
            message = recall
        }
        label.text = message
        UIImpactFeedbackGenerator(style: .light).impactOccurred()
        try? await Task.sleep(nanoseconds: 1_400_000_000)
        extensionContext?.completeRequest(returningItems: nil)
    }

    // MARK: inbox

    private func write(_ item: InboxItem) {
        guard let inbox = AppGroup.inboxURL else { return }
        let encoder = JSONEncoder()
        encoder.dateEncodingStrategy = .iso8601
        if let data = try? encoder.encode(item) {
            try? data.write(to: inbox.appendingPathComponent("\(item.id.uuidString).json"), options: .atomic)
        }
    }

    private func saveImage(_ value: NSSecureCoding, id: UUID) -> String? {
        guard let inbox = AppGroup.inboxURL else { return nil }
        var data: Data?
        if let url = value as? URL { data = try? Data(contentsOf: url) }
        else if let image = value as? UIImage { data = image.jpegData(compressionQuality: 0.8) }
        else if let d = value as? Data { data = d }
        guard let data, let source = CGImageSourceCreateWithData(data as CFData, nil) else { return nil }
        // Downscale: extensions have a small memory limit and the app only needs enough for text recognition.
        let options: [CFString: Any] = [kCGImageSourceCreateThumbnailFromImageAlways: true,
                                        kCGImageSourceCreateThumbnailWithTransform: true,
                                        kCGImageSourceThumbnailMaxPixelSize: 1600]
        guard let thumb = CGImageSourceCreateThumbnailAtIndex(source, 0, options as CFDictionary),
              let jpeg = UIImage(cgImage: thumb).jpegData(compressionQuality: 0.8) else { return nil }
        let name = "\(id.uuidString).jpg"
        try? jpeg.write(to: inbox.appendingPathComponent(name), options: .atomic)
        return name
    }

    // MARK: instant recall for a shared destination

    private func recall(for destination: SharedDestination) async -> String? {
        guard let url = AppGroup.snapshotURL, let data = try? Data(contentsOf: url),
              let snapshot = try? JSONDecoder().decode(ExtensionSnapshot.self, from: data) else { return nil }
        let matches = QuickDestinationMatcher.match(destination: destination, snapshot: snapshot)
        guard let summary = QuickDestinationMatcher.summary(for: matches, destinationName: destination.name) else { return nil }

        let content = UNMutableNotificationContent()
        content.title = "Before you go to \(destination.name ?? "your destination")"
        content.body = summary
        content.sound = .default
        let center = UNUserNotificationCenter.current()
        let settings = await center.notificationSettings()
        if settings.authorizationStatus == .authorized {
            try? await center.add(UNNotificationRequest(identifier: UUID().uuidString, content: content,
                                                        trigger: UNTimeIntervalNotificationTrigger(timeInterval: 1, repeats: false)))
        }
        return "Remember: " + summary
    }
}
