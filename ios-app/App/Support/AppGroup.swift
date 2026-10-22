import Foundation
import RecallCore

/// Shared with the Share extension.
enum AppGroup {
    static let identifier = "group.com.example.totalrecall"

    static var containerURL: URL? {
        FileManager.default.containerURL(forSecurityApplicationGroupIdentifier: identifier)
    }

    static var inboxURL: URL? {
        guard let base = containerURL else { return nil }
        let url = base.appendingPathComponent("Inbox", isDirectory: true)
        try? FileManager.default.createDirectory(at: url, withIntermediateDirectories: true)
        return url
    }

    static var snapshotURL: URL? { containerURL?.appendingPathComponent("snapshot.json") }
}

/// One thing shared into the app from the system share sheet.
struct InboxItem: Codable {
    var id: UUID
    var date: Date
    var text: String?
    var url: String?
    var imageFile: String?
}
