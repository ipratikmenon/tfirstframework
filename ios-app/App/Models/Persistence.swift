import Foundation
import SwiftData
import RecallCore
import os

/// One shared SwiftData container. The store is protected until the first unlock after boot,
/// so background captures (visits, refresh tasks) can still write while the phone is locked.
@MainActor
final class Persistence {
    static let shared = Persistence()

    let container: ModelContainer
    var context: ModelContext { container.mainContext }
    private let log = Logger(subsystem: "com.example.totalrecall", category: "persistence")

    private init() {
        let schema = Schema([MemoryEntry.self, TodoItem.self, ProfileModel.self, DestinationRecord.self])
        do {
            let support = try FileManager.default.url(for: .applicationSupportDirectory, in: .userDomainMask,
                                                      appropriateFor: nil, create: true)
            try? FileManager.default.setAttributes([.protectionKey: FileProtectionType.completeUntilFirstUserAuthentication],
                                                   ofItemAtPath: support.path)
            let config = ModelConfiguration(schema: schema, url: support.appendingPathComponent("TotalRecall.store"))
            container = try ModelContainer(for: schema, configurations: config)
        } catch {
            // Never lose the app to a storage failure: fall back to memory and say so in the log.
            log.error("Persistent store failed (\(error.localizedDescription)); using in-memory store")
            container = try! ModelContainer(for: schema, configurations: ModelConfiguration(isStoredInMemoryOnly: true))
        }
        seedDefaultProfile()
    }

    private func seedDefaultProfile() {
        let existing = (try? context.fetchCount(FetchDescriptor<ProfileModel>())) ?? 0
        guard existing == 0 else { return }
        context.insert(ProfileModel(profile: .defaultWork(), sortOrder: 0))
        try? context.save()
    }
}
