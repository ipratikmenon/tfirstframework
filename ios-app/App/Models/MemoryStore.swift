import Foundation
import SwiftData
import RecallCore

/// All writes go through here so profile assignment, de-duplication and the
/// extension snapshot stay consistent.
@MainActor
final class MemoryStore {
    static let shared = MemoryStore()
    private var context: ModelContext { Persistence.shared.context }

    // MARK: profiles

    func profiles() -> [Profile] {
        let models = (try? context.fetch(FetchDescriptor<ProfileModel>(sortBy: [SortDescriptor(\.sortOrder), SortDescriptor(\.createdAt)]))) ?? []
        return models.map(\.profile)
    }

    func lockedProfileIDs() -> Set<UUID> {
        let models = (try? context.fetch(FetchDescriptor<ProfileModel>(predicate: #Predicate { $0.requiresAuth }))) ?? []
        return Set(models.map(\.id))
    }

    // MARK: memories

    /// Inserts a memory unless one with the same external id exists. Returns the stored entry.
    @discardableResult
    func upsert(_ entry: MemoryEntry) -> MemoryEntry {
        if let ext = entry.externalID, let existing = existing(externalID: ext) {
            // Keep user edits (profiles) but refresh content that can change at the source.
            existing.title = entry.title
            existing.body = entry.body
            existing.date = entry.date
            existing.endDate = entry.endDate
            if entry.latitude != nil { existing.latitude = entry.latitude; existing.longitude = entry.longitude }
            if entry.placeName != nil { existing.placeName = entry.placeName }
            existing.people = entry.people
            assignProfiles(existing)
            return existing
        }
        context.insert(entry)
        assignProfiles(entry)
        return entry
    }

    func existing(externalID ext: String) -> MemoryEntry? {
        var d = FetchDescriptor<MemoryEntry>(predicate: #Predicate { $0.externalID == ext })
        d.fetchLimit = 1
        return try? context.fetch(d).first
    }

    func hasEntry(externalID ext: String) -> Bool { existing(externalID: ext) != nil }

    func assignProfiles(_ entry: MemoryEntry) {
        guard !entry.profilesLocked else { return }
        entry.profileIDs = ProfileAssigner().assign(entry.record, profiles: profiles())
    }

    func assignProfiles(_ todo: TodoItem) {
        guard !todo.profilesLocked else { return }
        // A to-do is assigned like a memory made at creation time from its text.
        let record = MemoryRecord(date: todo.createdAt, source: .voice, title: todo.title,
                                  coordinate: todo.coordinate, placeName: todo.placeName,
                                  people: todo.person.map { [$0] } ?? [])
        todo.profileIDs = ProfileAssigner().assign(record, profiles: profiles())
    }

    /// Re-run automatic assignment after profiles or their rules changed.
    func reassignAll() {
        let profiles = profiles()
        let assigner = ProfileAssigner()
        for entry in (try? context.fetch(FetchDescriptor<MemoryEntry>())) ?? [] where !entry.profilesLocked {
            entry.profileIDs = assigner.assign(entry.record, profiles: profiles)
        }
        for todo in (try? context.fetch(FetchDescriptor<TodoItem>())) ?? [] { assignProfiles(todo) }
        save()
    }

    func save() {
        try? context.save()
        ExtensionSnapshotWriter.write()
    }

    // MARK: to-dos

    @discardableResult
    func addTodo(from draft: TodoDraft, status: TodoStatus) -> TodoItem {
        let todo = TodoItem(title: draft.title, due: draft.due, dueHasTime: draft.dueHasTime, category: draft.category,
                            placeName: draft.placeName, person: draft.person, status: status, confidence: draft.confidence)
        context.insert(todo)
        assignProfiles(todo)
        return todo
    }

    func openTodos() -> [TodoItem] {
        let open = TodoStatus.open.rawValue
        return (try? context.fetch(FetchDescriptor<TodoItem>(predicate: #Predicate { $0.statusRaw == open }))) ?? []
    }

    // MARK: deleting

    func delete(source: MemorySource) {
        let raw = source.rawValue
        for e in (try? context.fetch(FetchDescriptor<MemoryEntry>(predicate: #Predicate { $0.sourceRaw == raw }))) ?? [] { context.delete(e) }
        save()
    }

    func deleteDay(_ day: Date) {
        let cal = Calendar.current
        let start = cal.startOfDay(for: day)
        let end = cal.date(byAdding: .day, value: 1, to: start) ?? start
        for e in (try? context.fetch(FetchDescriptor<MemoryEntry>(predicate: #Predicate { $0.date >= start && $0.date < end }))) ?? [] { context.delete(e) }
        save()
    }

    /// Removes everything assigned to a profile (items also in other profiles are kept and just unassigned from it).
    func deleteData(of profile: Profile) {
        let id = profile.id.uuidString
        for e in (try? context.fetch(FetchDescriptor<MemoryEntry>())) ?? [] where e.profileIDStrings.contains(id) {
            if e.profileIDStrings.count == 1 { context.delete(e) } else { e.profileIDStrings.removeAll { $0 == id } }
        }
        for t in (try? context.fetch(FetchDescriptor<TodoItem>())) ?? [] where t.profileIDStrings.contains(id) {
            if t.profileIDStrings.count == 1 { context.delete(t) } else { t.profileIDStrings.removeAll { $0 == id } }
        }
        save()
    }

    func deleteEverything() {
        try? context.delete(model: MemoryEntry.self)
        try? context.delete(model: TodoItem.self)
        try? context.delete(model: DestinationRecord.self)
        NotificationService.shared.cancelAll()
        save()
    }

    // MARK: export

    struct ExportFile: Codable {
        struct Memory: Codable {
            var id: UUID; var date: Date; var endDate: Date?; var source: String; var title: String; var body: String
            var latitude: Double?; var longitude: Double?; var place: String?; var people: [String]; var profiles: [String]
        }
        struct Todo: Codable {
            var id: UUID; var title: String; var due: Date?; var status: String; var place: String?; var person: String?
            var created: Date; var profiles: [String]
        }
        var exportedAt: Date
        var profiles: [Profile]
        var memories: [Memory]
        var todos: [Todo]
    }

    /// Writes all data as JSON to a temporary file and returns its URL. `profile` limits the export to one profile.
    func exportJSON(profile: Profile? = nil) throws -> URL {
        let allProfiles = profiles()
        let names = Dictionary(uniqueKeysWithValues: allProfiles.map { ($0.id.uuidString, $0.name) })
        let entries = ((try? context.fetch(FetchDescriptor<MemoryEntry>(sortBy: [SortDescriptor(\.date)]))) ?? [])
            .filter { profile == nil || $0.profileIDStrings.contains(profile!.id.uuidString) }
        let todos = ((try? context.fetch(FetchDescriptor<TodoItem>(sortBy: [SortDescriptor(\.createdAt)]))) ?? [])
            .filter { profile == nil || $0.profileIDStrings.contains(profile!.id.uuidString) }
        let file = ExportFile(
            exportedAt: Date(),
            profiles: profile.map { [$0] } ?? allProfiles,
            memories: entries.map { .init(id: $0.id, date: $0.date, endDate: $0.endDate, source: $0.sourceRaw, title: $0.title,
                                          body: $0.body, latitude: $0.latitude, longitude: $0.longitude, place: $0.placeName,
                                          people: $0.people, profiles: $0.profileIDStrings.compactMap { names[$0] }) },
            todos: todos.map { .init(id: $0.id, title: $0.title, due: $0.due, status: $0.statusRaw, place: $0.placeName,
                                     person: $0.person, created: $0.createdAt, profiles: $0.profileIDStrings.compactMap { names[$0] }) })
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        encoder.dateEncodingStrategy = .iso8601
        let stamp = ISO8601DateFormatter().string(from: Date()).replacingOccurrences(of: ":", with: "-")
        let url = FileManager.default.temporaryDirectory.appendingPathComponent("TotalRecall-\(profile?.name ?? "All")-\(stamp).json")
        try encoder.encode(file).write(to: url, options: [.atomic, .completeFileProtection])
        return url
    }
}
