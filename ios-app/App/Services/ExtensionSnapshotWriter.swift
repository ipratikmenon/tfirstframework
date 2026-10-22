import Foundation
import SwiftData
import RecallCore

/// Keeps the shared snapshot fresh so the Share extension can recall things instantly.
@MainActor
enum ExtensionSnapshotWriter {
    static func write() {
        guard let url = AppGroup.snapshotURL else { return }
        let context = Persistence.shared.context
        let open = TodoStatus.open.rawValue
        let todos = (try? context.fetch(FetchDescriptor<TodoItem>(predicate: #Predicate { $0.statusRaw == open }))) ?? []
        let locked = MemoryStore.shared.lockedProfileIDs()

        var items: [ExtensionSnapshot.Item] = todos
            .filter { Set($0.profileIDs).isDisjoint(with: locked) }
            .map { .init(id: $0.id, kind: .todo, title: $0.title, category: $0.category, coordinate: $0.coordinate, date: $0.due) }

        // Recent places and upcoming events with a location.
        let since = Calendar.current.date(byAdding: .day, value: -180, to: Date()) ?? Date()
        let memories = (try? context.fetch(FetchDescriptor<MemoryEntry>(predicate: #Predicate { $0.date >= since && $0.latitude != nil }))) ?? []
        for m in memories where m.source == .location || m.source == .calendar || m.source == .glasses || m.source == .manual {
            guard Set(m.profileIDs).isDisjoint(with: locked) else { continue }
            items.append(.init(id: m.id, kind: m.source == .calendar ? .event : .memory, title: m.title, coordinate: m.coordinate, date: m.date))
        }
        let snapshot = ExtensionSnapshot(items: Array(items.prefix(1500)))
        if let data = try? JSONEncoder().encode(snapshot) {
            try? data.write(to: url, options: [.atomic, .completeFileProtectionUntilFirstUserAuthentication])
        }
    }
}
