import Foundation
import SwiftData
import RecallCore

@MainActor
enum HistoryQueries {
    static func visits(since: Date) -> [Visit] {
        let location = MemorySource.location.rawValue
        let d = FetchDescriptor<MemoryEntry>(predicate: #Predicate { $0.sourceRaw == location && $0.date >= since && $0.latitude != nil })
        return ((try? Persistence.shared.context.fetch(d)) ?? []).compactMap { e in
            e.coordinate.map { Visit(coordinate: $0, arrival: e.date, departure: e.endDate) }
        }
    }

    static func records(since: Date? = nil) -> [MemoryRecord] {
        var d = FetchDescriptor<MemoryEntry>(sortBy: [SortDescriptor(\.date, order: .reverse)])
        if let since { d.predicate = #Predicate { $0.date >= since } }
        return ((try? Persistence.shared.context.fetch(d)) ?? []).map(\.record)
    }
}
