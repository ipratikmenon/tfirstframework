import Foundation

/// A compact copy of what is worth recalling, written by the app into the shared App Group
/// so the Share extension can answer "anything relevant near this destination?" instantly,
/// without opening the app or touching the database.
public struct ExtensionSnapshot: Codable, Sendable {
    public struct Item: Codable, Hashable, Sendable {
        public var id: UUID
        public var kind: CandidateKind
        public var title: String
        public var category: PlaceCategory?
        public var coordinate: Coordinate?
        public var date: Date?
        public init(id: UUID, kind: CandidateKind, title: String, category: PlaceCategory? = nil,
                    coordinate: Coordinate? = nil, date: Date? = nil) {
            self.id = id; self.kind = kind; self.title = title
            self.category = category; self.coordinate = coordinate; self.date = date
        }
    }

    public var generatedAt: Date
    public var items: [Item]
    public init(generatedAt: Date = Date(), items: [Item]) {
        self.generatedAt = generatedAt
        self.items = items
    }

    public var candidates: [RouteCandidate] {
        items.map { RouteCandidate(id: $0.id, kind: $0.kind, title: $0.title, category: $0.category,
                                   coordinate: $0.coordinate, date: $0.date) }
    }
}

/// Destination-only matching used when no route is available (for example from the Share extension).
public enum QuickDestinationMatcher {
    public static func match(destination: SharedDestination, snapshot: ExtensionSnapshot,
                             radius: Double = 600, now: Date = Date()) -> [RouteMatch] {
        let category = PlaceCategory.detect(in: [destination.name, destination.address].compactMap { $0 }.joined(separator: " "))
        guard let coordinate = destination.coordinate else {
            // Without coordinates only a category match is possible.
            guard let category else { return [] }
            return snapshot.candidates.filter { $0.category == category && $0.coordinate == nil }
                .map { RouteMatch(candidate: $0, reason: .categoryAtDestination(category), score: 0.9) }
        }
        return RouteMatcher.match(candidates: snapshot.candidates, route: [coordinate], destination: coordinate,
                                  destinationCategory: category, maxOffRoute: 0, destinationRadius: radius, now: now)
    }

    /// One-line text for a notification.
    public static func summary(for matches: [RouteMatch], destinationName: String?) -> String? {
        guard !matches.isEmpty else { return nil }
        let titles = matches.prefix(3).map { $0.candidate.title }
        let more = matches.count > 3 ? " and \(matches.count - 3) more" : ""
        return titles.joined(separator: ", ") + more
    }
}
