import Foundation

public enum CandidateKind: String, Codable, Sendable { case todo, memory, event }

/// Something that might be worth surfacing on a trip.
public struct RouteCandidate: Hashable, Identifiable, Sendable {
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

public struct NamedPlace: Hashable, Sendable {
    public var name: String
    public var coordinate: Coordinate
    public init(name: String, coordinate: Coordinate) { self.name = name; self.coordinate = coordinate }
}

public enum MatchReason: Hashable, Sendable {
    case atDestination
    case onRoute(distanceMeters: Double)
    case categoryAtDestination(PlaceCategory)
    case categoryAlongRoute(placeName: String, distanceMeters: Double)

    public var description: String {
        switch self {
        case .atDestination: return "At your destination"
        case .onRoute(let d): return "\(Self.metres(d)) from your route"
        case .categoryAtDestination(let c): return "Your destination is a \(c.displayName.lowercased())"
        case .categoryAlongRoute(let n, let d): return "\(n), \(Self.metres(d)) off your route"
        }
    }

    static func metres(_ d: Double) -> String { d < 1000 ? "\(Int((d / 10).rounded() * 10)) m" : String(format: "%.1f km", d / 1000) }
}

public struct RouteMatch: Hashable, Identifiable, Sendable {
    public var candidate: RouteCandidate
    public var reason: MatchReason
    public var score: Double
    public var id: UUID { candidate.id }
    public init(candidate: RouteCandidate, reason: MatchReason, score: Double) {
        self.candidate = candidate; self.reason = reason; self.score = score
    }
}

/// Decides which stored to-dos and memories are relevant to a planned trip.
public enum RouteMatcher {
    public static func match(candidates: [RouteCandidate], route: [Coordinate], destination: Coordinate,
                             destinationCategory: PlaceCategory? = nil,
                             nearbyPlaces: [PlaceCategory: [NamedPlace]] = [:],
                             maxOffRoute: Double = 400, destinationRadius: Double = 300,
                             now: Date = Date()) -> [RouteMatch] {
        var matches: [RouteMatch] = []
        let path = route.isEmpty ? [destination] : route
        for c in candidates {
            var found: (MatchReason, Double)?
            if let coord = c.coordinate {
                if Geo.distance(coord, destination) <= destinationRadius {
                    found = (.atDestination, 1.0)
                } else {
                    let d = Geo.distance(from: coord, toPath: path)
                    if d <= maxOffRoute { found = (.onRoute(distanceMeters: d), 0.9 - 0.4 * d / maxOffRoute) }
                }
            } else if let category = c.category {
                if destinationCategory == category {
                    found = (.categoryAtDestination(category), 0.95)
                } else if let places = nearbyPlaces[category] {
                    let best = places.map { ($0, Geo.distance(from: $0.coordinate, toPath: path)) }
                        .filter { $0.1 <= maxOffRoute }.min { $0.1 < $1.1 }
                    if let (place, d) = best {
                        found = (.categoryAlongRoute(placeName: place.name, distanceMeters: d), 0.8 - 0.4 * d / maxOffRoute)
                    }
                }
            }
            guard let (reason, initialScore) = found else { continue }
            var score = initialScore
            if c.kind == .todo { score += 0.05 }
            if let due = c.date, due.timeIntervalSince(now) < 48 * 3600, due.timeIntervalSince(now) > -86400 { score += 0.1 }
            matches.append(RouteMatch(candidate: c, reason: reason, score: score))
        }
        return matches.sorted { $0.score > $1.score }
    }
}
