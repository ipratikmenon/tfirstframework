import Foundation

public struct TodoRecord: Hashable, Identifiable, Sendable {
    public var id: UUID
    public var title: String
    public var due: Date?
    public var dueHasTime: Bool
    public var category: PlaceCategory?
    public var coordinate: Coordinate?
    public var placeName: String?
    public var profileIDs: [UUID]
    public var lastAlerted: Date?

    public init(id: UUID = UUID(), title: String, due: Date? = nil, dueHasTime: Bool = false,
                category: PlaceCategory? = nil, coordinate: Coordinate? = nil, placeName: String? = nil,
                profileIDs: [UUID] = [], lastAlerted: Date? = nil) {
        self.id = id; self.title = title; self.due = due; self.dueHasTime = dueHasTime
        self.category = category; self.coordinate = coordinate; self.placeName = placeName
        self.profileIDs = profileIDs; self.lastAlerted = lastAlerted
    }
}

public enum AlertReason: Hashable, Sendable {
    case due
    case likelyPlace(String)
}

public struct PlannedAlert: Hashable, Identifiable, Sendable {
    public var id: UUID            // the to-do id
    public var title: String
    public var body: String
    public var fireDate: Date
    public var reason: AlertReason
    public var score: Double
}

/// Chooses which alerts to schedule: due dates plus places the routine model expects the user to be.
public enum AlertPlanner {
    public static func plan(todos: [TodoRecord], routine: RoutineModel?, knownPlaces: [PlaceCategory: [NamedPlace]],
                            now: Date, horizonHours: Int = 6, dailyCap: Int = 5, alreadyAlertedToday: Int = 0,
                            filter: ProfileFilter = .all, cooldown: TimeInterval = 4 * 3600) -> [PlannedAlert] {
        var alerts: [UUID: PlannedAlert] = [:]
        let horizon = now.addingTimeInterval(Double(horizonHours) * 3600)
        let predictions = routine?.predictions(from: now, hours: horizonHours) ?? []

        for todo in todos where filter.includes(todo.profileIDs) {
            if let last = todo.lastAlerted, now.timeIntervalSince(last) < cooldown { continue }

            if let due = todo.due, due <= horizon, due > now.addingTimeInterval(-3600) {
                let lead: TimeInterval = todo.dueHasTime ? 30 * 60 : 0
                let fire = max(now.addingTimeInterval(60), due.addingTimeInterval(-lead))
                alerts[todo.id] = PlannedAlert(id: todo.id, title: todo.title,
                                               body: todo.dueHasTime ? "Due soon" : "Due today",
                                               fireDate: fire, reason: .due, score: 0.9)
            }

            for p in predictions {
                var placeName: String?
                if let c = todo.coordinate, Geo.distance(c, p.coordinate) <= 500 {
                    placeName = todo.placeName ?? "this area"
                } else if let cat = todo.category, todo.coordinate == nil,
                          let place = knownPlaces[cat]?.first(where: { Geo.distance($0.coordinate, p.coordinate) <= 600 }) {
                    placeName = place.name
                }
                guard let name = placeName else { continue }
                let score = 0.5 + 0.4 * p.probability
                if let existing = alerts[todo.id], existing.score >= score { continue }
                alerts[todo.id] = PlannedAlert(id: todo.id, title: todo.title, body: "You usually pass \(name) around now",
                                               fireDate: max(now.addingTimeInterval(60), p.expectedArrival.addingTimeInterval(-15 * 60)),
                                               reason: .likelyPlace(name), score: score)
            }
        }
        let room = max(0, dailyCap - alreadyAlertedToday)
        return Array(alerts.values.sorted { $0.score > $1.score }.prefix(room))
    }
}
