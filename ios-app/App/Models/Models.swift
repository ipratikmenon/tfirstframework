import Foundation
import SwiftData
import RecallCore

@Model
final class MemoryEntry {
    @Attribute(.unique) var id: UUID
    var date: Date
    var endDate: Date?
    var sourceRaw: String
    var title: String
    var body: String
    var latitude: Double?
    var longitude: Double?
    var placeName: String?
    var people: [String]
    var accountLabel: String?
    var calendarName: String?
    var profileIDStrings: [String]
    /// True once the user has chosen profiles by hand; automatic assignment then leaves it alone.
    var profilesLocked: Bool
    var externalID: String?
    var createdAt: Date

    init(id: UUID = UUID(), date: Date, endDate: Date? = nil, source: MemorySource, title: String, body: String = "",
         coordinate: Coordinate? = nil, placeName: String? = nil, people: [String] = [],
         accountLabel: String? = nil, calendarName: String? = nil, externalID: String? = nil) {
        self.id = id
        self.date = date
        self.endDate = endDate
        self.sourceRaw = source.rawValue
        self.title = title
        self.body = body
        self.latitude = coordinate?.latitude
        self.longitude = coordinate?.longitude
        self.placeName = placeName
        self.people = people
        self.accountLabel = accountLabel
        self.calendarName = calendarName
        self.profileIDStrings = []
        self.profilesLocked = false
        self.externalID = externalID
        self.createdAt = Date()
    }

    var source: MemorySource {
        get { MemorySource(rawValue: sourceRaw) ?? .manual }
        set { sourceRaw = newValue.rawValue }
    }

    var coordinate: Coordinate? {
        get {
            guard let latitude, let longitude else { return nil }
            return Coordinate(latitude: latitude, longitude: longitude)
        }
        set {
            latitude = newValue?.latitude
            longitude = newValue?.longitude
        }
    }

    var profileIDs: [UUID] {
        get { profileIDStrings.compactMap(UUID.init(uuidString:)) }
        set { profileIDStrings = newValue.map(\.uuidString) }
    }

    var record: MemoryRecord {
        MemoryRecord(id: id, date: date, source: source, title: title, body: body, coordinate: coordinate,
                     placeName: placeName, people: people, accountLabel: accountLabel,
                     calendarName: calendarName, profileIDs: profileIDs)
    }
}

enum TodoStatus: String, Codable {
    case draft      // heard, waiting for the user to confirm
    case open
    case done
    case dismissed
}

@Model
final class TodoItem {
    @Attribute(.unique) var id: UUID
    var title: String
    var due: Date?
    var dueHasTime: Bool
    var categoryRaw: String?
    var placeName: String?
    var latitude: Double?
    var longitude: Double?
    var person: String?
    var statusRaw: String
    var confidence: Double
    var createdAt: Date
    var completedAt: Date?
    var lastAlerted: Date?
    var profileIDStrings: [String]
    var profilesLocked: Bool

    init(id: UUID = UUID(), title: String, due: Date? = nil, dueHasTime: Bool = false, category: PlaceCategory? = nil,
         placeName: String? = nil, coordinate: Coordinate? = nil, person: String? = nil,
         status: TodoStatus = .open, confidence: Double = 1) {
        self.id = id
        self.title = title
        self.due = due
        self.dueHasTime = dueHasTime
        self.categoryRaw = category?.rawValue
        self.placeName = placeName
        self.latitude = coordinate?.latitude
        self.longitude = coordinate?.longitude
        self.person = person
        self.statusRaw = status.rawValue
        self.confidence = confidence
        self.createdAt = Date()
        self.profileIDStrings = []
        self.profilesLocked = false
    }

    var status: TodoStatus {
        get { TodoStatus(rawValue: statusRaw) ?? .open }
        set { statusRaw = newValue.rawValue }
    }

    var category: PlaceCategory? {
        get { categoryRaw.flatMap(PlaceCategory.init(rawValue:)) }
        set { categoryRaw = newValue?.rawValue }
    }

    var coordinate: Coordinate? {
        get {
            guard let latitude, let longitude else { return nil }
            return Coordinate(latitude: latitude, longitude: longitude)
        }
        set {
            latitude = newValue?.latitude
            longitude = newValue?.longitude
        }
    }

    var profileIDs: [UUID] {
        get { profileIDStrings.compactMap(UUID.init(uuidString:)) }
        set { profileIDStrings = newValue.map(\.uuidString) }
    }

    var record: TodoRecord {
        TodoRecord(id: id, title: title, due: due, dueHasTime: dueHasTime, category: category,
                   coordinate: coordinate, placeName: placeName, profileIDs: profileIDs, lastAlerted: lastAlerted)
    }
}

@Model
final class ProfileModel {
    @Attribute(.unique) var id: UUID
    var name: String
    var kindRaw: String
    var symbolName: String
    var rulesData: Data
    var threshold: Double
    var requiresAuth: Bool
    var sortOrder: Int
    var createdAt: Date

    init(profile: Profile, requiresAuth: Bool = false, sortOrder: Int = 0) {
        self.id = profile.id
        self.name = profile.name
        self.kindRaw = profile.kind.rawValue
        self.symbolName = profile.symbolName
        self.rulesData = (try? JSONEncoder().encode(profile.rules)) ?? Data()
        self.threshold = profile.threshold
        self.requiresAuth = requiresAuth
        self.sortOrder = sortOrder
        self.createdAt = Date()
    }

    var profile: Profile {
        get {
            let rules = (try? JSONDecoder().decode([ProfileRule].self, from: rulesData)) ?? []
            return Profile(id: id, name: name, kind: ProfileKind(rawValue: kindRaw) ?? .custom,
                           symbolName: symbolName, rules: rules, threshold: threshold)
        }
        set {
            name = newValue.name
            kindRaw = newValue.kind.rawValue
            symbolName = newValue.symbolName
            rulesData = (try? JSONEncoder().encode(newValue.rules)) ?? Data()
            threshold = newValue.threshold
        }
    }
}

/// Places the user has asked to be reminded about before, used for suggestions.
@Model
final class DestinationRecord {
    @Attribute(.unique) var id: UUID
    var name: String
    var latitude: Double
    var longitude: Double
    var lastUsed: Date
    var useCount: Int

    init(name: String, coordinate: Coordinate) {
        self.id = UUID()
        self.name = name
        self.latitude = coordinate.latitude
        self.longitude = coordinate.longitude
        self.lastUsed = Date()
        self.useCount = 1
    }

    var coordinate: Coordinate { Coordinate(latitude: latitude, longitude: longitude) }
}
