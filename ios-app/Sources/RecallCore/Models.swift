import Foundation

/// Where a memory came from. Every source is opt-in in the app.
public enum MemorySource: String, Codable, CaseIterable, Sendable {
    case location, photo, calendar, health, mail, share, voice, glasses, manual

    public var displayName: String {
        switch self {
        case .location: return "Place"
        case .photo: return "Photo"
        case .calendar: return "Event"
        case .health: return "Health"
        case .mail: return "Mail"
        case .share: return "Shared"
        case .voice: return "Voice"
        case .glasses: return "Glasses"
        case .manual: return "Note"
        }
    }

    /// SF Symbol name used by the app.
    public var symbolName: String {
        switch self {
        case .location: return "location.fill"
        case .photo: return "photo.fill"
        case .calendar: return "calendar"
        case .health: return "heart.fill"
        case .mail: return "envelope.fill"
        case .share: return "square.and.arrow.down.fill"
        case .voice: return "waveform"
        case .glasses: return "eyeglasses"
        case .manual: return "note.text"
        }
    }
}

public struct Coordinate: Codable, Hashable, Sendable {
    public var latitude: Double
    public var longitude: Double
    public init(latitude: Double, longitude: Double) {
        self.latitude = latitude
        self.longitude = longitude
    }
}

/// Plain value view of a stored memory. The app maps its SwiftData models to
/// this type so the logic in RecallCore stays testable without Apple frameworks.
public struct MemoryRecord: Codable, Hashable, Identifiable, Sendable {
    public var id: UUID
    public var date: Date
    public var source: MemorySource
    public var title: String
    public var body: String
    public var coordinate: Coordinate?
    public var placeName: String?
    public var people: [String]
    public var accountLabel: String?
    public var calendarName: String?
    public var profileIDs: [UUID]

    public init(id: UUID = UUID(), date: Date, source: MemorySource, title: String,
                body: String = "", coordinate: Coordinate? = nil, placeName: String? = nil,
                people: [String] = [], accountLabel: String? = nil, calendarName: String? = nil,
                profileIDs: [UUID] = []) {
        self.id = id
        self.date = date
        self.source = source
        self.title = title
        self.body = body
        self.coordinate = coordinate
        self.placeName = placeName
        self.people = people
        self.accountLabel = accountLabel
        self.calendarName = calendarName
        self.profileIDs = profileIDs
    }

    public var searchableText: String {
        [title, body, placeName ?? "", people.joined(separator: " ")].joined(separator: " ")
    }
}
