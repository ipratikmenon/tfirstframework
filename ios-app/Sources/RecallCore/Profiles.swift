import Foundation

/// A signal that suggests a memory belongs to a profile.
public enum ProfileRule: Codable, Hashable, Sendable {
    case source(MemorySource)
    case accountContains(String)
    case calendarNamed(String)
    case place(center: Coordinate, radiusMeters: Double, name: String)
    /// `weekdays` use Calendar numbering (1 = Sunday). Minutes are minutes after midnight.
    case timeWindow(weekdays: [Int], startMinute: Int, endMinute: Int)
    case person(String)
    case keyword(String)

    /// Strong signals (an account or calendar) decide alone; weak ones (time of day) only add up.
    public var weight: Double {
        switch self {
        case .accountContains, .calendarNamed: return 1.0
        case .place: return 0.8
        case .person, .keyword: return 0.6
        case .source: return 0.5
        case .timeWindow: return 0.4
        }
    }

    public var summary: String {
        switch self {
        case .source(let s): return "Source is \(s.displayName)"
        case .accountContains(let a): return "Account contains \u{201C}\(a)\u{201D}"
        case .calendarNamed(let c): return "Calendar \u{201C}\(c)\u{201D}"
        case .place(_, let r, let n): return "Within \(Int(r)) m of \(n)"
        case .timeWindow(let d, let s, let e): return "\(Self.dayNames(d)) \(Self.clock(s))\u{2013}\(Self.clock(e))"
        case .person(let p): return "Involves \(p)"
        case .keyword(let k): return "Mentions \u{201C}\(k)\u{201D}"
        }
    }

    static func clock(_ minute: Int) -> String {
        String(format: "%02d:%02d", minute / 60, minute % 60)
    }

    static func dayNames(_ days: [Int]) -> String {
        let names = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"]
        let set = Set(days)
        if set == [2, 3, 4, 5, 6] { return "Weekdays" }
        if set.count == 7 { return "Every day" }
        return days.sorted().compactMap { (1...7).contains($0) ? names[$0 - 1] : nil }.joined(separator: " ")
    }
}

public enum ProfileKind: String, Codable, Sendable {
    case work, custom
}

public struct Profile: Codable, Hashable, Identifiable, Sendable {
    public var id: UUID
    public var name: String
    public var kind: ProfileKind
    public var symbolName: String
    public var rules: [ProfileRule]
    /// Items are assigned when the combined weight of matching rules reaches this value.
    public var threshold: Double

    public init(id: UUID = UUID(), name: String, kind: ProfileKind = .custom, symbolName: String = "person.2.fill",
                rules: [ProfileRule] = [], threshold: Double = 0.8) {
        self.id = id
        self.name = name
        self.kind = kind
        self.symbolName = symbolName
        self.rules = rules
        self.threshold = threshold
    }

    /// The profile every install starts with.
    public static func defaultWork() -> Profile {
        Profile(name: "Work", kind: .work, symbolName: "briefcase.fill",
                rules: [.timeWindow(weekdays: [2, 3, 4, 5, 6], startMinute: 9 * 60, endMinute: 17 * 60),
                        .keyword("meeting"), .keyword("client"), .keyword("deadline"), .keyword("invoice"),
                        .calendarNamed("Work")],
                threshold: 0.8)
    }
}

public struct ProfileMatch: Hashable, Sendable {
    public var profile: Profile
    public var score: Double
    public var matchedRules: [ProfileRule]
}

public struct ProfileAssigner: Sendable {
    public var calendar: Calendar
    public init(calendar: Calendar = .current) { self.calendar = calendar }

    /// Profiles whose rules score at or above their threshold, best first.
    public func matches(for record: MemoryRecord, profiles: [Profile]) -> [ProfileMatch] {
        profiles.compactMap { profile -> ProfileMatch? in
            let hit = profile.rules.filter { matches($0, record) }
            let score = hit.reduce(0) { $0 + $1.weight }
            guard !hit.isEmpty, score >= profile.threshold else { return nil }
            return ProfileMatch(profile: profile, score: score, matchedRules: hit)
        }.sorted { $0.score > $1.score }
    }

    public func assign(_ record: MemoryRecord, profiles: [Profile]) -> [UUID] {
        matches(for: record, profiles: profiles).map { $0.profile.id }
    }

    func matches(_ rule: ProfileRule, _ r: MemoryRecord) -> Bool {
        switch rule {
        case .source(let s):
            return r.source == s
        case .accountContains(let a):
            return (r.accountLabel ?? "").localizedCaseInsensitiveContains(a)
        case .calendarNamed(let c):
            return (r.calendarName ?? "").localizedCaseInsensitiveCompare(c) == .orderedSame
        case .place(let center, let radius, _):
            guard let c = r.coordinate else { return false }
            return Geo.distance(c, center) <= radius
        case .timeWindow(let days, let start, let end):
            let comps = calendar.dateComponents([.weekday, .hour, .minute], from: r.date)
            guard let wd = comps.weekday, days.contains(wd) else { return false }
            let minute = (comps.hour ?? 0) * 60 + (comps.minute ?? 0)
            return minute >= start && minute < end
        case .person(let p):
            return r.people.contains { $0.localizedCaseInsensitiveCompare(p) == .orderedSame }
                || r.body.localizedCaseInsensitiveContains(p)
        case .keyword(let k):
            return r.searchableText.localizedCaseInsensitiveContains(k)
        }
    }

    /// Profile whose time windows cover `date`, used for "switch automatically".
    public func scheduledProfile(at date: Date, profiles: [Profile]) -> Profile? {
        let comps = calendar.dateComponents([.weekday, .hour, .minute], from: date)
        guard let wd = comps.weekday else { return nil }
        let minute = (comps.hour ?? 0) * 60 + (comps.minute ?? 0)
        return profiles.first { profile in
            profile.rules.contains { rule in
                if case .timeWindow(let days, let s, let e) = rule { return days.contains(wd) && minute >= s && minute < e }
                return false
            }
        }
    }
}

/// Which slice of memory a screen shows.
public enum ProfileFilter: Hashable, Sendable {
    case all
    case profile(UUID)
    case unassigned

    public func includes(_ profileIDs: [UUID]) -> Bool {
        switch self {
        case .all: return true
        case .profile(let id): return profileIDs.contains(id)
        case .unassigned: return profileIDs.isEmpty
        }
    }
}
