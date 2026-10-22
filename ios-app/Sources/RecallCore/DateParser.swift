import Foundation

/// Small natural-language date and time parser for spoken to-dos
/// ("tomorrow at 5pm", "by Friday", "in 2 hours"). English only.
public struct DateParser: Sendable {
    public struct Result: Sendable, Equatable {
        public var date: Date
        public var hasTime: Bool
        /// The exact phrases that were understood, so callers can strip them from a title.
        public var phrases: [String]
    }

    public var calendar: Calendar
    public var defaultHour = 9

    public init(calendar: Calendar = .current) {
        self.calendar = calendar
    }

    static let weekdays = ["sunday": 1, "monday": 2, "tuesday": 3, "wednesday": 4,
                           "thursday": 5, "friday": 6, "saturday": 7]
    static let weekdayPattern = "(sunday|monday|tuesday|wednesday|thursday|friday|saturday)"

    public func parse(_ text: String, relativeTo now: Date = Date()) -> Result? {
        let l = text.lowercased()
        var phrases: [String] = []

        // 1. Relative offsets: "in 2 hours", "in 30 minutes", "in 3 days"
        if let m = firstMatch(#"\bin (\d+|an|a|one|two|three|four|five|six|seven|eight|nine|ten) (minute|hour|day|week)s?\b"#, in: l) {
            let n = Self.number(m[1]) ?? 1
            let unit = m[2]
            let comp: Calendar.Component = unit == "minute" ? .minute : unit == "hour" ? .hour : unit == "day" ? .day : .weekOfYear
            if let d = calendar.date(byAdding: comp, value: n, to: now) {
                let hasTime = unit == "minute" || unit == "hour"
                return Result(date: hasTime ? d : atHour(defaultHour, on: d), hasTime: hasTime, phrases: [m[0]])
            }
        }

        // 2. Day part
        var day: Date?
        let startToday = calendar.startOfDay(for: now)
        var impliedHour: Int?

        if let m = firstMatch(#"\bday after tomorrow\b"#, in: l) {
            day = calendar.date(byAdding: .day, value: 2, to: startToday); phrases.append(m[0])
        } else if let m = firstMatch(#"\btomorrow(?: (morning|afternoon|evening|night))?\b"#, in: l) {
            day = calendar.date(byAdding: .day, value: 1, to: startToday); phrases.append(m[0])
            impliedHour = Self.partOfDayHour(m[1])
        } else if let m = firstMatch(#"\btonight\b"#, in: l) {
            day = startToday; impliedHour = 20; phrases.append(m[0])
        } else if let m = firstMatch(#"\b(?:this|today)(?: (morning|afternoon|evening))?\b"#, in: l), m[0] != "this" {
            if m[0].hasPrefix("today") || !m[1].isEmpty {
                day = startToday; phrases.append(m[0]); impliedHour = Self.partOfDayHour(m[1])
            }
        }
        if day == nil, let m = firstMatch(#"\b(?:end of (?:the )?day|eod)\b"#, in: l) {
            day = startToday; impliedHour = 17; phrases.append(m[0])
        }
        if day == nil, let m = firstMatch(#"\bnext week\b"#, in: l) {
            let next = calendar.date(byAdding: .weekOfYear, value: 1, to: startToday) ?? startToday
            day = mondayOfWeek(containing: next); phrases.append(m[0])
        }
        if day == nil, let m = firstMatch(#"\b(?:this )?weekend\b"#, in: l) {
            day = next(weekday: 7, after: startToday, allowToday: true); phrases.append(m[0])
        }
        if day == nil, let m = firstMatch(#"\b(?:next |this |on |by |before |until )?"# + Self.weekdayPattern + #"(?: (morning|afternoon|evening))?\b"#, in: l) {
            let wd = Self.weekdays[m[1]] ?? 2
            let isNext = m[0].hasPrefix("next ")
            var d = next(weekday: wd, after: startToday, allowToday: false)
            if isNext, let plus = calendar.date(byAdding: .day, value: 7, to: d) {
                // "next friday" means the Friday of next week when this week's is only a day or two away.
                let gap = calendar.dateComponents([.day], from: startToday, to: d).day ?? 7
                if gap <= 2 { d = plus }
            }
            day = d; phrases.append(m[0]); impliedHour = Self.partOfDayHour(m[2])
        }
        if day == nil, let m = firstMatch(#"\bnext month\b"#, in: l) {
            day = calendar.date(byAdding: .month, value: 1, to: startToday); phrases.append(m[0])
        }

        // 3. Time part
        var hour: Int?
        var minute = 0
        if let m = firstMatch(#"\b(?:at |around |by )?(\d{1,2})(?::(\d{2}))?\s?(am|pm|a\.m\.|p\.m\.)(?!\w)"#, in: l) {
            var h = Int(m[1]) ?? 0
            minute = Int(m[2]) ?? 0
            let pm = m[3].hasPrefix("p")
            if h == 12 { h = pm ? 12 : 0 } else if pm { h += 12 }
            if h < 24 { hour = h; phrases.append(m[0]) }
        } else if let m = firstMatch(#"\bat (\d{1,2})(?::(\d{2}))?\b(?! ?(?:minutes|hours|days|people|items))"#, in: l) {
            // Bare "at 5": assume daytime (1-7 -> pm, 8-11 -> am, 12 -> noon).
            var h = Int(m[1]) ?? 0
            minute = Int(m[2]) ?? 0
            if (1...7).contains(h) { h += 12 }
            if h < 24 { hour = h; phrases.append(m[0]) }
        } else if let m = firstMatch(#"\b(noon|midnight)\b"#, in: l) {
            hour = m[1] == "noon" ? 12 : 0; phrases.append(m[0])
        }

        if day == nil && hour == nil { return nil }

        if let h = hour ?? impliedHour {
            let base = day ?? startToday
            var date = atHour(h, minute: minute, on: base)
            if day == nil, date <= now {   // time only and already passed -> tomorrow
                date = calendar.date(byAdding: .day, value: 1, to: date) ?? date
            }
            return Result(date: date, hasTime: true, phrases: phrases)
        }
        return Result(date: atHour(defaultHour, on: day ?? startToday), hasTime: false, phrases: phrases)
    }

    // MARK: helpers

    func atHour(_ hour: Int, minute: Int = 0, on day: Date) -> Date {
        calendar.date(bySettingHour: hour, minute: minute, second: 0, of: day) ?? day
    }

    func next(weekday: Int, after start: Date, allowToday: Bool) -> Date {
        for offset in (allowToday ? 0 : 1)...7 {
            if let d = calendar.date(byAdding: .day, value: offset, to: start),
               calendar.component(.weekday, from: d) == weekday { return d }
        }
        return start
    }

    func mondayOfWeek(containing date: Date) -> Date {
        let wd = calendar.component(.weekday, from: date)
        let back = (wd + 5) % 7   // Monday -> 0
        return calendar.date(byAdding: .day, value: -back, to: calendar.startOfDay(for: date)) ?? date
    }

    static func partOfDayHour(_ s: String) -> Int? {
        switch s {
        case "morning": return 9
        case "afternoon": return 15
        case "evening": return 18
        case "night": return 20
        default: return nil
        }
    }

    static func number(_ s: String) -> Int? {
        if let n = Int(s) { return n }
        let words = ["a": 1, "an": 1, "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
                     "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10]
        return words[s]
    }
}

/// Returns capture groups (index 0 is the whole match, missing groups are empty strings).
func firstMatch(_ pattern: String, in text: String) -> [String]? {
    guard let regex = try? NSRegularExpression(pattern: pattern, options: []) else { return nil }
    let ns = text as NSString
    guard let m = regex.firstMatch(in: text, options: [], range: NSRange(location: 0, length: ns.length)) else { return nil }
    return (0..<m.numberOfRanges).map { i in
        let r = m.range(at: i)
        return r.location == NSNotFound ? "" : ns.substring(with: r)
    }
}
