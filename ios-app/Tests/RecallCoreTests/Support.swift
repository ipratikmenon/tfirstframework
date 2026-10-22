import Foundation
@testable import RecallCore

/// Fixed clock so tests are deterministic: Wednesday 21 Oct 2026, 10:00 UTC.
enum Fixture {
    static var calendar: Calendar = {
        var c = Calendar(identifier: .gregorian)
        c.timeZone = TimeZone(identifier: "UTC")!
        c.firstWeekday = 2
        return c
    }()
    static let now: Date = calendar.date(from: DateComponents(year: 2026, month: 10, day: 21, hour: 10))!

    static func date(_ y: Int, _ m: Int, _ d: Int, _ h: Int = 0, _ min: Int = 0) -> Date {
        calendar.date(from: DateComponents(year: y, month: m, day: d, hour: h, minute: min))!
    }
}
