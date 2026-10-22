import Foundation

public struct Visit: Hashable, Sendable {
    public var coordinate: Coordinate
    public var arrival: Date
    public var departure: Date?
    public init(coordinate: Coordinate, arrival: Date, departure: Date? = nil) {
        self.coordinate = coordinate; self.arrival = arrival; self.departure = departure
    }
}

public struct RoutinePrediction: Hashable, Sendable {
    public var coordinate: Coordinate
    public var probability: Double
    public var expectedArrival: Date
    public init(coordinate: Coordinate, probability: Double, expectedArrival: Date) {
        self.coordinate = coordinate; self.probability = probability; self.expectedArrival = expectedArrival
    }
}

/// Learns where the user usually is by weekday and hour from past visits.
/// Everything stays on device; the model is rebuilt from the stored visit history.
public struct RoutineModel: Sendable {
    struct Cell: Hashable { var lat: Int; var lon: Int }
    struct Slot: Hashable { var weekday: Int; var hour: Int; var cell: Cell }

    let calendar: Calendar
    /// Distinct dates on which each slot's place was visited.
    private var slotDays: [Slot: Set<Date>] = [:]
    private var cellCenters: [Cell: Coordinate] = [:]
    private var observedDaysByWeekday: [Int: Int] = [:]

    static let cellSize = 0.0015   // roughly 150 m

    public init(visits: [Visit], calendar: Calendar = .current) {
        self.calendar = calendar
        guard let first = visits.map(\.arrival).min(), let last = visits.map(\.arrival).max() else { return }

        for v in visits {
            let cell = Self.cell(v.coordinate)
            cellCenters[cell] = Coordinate(latitude: (Double(cell.lat) + 0.5) * Self.cellSize,
                                           longitude: (Double(cell.lon) + 0.5) * Self.cellSize)
            let end = max(v.departure ?? v.arrival.addingTimeInterval(3600), v.arrival)
            var t = v.arrival
            while t <= end {
                let c = calendar.dateComponents([.weekday, .hour], from: t)
                let day = calendar.startOfDay(for: t)
                slotDays[Slot(weekday: c.weekday ?? 1, hour: c.hour ?? 0, cell: cell), default: []].insert(day)
                t = t.addingTimeInterval(3600)
            }
        }
        // Count how many of each weekday fall inside the observed range.
        var day = calendar.startOfDay(for: first)
        let lastDay = calendar.startOfDay(for: last)
        while day <= lastDay {
            observedDaysByWeekday[calendar.component(.weekday, from: day), default: 0] += 1
            guard let n = calendar.date(byAdding: .day, value: 1, to: day) else { break }
            day = n
        }
    }

    static func cell(_ c: Coordinate) -> Cell {
        Cell(lat: Int((c.latitude / cellSize).rounded(.down)), lon: Int((c.longitude / cellSize).rounded(.down)))
    }

    /// Places the user is likely to be during the next `hours`, most probable first.
    public func predictions(from now: Date, hours: Int = 6, minimumProbability: Double = 0.5,
                            minimumObservations: Int = 3) -> [RoutinePrediction] {
        var best: [Cell: RoutinePrediction] = [:]
        for h in 0..<hours {
            guard let t = calendar.date(byAdding: .hour, value: h, to: now) else { continue }
            let c = calendar.dateComponents([.weekday, .hour], from: t)
            guard let wd = c.weekday, let hour = c.hour, let observed = observedDaysByWeekday[wd], observed >= minimumObservations else { continue }
            for (slot, days) in slotDays where slot.weekday == wd && slot.hour == hour {
                let p = Double(days.count) / Double(observed)
                guard p >= minimumProbability, let center = cellCenters[slot.cell] else { continue }
                let arrival = calendar.date(bySettingHour: hour, minute: 0, second: 0, of: t) ?? t
                let prediction = RoutinePrediction(coordinate: center, probability: p, expectedArrival: max(arrival, now))
                if let existing = best[slot.cell], existing.probability >= p { continue }
                best[slot.cell] = prediction
            }
        }
        return best.values.sorted { $0.probability > $1.probability }
    }
}
