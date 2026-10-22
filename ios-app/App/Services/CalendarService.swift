import Foundation
import EventKit
import RecallCore

@MainActor
final class CalendarService: ObservableObject {
    static let shared = CalendarService()
    private let store = EKEventStore()

    @Published private(set) var authorization = EKEventStore.authorizationStatus(for: .event)
    var isEnabled: Bool { UserDefaults.standard.bool(forKey: Prefs.enabled("calendar")) }
    var isAuthorized: Bool { authorization == .fullAccess }

    func requestAuthorization() async {
        _ = try? await store.requestFullAccessToEvents()
        authorization = EKEventStore.authorizationStatus(for: .event)
    }

    /// Past two weeks and next month, so upcoming events with locations can drive trip recall.
    func sync() {
        guard isEnabled, isAuthorized else { return }
        let cal = Calendar.current
        let start = cal.date(byAdding: .day, value: -14, to: Date()) ?? Date()
        let end = cal.date(byAdding: .day, value: 30, to: Date()) ?? Date()
        let events = store.events(matching: store.predicateForEvents(withStart: start, end: end, calendars: nil))
        for e in events {
            let attendees = (e.attendees ?? []).filter { !$0.isCurrentUser }.compactMap { $0.name }
            var body = [e.location, e.notes].compactMap { $0 }.filter { !$0.isEmpty }.joined(separator: "\n")
            if body.count > 600 { body = String(body.prefix(600)) }
            let geo = e.structuredLocation?.geoLocation?.coordinate
            let entry = MemoryEntry(date: e.startDate, endDate: e.endDate, source: .calendar,
                                    title: e.title ?? "Event", body: body,
                                    coordinate: geo.map { Coordinate(latitude: $0.latitude, longitude: $0.longitude) },
                                    placeName: e.structuredLocation?.title ?? e.location,
                                    people: attendees, accountLabel: e.calendar.source.title, calendarName: e.calendar.title,
                                    externalID: "event-\(e.eventIdentifier ?? e.title ?? "x")-\(Int(e.startDate.timeIntervalSince1970))")
            MemoryStore.shared.upsert(entry)
        }
        MemoryStore.shared.save()
    }
}
