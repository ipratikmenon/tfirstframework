import Foundation
import RecallCore

/// Turns stored to-dos plus learned routines into scheduled alerts.
@MainActor
final class AlertEngine {
    static let shared = AlertEngine()

    func replan(now: Date = Date()) async {
        guard UserDefaults.standard.bool(forKey: Prefs.alertsEnabled) else { return }
        let store = MemoryStore.shared
        let locked = store.lockedProfileIDs()
        let todos = store.openTodos().filter { Set($0.profileIDs).isDisjoint(with: locked) }.map(\.record)
        guard !todos.isEmpty else { return }

        // Routine model from the last 10 weeks of visits.
        let since = Calendar.current.date(byAdding: .day, value: -70, to: now) ?? now
        let visits = HistoryQueries.visits(since: since)
        let routine = visits.count >= 6 ? RoutineModel(visits: visits) : nil

        // Places for category-only to-dos, near where the user usually is.
        var known: [PlaceCategory: [NamedPlace]] = [:]
        if let here = LocationService.shared.lastKnown {
            let c = Coordinate(latitude: here.coordinate.latitude, longitude: here.coordinate.longitude)
            await PlaceResolver.shared.warm(near: c)
            for cat in PlaceCategory.allCases {
                if let p = PlaceResolver.shared.cachedNearest(cat, near: c) { known[cat] = [p] }
            }
        }

        let cap = max(1, UserDefaults.standard.integer(forKey: Prefs.dailyAlertCap))
        let planned = AlertPlanner.plan(todos: todos, routine: routine, knownPlaces: known, now: now,
                                        dailyCap: cap, alreadyAlertedToday: NotificationService.shared.alertsSentToday())
        await NotificationService.shared.schedule(planned)
        for alert in planned {
            if let todo = store.openTodos().first(where: { $0.id == alert.id }) { todo.lastAlerted = now }
        }
        store.save()
        GeofenceRefresher.refresh()
    }
}
