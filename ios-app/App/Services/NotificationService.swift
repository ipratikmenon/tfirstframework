import Foundation
import UserNotifications
import SwiftData
import RecallCore

/// Local notifications only. Nothing is sent to a server.
@MainActor
final class NotificationService: NSObject, UNUserNotificationCenterDelegate {
    static let shared = NotificationService()
    private let center = UNUserNotificationCenter.current()

    static let todoCategory = "TODO"
    static let doneAction = "DONE"
    static let snoozeAction = "SNOOZE"

    func configure() {
        center.delegate = self
        let done = UNNotificationAction(identifier: Self.doneAction, title: "Done", options: [])
        let snooze = UNNotificationAction(identifier: Self.snoozeAction, title: "Snooze 1 hour", options: [])
        center.setNotificationCategories([UNNotificationCategory(identifier: Self.todoCategory, actions: [done, snooze],
                                                                 intentIdentifiers: [], options: [])])
    }

    @discardableResult
    func requestAuthorization() async -> Bool {
        (try? await center.requestAuthorization(options: [.alert, .sound, .badge])) ?? false
    }

    func isAuthorized() async -> Bool {
        let s = await center.notificationSettings()
        return s.authorizationStatus == .authorized || s.authorizationStatus == .provisional
    }

    func cancelAll() { center.removeAllPendingNotificationRequests() }

    func cancel(todoID: UUID) {
        center.removePendingNotificationRequests(withIdentifiers: ["todo-\(todoID)-plan", "todo-\(todoID)-due"])
    }

    /// Replaces scheduled alerts with the planned set.
    func schedule(_ alerts: [PlannedAlert]) async {
        guard UserDefaults.standard.bool(forKey: Prefs.alertsEnabled), await isAuthorized() else { return }
        for alert in alerts {
            let content = UNMutableNotificationContent()
            content.title = alert.title
            content.body = alert.body
            content.sound = .default
            content.categoryIdentifier = Self.todoCategory
            content.userInfo = ["todoID": alert.id.uuidString]
            let interval = max(5, alert.fireDate.timeIntervalSinceNow)
            let trigger = UNTimeIntervalNotificationTrigger(timeInterval: interval, repeats: false)
            let id = "todo-\(alert.id)-plan"
            center.removePendingNotificationRequests(withIdentifiers: [id])
            try? await center.add(UNNotificationRequest(identifier: id, content: content, trigger: trigger))
            noteAlertScheduled()
        }
    }

    /// Fires right away (used for geofence entry and trip recall).
    func postNow(title: String, body: String, todoID: UUID? = nil) async {
        guard UserDefaults.standard.bool(forKey: Prefs.alertsEnabled), await isAuthorized() else { return }
        let content = UNMutableNotificationContent()
        content.title = title
        content.body = body
        content.sound = .default
        if let todoID {
            content.categoryIdentifier = Self.todoCategory
            content.userInfo = ["todoID": todoID.uuidString]
        }
        try? await center.add(UNNotificationRequest(identifier: UUID().uuidString, content: content,
                                                    trigger: UNTimeIntervalNotificationTrigger(timeInterval: 1, repeats: false)))
    }

    // MARK: daily cap bookkeeping

    func alertsSentToday() -> Int {
        let d = UserDefaults.standard
        let today = Calendar.current.startOfDay(for: Date()).timeIntervalSince1970
        return d.double(forKey: Prefs.alertCounterDay) == today ? d.integer(forKey: Prefs.alertCounter) : 0
    }

    private func noteAlertScheduled() {
        let d = UserDefaults.standard
        let today = Calendar.current.startOfDay(for: Date()).timeIntervalSince1970
        if d.double(forKey: Prefs.alertCounterDay) != today {
            d.set(today, forKey: Prefs.alertCounterDay)
            d.set(0, forKey: Prefs.alertCounter)
        }
        d.set(d.integer(forKey: Prefs.alertCounter) + 1, forKey: Prefs.alertCounter)
    }

    // MARK: delegate

    nonisolated func userNotificationCenter(_ center: UNUserNotificationCenter, willPresent notification: UNNotification) async
        -> UNNotificationPresentationOptions {
        [.banner, .list, .sound]
    }

    nonisolated func userNotificationCenter(_ center: UNUserNotificationCenter, didReceive response: UNNotificationResponse) async {
        guard let idString = response.notification.request.content.userInfo["todoID"] as? String,
              let id = UUID(uuidString: idString) else { return }
        let action = response.actionIdentifier
        await MainActor.run {
            let context = Persistence.shared.context
            guard let todo = try? context.fetch(FetchDescriptor<TodoItem>(predicate: #Predicate { $0.id == id })).first else { return }
            switch action {
            case Self.doneAction:
                todo.status = .done
                todo.completedAt = Date()
            case Self.snoozeAction:
                todo.lastAlerted = nil
                Task {
                    await NotificationService.shared.schedule([PlannedAlert(id: id, title: todo.title, body: "Snoozed reminder",
                                                                              fireDate: Date().addingTimeInterval(3600),
                                                                              reason: .due, score: 1)])
                }
            default:
                break
            }
            MemoryStore.shared.save()
        }
    }
}
