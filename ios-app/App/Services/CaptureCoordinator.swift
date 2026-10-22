import Foundation
import BackgroundTasks
import UIKit

/// Runs every enabled source's sync. Called on launch, on foreground, from pull-to-refresh
/// and from the system's background refresh.
@MainActor
final class CaptureCoordinator: ObservableObject {
    static let shared = CaptureCoordinator()
    static let refreshTaskID = "com.example.totalrecall.refresh"

    @Published private(set) var isSyncing = false
    @Published private(set) var lastSync: Date?

    func syncAll() async {
        guard !isSyncing else { return }
        isSyncing = true
        defer { isSyncing = false; lastSync = Date() }
        LocationService.shared.applySettings()
        await ShareInbox.ingest()
        CalendarService.shared.sync()
        await PhotoService.shared.sync()
        await HealthService.shared.sync()
        await GmailService.shared.sync()
        MemoryStore.shared.save()
        await AlertEngine.shared.replan()
    }

    // MARK: background refresh

    /// Call from application(_:didFinishLaunchingWithOptions:) — registration must happen before launch finishes.
    func registerBackgroundTask() {
        BGTaskScheduler.shared.register(forTaskWithIdentifier: Self.refreshTaskID, using: nil) { task in
            guard let task = task as? BGAppRefreshTask else { task.setTaskCompleted(success: false); return }
            let work = Task { @MainActor in
                await CaptureCoordinator.shared.syncAll()
                CaptureCoordinator.shared.scheduleBackgroundRefresh()
                task.setTaskCompleted(success: true)
            }
            task.expirationHandler = { work.cancel() }
        }
    }

    func scheduleBackgroundRefresh() {
        let request = BGAppRefreshTaskRequest(identifier: Self.refreshTaskID)
        request.earliestBeginDate = Date(timeIntervalSinceNow: 60 * 60)
        try? BGTaskScheduler.shared.submit(request)
    }
}
