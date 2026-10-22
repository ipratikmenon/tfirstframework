import SwiftUI
import SwiftData
import UIKit
import RecallCore

@MainActor
final class AppDelegate: NSObject, UIApplicationDelegate {
    func application(_ application: UIApplication,
                     didFinishLaunchingWithOptions launchOptions: [UIApplication.LaunchOptionsKey: Any]? = nil) -> Bool {
        Prefs.register()
        // iOS can relaunch the app in the background for a visit or a region entry; everything the
        // system needs to deliver those events is attached here, before any UI exists.
        NotificationService.shared.configure()
        CaptureCoordinator.shared.registerBackgroundTask()
        CaptureCoordinator.shared.scheduleBackgroundRefresh()
        LocationService.shared.applySettings()
        return true
    }
}

@main
struct TotalRecallApp: App {
    @UIApplicationDelegateAdaptor(AppDelegate.self) private var delegate
    @StateObject private var profileState = ProfileState.shared
    @StateObject private var lock = AppLock.shared
    @StateObject private var voice = VoiceSession.shared

    var body: some Scene {
        WindowGroup {
            RootView()
                .environmentObject(profileState)
                .environmentObject(lock)
                .environmentObject(voice)
                .modelContainer(Persistence.shared.container)
        }
    }
}
