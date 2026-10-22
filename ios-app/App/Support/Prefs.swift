import Foundation

/// UserDefaults keys used across the app (and the matching @AppStorage names).
enum Prefs {
    static let onboarded = "onboarded"
    static let appLock = "appLockEnabled"
    static let autoAddEnabled = "autoAddTodos"
    static let autoAddThreshold = "autoAddThreshold"
    static let alertsEnabled = "alertsEnabled"
    static let dailyAlertCap = "dailyAlertCap"
    static let autoSwitchProfiles = "autoSwitchProfiles"
    static let activeProfile = "activeProfileID"          // "all", "unassigned" or a UUID
    static let pendingListen = "pendingListen"
    static let glassesAlbumID = "glassesAlbumID"
    static let lastPhotoSync = "lastPhotoSync"
    static let alertCounterDay = "alertCounterDay"
    static let alertCounter = "alertCounter"

    static func enabled(_ source: String) -> String { "source.\(source).enabled" }

    static func register() {
        UserDefaults.standard.register(defaults: [
            autoAddEnabled: false,
            autoAddThreshold: 0.85,
            alertsEnabled: true,
            dailyAlertCap: 5,
            autoSwitchProfiles: false,
            activeProfile: "all",
        ])
    }
}
