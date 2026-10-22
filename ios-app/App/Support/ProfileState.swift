import Foundation
import SwiftUI
import LocalAuthentication
import RecallCore

/// Which profile the user is looking at, and which locked profiles are open right now.
@MainActor
final class ProfileState: ObservableObject {
    static let shared = ProfileState()

    @Published private(set) var filter: ProfileFilter
    /// Locked profiles the user has authenticated for in this session.
    @Published private(set) var unlocked: Set<UUID> = []

    private init() {
        filter = Self.decode(UserDefaults.standard.string(forKey: Prefs.activeProfile) ?? "all")
    }

    static func decode(_ s: String) -> ProfileFilter {
        switch s {
        case "all": return .all
        case "unassigned": return .unassigned
        default: return UUID(uuidString: s).map { .profile($0) } ?? .all
        }
    }

    static func encode(_ f: ProfileFilter) -> String {
        switch f {
        case .all: return "all"
        case .unassigned: return "unassigned"
        case .profile(let id): return id.uuidString
        }
    }

    /// Switch to a filter; locked profiles ask for Face ID / passcode first.
    func select(_ new: ProfileFilter, lockedIDs: Set<UUID>) async {
        if case .profile(let id) = new, lockedIDs.contains(id), !unlocked.contains(id) {
            guard await Self.authenticate(reason: "Open this profile") else { return }
            unlocked.insert(id)
        }
        apply(new)
    }

    /// Used by Focus filters and the schedule, which must not prompt.
    func apply(_ new: ProfileFilter) {
        filter = new
        UserDefaults.standard.set(Self.encode(new), forKey: Prefs.activeProfile)
    }

    /// Whether an item with these profiles may be shown now.
    func canShow(_ profileIDs: [UUID], lockedIDs: Set<UUID>) -> Bool {
        let hidden = Set(profileIDs).intersection(lockedIDs).subtracting(unlocked)
        return hidden.isEmpty && filter.includes(profileIDs)
    }

    func relock() {
        unlocked = []
        if case .profile = filter { apply(.all) }
    }

    static func authenticate(reason: String) async -> Bool {
        let context = LAContext()
        var error: NSError?
        guard context.canEvaluatePolicy(.deviceOwnerAuthentication, error: &error) else { return false }
        return (try? await context.evaluatePolicy(.deviceOwnerAuthentication, localizedReason: reason)) ?? false
    }
}

/// Face ID / passcode gate for the whole app.
@MainActor
final class AppLock: ObservableObject {
    static let shared = AppLock()
    @Published private(set) var isLocked: Bool

    private init() { isLocked = UserDefaults.standard.bool(forKey: Prefs.appLock) }

    var isEnabled: Bool { UserDefaults.standard.bool(forKey: Prefs.appLock) }

    func lockIfNeeded() { if isEnabled { isLocked = true } }

    func unlock() async {
        guard isLocked else { return }
        if await ProfileState.authenticate(reason: "Unlock Total Recall") { isLocked = false }
    }

    func setEnabled(_ on: Bool) async {
        if on {
            // Prove it works before turning it on, so the user can't lock themselves out.
            guard await ProfileState.authenticate(reason: "Turn on app lock") else { return }
        }
        UserDefaults.standard.set(on, forKey: Prefs.appLock)
        isLocked = false
        objectWillChange.send()
    }
}
