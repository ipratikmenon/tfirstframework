import AppIntents
import Foundation
import SwiftData
import RecallCore

// MARK: Action Button

/// Assign this to the Action Button (Settings > Action Button > Shortcut) to start capturing to-dos by voice.
struct CaptureTodosIntent: AppIntent {
    static var title: LocalizedStringResource = "Capture To-Dos"
    static var description = IntentDescription("Listens for tasks and commitments in what you say and saves them as to-dos. Audio is never stored.")
    static var openAppWhenRun: Bool = true

    @MainActor
    func perform() async throws -> some IntentResult {
        UserDefaults.standard.set(true, forKey: Prefs.pendingListen)
        return .result()
    }
}

struct TotalRecallShortcuts: AppShortcutsProvider {
    static var appShortcuts: [AppShortcut] {
        AppShortcut(intent: CaptureTodosIntent(), phrases: ["Capture to-dos with \(.applicationName)", "Listen with \(.applicationName)"],
                    shortTitle: "Capture To-Dos", systemImageName: "waveform")
    }
}

// MARK: Focus filter (switch profile with a Focus)

struct ProfileEntity: AppEntity {
    static var typeDisplayRepresentation: TypeDisplayRepresentation = "Profile"
    static var defaultQuery = ProfileQuery()

    var id: String
    var name: String
    var displayRepresentation: DisplayRepresentation { DisplayRepresentation(title: "\(name)") }
}

struct ProfileQuery: EntityQuery {
    @MainActor
    func entities(for identifiers: [String]) async throws -> [ProfileEntity] {
        MemoryStore.shared.profiles().filter { identifiers.contains($0.id.uuidString) }
            .map { ProfileEntity(id: $0.id.uuidString, name: $0.name) }
    }

    @MainActor
    func suggestedEntities() async throws -> [ProfileEntity] {
        MemoryStore.shared.profiles().map { ProfileEntity(id: $0.id.uuidString, name: $0.name) }
    }
}

struct ProfileFocusFilter: SetFocusFilterIntent {
    static var title: LocalizedStringResource = "Set Profile"
    static var description: IntentDescription? = IntentDescription("Shows only one Total Recall profile while this Focus is on.")

    @Parameter(title: "Profile")
    var profile: ProfileEntity?

    var displayRepresentation: DisplayRepresentation {
        DisplayRepresentation(title: "\(profile?.name ?? "All memories")")
    }

    @MainActor
    func perform() async throws -> some IntentResult {
        if let profile, let id = UUID(uuidString: profile.id) {
            let locked = MemoryStore.shared.lockedProfileIDs()
            // A locked profile can't be opened without authentication, so Focus leaves the view alone.
            if !locked.contains(id) { ProfileState.shared.apply(.profile(id)) }
        } else {
            ProfileState.shared.apply(.all)
        }
        return .result()
    }
}
