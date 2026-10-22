import SwiftUI
import SwiftData
import RecallCore

struct RootView: View {
    @Environment(\.scenePhase) private var scenePhase
    @EnvironmentObject private var lock: AppLock
    @EnvironmentObject private var profileState: ProfileState
    @EnvironmentObject private var voice: VoiceSession
    @StateObject private var router = Router.shared
    @AppStorage(Prefs.onboarded) private var onboarded = false
    @AppStorage(Prefs.autoSwitchProfiles) private var autoSwitch = false
    @Query private var profiles: [ProfileModel]
    @Query(filter: #Predicate<TodoItem> { $0.statusRaw == "draft" }) private var drafts: [TodoItem]
    @State private var lastScheduled: UUID?

    var body: some View {
        Group {
            if onboarded { tabs } else { OnboardingView() }
        }
        .overlay {
            if onboarded && (lock.isLocked || (lock.isEnabled && scenePhase != .active)) { LockView() }
        }
        .sheet(isPresented: $router.showListening) {
            ListeningView()
                .presentationDetents([.medium, .large])
                .presentationDragIndicator(.visible)
        }
        .onChange(of: scenePhase) { _, phase in
            switch phase {
            case .active:
                handleActivation()
            case .background:
                lock.lockIfNeeded()
                profileState.relock()
                CaptureCoordinator.shared.scheduleBackgroundRefresh()
                MemoryStore.shared.save()
            default:
                break
            }
        }
        .task {
            handleActivation()
            while !Task.isCancelled {
                applySchedule()
                try? await Task.sleep(nanoseconds: 300 * 1_000_000_000)
            }
        }
    }

    private var tabs: some View {
        TabView(selection: $router.selectedTab) {
            TimelineScreen()
                .tabItem { Label("Timeline", systemImage: "clock") }
                .tag(AppTab.timeline)
            SearchScreen()
                .tabItem { Label("Ask", systemImage: "magnifyingglass") }
                .tag(AppTab.search)
            TodosScreen()
                .tabItem { Label("To-Dos", systemImage: "checklist") }
                .badge(drafts.count)
                .tag(AppTab.todos)
            TripsScreen()
                .tabItem { Label("Trips", systemImage: "map") }
                .tag(AppTab.trips)
        }
        .environmentObject(router)
    }

    private func handleActivation() {
        guard onboarded else { return }
        if lock.isLocked { Task { await lock.unlock() } }
        // The Action Button intent sets this flag and opens the app.
        if UserDefaults.standard.bool(forKey: Prefs.pendingListen) {
            UserDefaults.standard.set(false, forKey: Prefs.pendingListen)
            if !lock.isLocked { router.showListening = true } else { UserDefaults.standard.set(true, forKey: Prefs.pendingListen) }
        }
        Task { await CaptureCoordinator.shared.syncAll() }
    }

    /// "Switch automatically": follow the schedule of the profile whose time window is active.
    private func applySchedule() {
        guard autoSwitch, onboarded else { return }
        let scheduled = ProfileAssigner().scheduledProfile(at: Date(), profiles: profiles.map(\.profile))
        guard scheduled?.id != lastScheduled else { return }
        lastScheduled = scheduled?.id
        if let scheduled, !profiles.lockedIDs.contains(scheduled.id) { profileState.apply(.profile(scheduled.id)) }
        else if scheduled == nil { profileState.apply(.all) }
    }
}

struct LockView: View {
    @EnvironmentObject private var lock: AppLock

    var body: some View {
        ZStack {
            Rectangle().fill(.regularMaterial).ignoresSafeArea()
            VStack(spacing: 16) {
                Image(systemName: "lock.fill").font(.system(size: 44)).foregroundStyle(.secondary)
                Text("Total Recall is locked").font(.title2.weight(.semibold))
                Button("Unlock") { Task { await lock.unlock() } }
                    .buttonStyle(.borderedProminent)
                    .controlSize(.large)
            }
        }
        .transition(.opacity)
    }
}
