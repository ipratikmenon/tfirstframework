import SwiftUI
import SwiftData
import CoreLocation
import RecallCore

struct SettingsScreen: View {
    @Environment(\.dismiss) private var dismiss
    @EnvironmentObject private var lock: AppLock
    @ObservedObject private var location = LocationService.shared
    @ObservedObject private var photos = PhotoService.shared
    @ObservedObject private var calendar = CalendarService.shared
    @ObservedObject private var gmail = GmailService.shared

    @AppStorage(Prefs.enabled("location")) private var locationOn = false
    @AppStorage(Prefs.enabled("photo")) private var photoOn = false
    @AppStorage(Prefs.enabled("calendar")) private var calendarOn = false
    @AppStorage(Prefs.enabled("health")) private var healthOn = false
    @AppStorage(Prefs.enabled("mail")) private var mailOn = false
    @AppStorage(Prefs.alertsEnabled) private var alertsOn = true
    @AppStorage(Prefs.dailyAlertCap) private var alertCap = 5
    @AppStorage(Prefs.autoSwitchProfiles) private var autoSwitch = false

    @State private var exportURL: URL?
    @State private var confirmWipe = false
    @State private var sourceToDelete: MemorySource?

    var body: some View {
        NavigationStack {
            Form {
                captureSection
                voiceSection
                profilesSection
                alertsSection
                privacySection
                dataSection
                Section {
                    LabeledContent("Version", value: Bundle.main.object(forInfoDictionaryKey: "CFBundleShortVersionString") as? String ?? "1.0")
                }
            }
            .navigationTitle("Settings")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar { ToolbarItem(placement: .confirmationAction) { Button("Done") { dismiss() } } }
            .sheet(item: Binding(get: { exportURL.map(IdentifiedURL.init) }, set: { exportURL = $0?.url })) { ShareSheet(items: [$0.url]) }
            .confirmationDialog("Delete everything?", isPresented: $confirmWipe, titleVisibility: .visible) {
                Button("Delete Everything", role: .destructive) { MemoryStore.shared.deleteEverything() }
            } message: {
                Text("All memories, to-dos and saved trips are removed from this iPhone. This can\u{2019}t be undone.")
            }
            .confirmationDialog("Delete \(sourceToDelete?.displayName.lowercased() ?? "") items?", isPresented: Binding(get: { sourceToDelete != nil }, set: { if !$0 { sourceToDelete = nil } }),
                                titleVisibility: .visible) {
                Button("Delete", role: .destructive) { if let s = sourceToDelete { MemoryStore.shared.delete(source: s) }; sourceToDelete = nil }
            }
        }
    }

    // MARK: capture

    private var captureSection: some View {
        Section {
            sourceToggle("Places", "location.fill", .blue, $locationOn, status: locationStatus) { on in
                if on { location.requestAuthorization() }
                location.applySettings()
            }
            sourceToggle("Calendar", "calendar", .red, $calendarOn, status: calendar.isAuthorized ? nil : "Needs permission") { on in
                if on { Task { await calendar.requestAuthorization(); calendar.sync() } }
            }
            sourceToggle("Photos", "photo.fill", .orange, $photoOn, status: photos.isAuthorized ? nil : "Needs permission") { on in
                if on { Task { await photos.requestAuthorization(); await photos.sync() } }
            }
            sourceToggle("Health", "heart.fill", .pink, $healthOn, status: nil) { on in
                if on { Task { await HealthService.shared.requestAuthorization(); await HealthService.shared.sync() } }
            }
            NavigationLink { MailSettings() } label: {
                HStack(spacing: 12) {
                    IconTile(symbol: "envelope.fill", color: .cyan)
                    VStack(alignment: .leading) {
                        Text("Mail")
                        Text(gmail.isConnected ? (gmail.accountEmail ?? "Connected") : "Not connected").font(.caption).foregroundStyle(.secondary)
                    }
                }
            }
            NavigationLink { GlassesSettings() } label: {
                HStack(spacing: 12) {
                    IconTile(symbol: "eyeglasses", color: .teal)
                    Text("Ray-Ban Meta Glasses")
                }
            }
        } header: { Text("Capture") } footer: {
            Text("Each source is optional and can be paused here. Turning one off stops capture; use Data below to remove what was already saved.")
        }
    }

    private var locationStatus: String? {
        switch location.authorization {
        case .authorizedAlways: return nil
        case .authorizedWhenInUse: return "Set to Always for background capture"
        case .denied, .restricted: return "Off in Settings"
        default: return "Needs permission"
        }
    }

    private func sourceToggle(_ title: String, _ symbol: String, _ color: Color, _ binding: Binding<Bool>, status: String?,
                              onChange: @escaping (Bool) -> Void) -> some View {
        Toggle(isOn: binding) {
            HStack(spacing: 12) {
                IconTile(symbol: symbol, color: color)
                VStack(alignment: .leading, spacing: 1) {
                    Text(title)
                    if binding.wrappedValue, let status { Text(status).font(.caption).foregroundStyle(.orange) }
                }
            }
        }
        .onChange(of: binding.wrappedValue) { _, new in onChange(new) }
    }

    // MARK: voice

    private var voiceSection: some View {
        Section {
            NavigationLink { VoiceSettings() } label: {
                HStack(spacing: 12) {
                    IconTile(symbol: "waveform", color: .purple)
                    Text("Voice and Action Button")
                }
            }
        }
    }

    private var profilesSection: some View {
        Section {
            NavigationLink { ProfilesScreen() } label: {
                HStack(spacing: 12) {
                    IconTile(symbol: "person.2.fill", color: .green)
                    Text("Profiles")
                }
            }
            Toggle("Switch Automatically", isOn: $autoSwitch)
        } footer: {
            Text("Follows the hours of your Work profile. You can also change profile with a Focus: Settings > Focus > add the Total Recall filter.")
        }
    }

    private var alertsSection: some View {
        Section {
            Toggle("Alerts", isOn: $alertsOn)
            if alertsOn {
                Stepper("Up to \(alertCap) a day", value: $alertCap, in: 1...20)
            }
        } header: { Text("Reminders") } footer: {
            Text("Total Recall alerts you when a to-do is due, when a trip passes something relevant, and before places you usually go.")
        }
    }

    // MARK: privacy and data

    private var privacySection: some View {
        Section {
            Toggle("Require Face ID", isOn: Binding(get: { lock.isEnabled }, set: { on in Task { await lock.setEnabled(on) } }))
            NavigationLink("What Total Recall Can\u{2019}t See") { LimitsView() }
        } header: { Text("Privacy") } footer: {
            Text("Data stays on this iPhone, protected by the device passcode. Nothing is sent to a server.")
        }
    }

    private var dataSection: some View {
        Section {
            Button("Export All Data\u{2026}") { exportURL = try? MemoryStore.shared.exportJSON() }
            Menu("Delete by Type\u{2026}") {
                ForEach(MemorySource.allCases, id: \.self) { s in
                    Button(s.displayName, role: .destructive) { sourceToDelete = s }
                }
            }
            Button("Delete Everything\u{2026}", role: .destructive) { confirmWipe = true }
        } header: { Text("Your Data") }
    }
}

// MARK: - detail screens

struct MailSettings: View {
    @ObservedObject private var gmail = GmailService.shared
    @AppStorage(Prefs.enabled("mail")) private var mailOn = false

    var body: some View {
        Form {
            Section {
                if gmail.isConnected {
                    LabeledContent("Account", value: gmail.accountEmail ?? "Connected")
                    Toggle("Capture Mail", isOn: $mailOn)
                    Button("Sync Now") { Task { await gmail.sync() } }
                    Button("Disconnect", role: .destructive) { gmail.disconnect(); mailOn = false }
                } else {
                    Button("Connect Gmail") { Task { await gmail.connect(); if gmail.isConnected { mailOn = true; await gmail.sync() } } }
                        .disabled(!gmail.isConfigured)
                }
                if let error = gmail.lastError { Text(error).font(.footnote).foregroundStyle(.red) }
            } footer: {
                Text(gmail.isConfigured
                     ? "Total Recall reads the sender, subject and a short preview of recent messages, never full bodies or attachments. Disconnect at any time."
                     : "Gmail isn\u{2019}t set up in this build. Add a Google OAuth client id as GMAIL_CLIENT_ID in project.yml to enable it.")
            }
            Section {
                Label("Apple Mail has no interface for other apps, so it can\u{2019}t be read. Forward or share a message to Total Recall to save it.", systemImage: "info.circle")
                    .font(.footnote).foregroundStyle(.secondary)
            }
        }
        .navigationTitle("Mail")
        .navigationBarTitleDisplayMode(.inline)
    }
}

struct VoiceSettings: View {
    @AppStorage(Prefs.autoAddEnabled) private var autoAdd = false
    @AppStorage(Prefs.autoAddThreshold) private var threshold = 0.85

    var body: some View {
        Form {
            Section {
                Button { Router.shared.showListening = true } label: { Label("Try It Now", systemImage: "waveform") }
            } footer: {
                Text("Listening ends when you tap Done, after 30 seconds of silence, or after three minutes.")
            }
            Section {
                step(1, "Open Settings > Action Button.")
                step(2, "Swipe to Shortcut, then tap Choose a Shortcut.")
                step(3, "Pick Total Recall > Capture To-Dos.")
            } header: { Text("Set Up the Action Button") } footer: {
                Text("Requires iPhone 15 Pro or later. On other iPhones, add the Capture To-Dos shortcut to the Lock Screen or Control Center.")
            }
            Section {
                Toggle("Add Confident To-Dos Automatically", isOn: $autoAdd)
                if autoAdd {
                    VStack(alignment: .leading) {
                        Text("Confidence: \(Int(threshold * 100))%").font(.subheadline)
                        Slider(value: $threshold, in: 0.7...0.95, step: 0.05)
                    }
                }
            } header: { Text("Review") } footer: {
                Text(autoAdd ? "To-dos heard with at least this confidence skip the review list. Everything else waits for you."
                             : "Every to-do waits in To-Dos > To Review until you keep it.")
            }
            Section {
                Label("Speech is recognised on your iPhone. Audio is discarded as you speak and a transcript is never saved.", systemImage: "lock.shield")
                    .font(.footnote).foregroundStyle(.secondary)
            }
        }
        .navigationTitle("Voice")
        .navigationBarTitleDisplayMode(.inline)
    }

    private func step(_ n: Int, _ text: String) -> some View {
        HStack(alignment: .firstTextBaseline, spacing: 12) {
            Text("\(n)").font(.headline).foregroundStyle(.tint).frame(width: 20)
            Text(text)
        }
    }
}

struct GlassesSettings: View {
    @AppStorage(Prefs.glassesAlbumID) private var albumID = ""
    @ObservedObject private var photos = PhotoService.shared

    var body: some View {
        Form {
            Section {
                Picker("Glasses Album", selection: $albumID) {
                    Text("None").tag("")
                    ForEach(photos.albums()) { album in Text(album.title).tag(album.id) }
                }
            } header: { Text("Photos From Your Glasses") } footer: {
                Text("Pick the album where the Meta AI app saves photos from your Ray-Ban Meta glasses. They are added to your timeline as Glasses items, with text and objects found on device.")
            }
            if !photos.isAuthorized {
                Section {
                    Button("Allow Photos Access") { Task { await photos.requestAuthorization() } }
                }
            }
            Section {
                Label("Live capture straight from the glasses needs Meta\u{2019}s developer toolkit, which isn\u{2019}t available in this build.", systemImage: "info.circle")
                    .font(.footnote).foregroundStyle(.secondary)
            }
        }
        .navigationTitle("Ray-Ban Meta")
        .navigationBarTitleDisplayMode(.inline)
    }
}

struct LimitsView: View {
    var body: some View {
        List {
            Section {
                row("Apple Mail", "iOS doesn\u{2019}t let apps read it. Connect Gmail, or share a message to Total Recall.")
                row("Messages and iMessage", "No app can read them. Share a message to Total Recall to keep it.")
                row("Other apps\u{2019} notifications and screens", "Not available to any app.")
                row("Apple Maps destinations", "Not shared with other apps. Enter a destination in Trips, or share it from Maps.")
                row("Safari history", "Not available. Share a page to save it.")
            } header: { Text("What iOS keeps private") } footer: {
                Text("These limits protect you. Total Recall tells you plainly where its memory has gaps instead of pretending to be complete.")
            }
            Section("What Total Recall does") {
                row("Microphone", "Only when you start a session. Audio is understood on device and discarded.")
                row("Location", "Visits to places, using the system\u{2019}s low-power detection.")
                row("Photos", "Dates, places and text in screenshots. Your photos are never copied.")
            }
        }
        .navigationTitle("Privacy Limits")
        .navigationBarTitleDisplayMode(.inline)
    }

    private func row(_ title: String, _ detail: String) -> some View {
        VStack(alignment: .leading, spacing: 2) {
            Text(title)
            Text(detail).font(.subheadline).foregroundStyle(.secondary)
        }
        .accessibilityElement(children: .combine)
    }
}
