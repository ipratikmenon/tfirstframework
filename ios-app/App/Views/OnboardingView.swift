import SwiftUI
import RecallCore

struct OnboardingView: View {
    @AppStorage(Prefs.onboarded) private var onboarded = false
    @State private var step = 0
    @State private var working = false

    @AppStorage(Prefs.enabled("location")) private var location = true
    @AppStorage(Prefs.enabled("calendar")) private var calendar = true
    @AppStorage(Prefs.enabled("photo")) private var photos = true
    @AppStorage(Prefs.enabled("health")) private var health = false
    @AppStorage(Prefs.enabled("mail")) private var mail = false

    var body: some View {
        if step == 0 { welcome } else { sources }
    }

    // MARK: welcome

    private var welcome: some View {
        ScrollView {
            VStack(spacing: 28) {
                Image(systemName: "brain.head.profile")
                    .font(.system(size: 56))
                    .foregroundStyle(.tint)
                    .padding(.top, 48)
                    .accessibilityHidden(true)
                Text("Welcome to\nTotal Recall").font(.largeTitle.bold()).multilineTextAlignment(.center)
                VStack(alignment: .leading, spacing: 22) {
                    feature("lock.shield.fill", "Private by design", "Everything is stored and searched on your iPhone. Nothing is uploaded.")
                    feature("clock.fill", "A quiet record of your day", "Places, photos, events and mail you choose, gathered in the background.")
                    feature("waveform", "Say it once", "Press the Action Button and speak. Tasks you mention become to-dos. Audio is never saved.")
                    feature("map.fill", "Reminded at the right moment", "Heading somewhere? Total Recall brings back what matters on the way.")
                    feature("person.2.fill", "Work and personal, apart", "Profiles keep each part of your life in its own space.")
                }
                .padding(.horizontal, 28)
            }
        }
        .safeAreaInset(edge: .bottom) {
            Button { withAnimation { step = 1 } } label: { Text("Continue").frame(maxWidth: .infinity) }
                .buttonStyle(.borderedProminent).controlSize(.large)
                .padding()
                .background(.bar)
        }
    }

    private func feature(_ symbol: String, _ title: String, _ detail: String) -> some View {
        HStack(alignment: .top, spacing: 16) {
            Image(systemName: symbol).font(.title).foregroundStyle(.tint).frame(width: 40)
                .accessibilityHidden(true)
            VStack(alignment: .leading, spacing: 2) {
                Text(title).font(.headline)
                Text(detail).font(.subheadline).foregroundStyle(.secondary)
            }
        }
        .accessibilityElement(children: .combine)
    }

    // MARK: sources

    private var sources: some View {
        NavigationStack {
            Form {
                Section {
                    toggle("Places", "location.fill", .blue, $location, "Where you spend time, using low-power visit detection.")
                    toggle("Calendar", "calendar", .red, $calendar, "Events and who was there.")
                    toggle("Photos", "photo.fill", .orange, $photos, "When and where photos were taken, and text in screenshots.")
                    toggle("Health", "heart.fill", .pink, $health, "Daily steps and workouts.")
                    toggle("Mail", "envelope.fill", .cyan, $mail, "Gmail subjects and previews after you sign in. Apple Mail can\u{2019}t be read by any app.")
                } header: {
                    Text("Choose what to remember")
                } footer: {
                    Text("You can change this at any time in Settings. iOS asks for permission for each one.")
                }
            }
            .navigationTitle("Your Memory")
            .safeAreaInset(edge: .bottom) {
                Button { Task { await finish() } } label: {
                    Group {
                        if working { ProgressView() } else { Text("Turn On").frame(maxWidth: .infinity) }
                    }
                }
                .buttonStyle(.borderedProminent).controlSize(.large).disabled(working)
                .padding()
                .background(.bar)
            }
        }
    }

    private func toggle(_ title: String, _ symbol: String, _ color: Color, _ binding: Binding<Bool>, _ detail: String) -> some View {
        Toggle(isOn: binding) {
            HStack(spacing: 12) {
                IconTile(symbol: symbol, color: color)
                VStack(alignment: .leading, spacing: 2) {
                    Text(title)
                    Text(detail).font(.caption).foregroundStyle(.secondary)
                }
            }
        }
    }

    private func finish() async {
        working = true
        // @AppStorage defaults are not written until changed, so persist the choices explicitly.
        let defaults = UserDefaults.standard
        defaults.set(location, forKey: Prefs.enabled("location"))
        defaults.set(calendar, forKey: Prefs.enabled("calendar"))
        defaults.set(photos, forKey: Prefs.enabled("photo"))
        defaults.set(health, forKey: Prefs.enabled("health"))
        defaults.set(mail, forKey: Prefs.enabled("mail"))
        _ = await NotificationService.shared.requestAuthorization()
        if location { LocationService.shared.requestAuthorization() }
        if calendar { await CalendarService.shared.requestAuthorization() }
        if photos { await PhotoService.shared.requestAuthorization() }
        if health { await HealthService.shared.requestAuthorization() }
        working = false
        onboarded = true
        await CaptureCoordinator.shared.syncAll()
        if mail { await GmailService.shared.connect() }
    }
}
