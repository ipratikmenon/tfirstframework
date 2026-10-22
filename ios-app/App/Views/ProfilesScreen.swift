import SwiftUI
import SwiftData
import RecallCore

struct ProfilesScreen: View {
    @Query(sort: [SortDescriptor(\ProfileModel.sortOrder), SortDescriptor(\ProfileModel.createdAt)]) private var profiles: [ProfileModel]
    @State private var creating = false

    var body: some View {
        List {
            Section {
                ForEach(profiles) { p in
                    NavigationLink { ProfileEditor(model: p) } label: {
                        HStack(spacing: 12) {
                            IconTile(symbol: p.requiresAuth ? "lock.fill" : p.symbolName, color: p.kindRaw == ProfileKind.work.rawValue ? .blue : .green)
                            VStack(alignment: .leading) {
                                Text(p.name)
                                Text(p.profile.rules.isEmpty ? "No rules yet" : "\(p.profile.rules.count) rule\(p.profile.rules.count == 1 ? "" : "s")")
                                    .font(.caption).foregroundStyle(.secondary)
                            }
                        }
                    }
                }
                .onDelete(perform: delete)
            } footer: {
                Text("A profile collects the memories and to-dos that belong to one part of your life. Rules decide what goes where, and you can always override them on an item.")
            }
            Section {
                Button { creating = true } label: { Label("Add Profile", systemImage: "plus.circle.fill") }
            }
        }
        .navigationTitle("Profiles")
        .navigationBarTitleDisplayMode(.inline)
        .sheet(isPresented: $creating) { NavigationStack { ProfileEditor(model: nil) } }
    }

    private func delete(at offsets: IndexSet) {
        for i in offsets {
            let model = profiles[i]
            guard model.kindRaw != ProfileKind.work.rawValue else { continue }   // Work is built in
            Persistence.shared.context.delete(model)
        }
        MemoryStore.shared.reassignAll()
    }
}

struct ProfileEditor: View {
    let model: ProfileModel?
    @Environment(\.dismiss) private var dismiss

    @State private var name = ""
    @State private var symbol = "person.2.fill"
    @State private var rules: [ProfileRule] = []
    @State private var threshold = 0.8
    @State private var requiresAuth = false
    @State private var addingRule = false
    @State private var loaded = false
    @State private var exportURL: URL?
    @State private var confirmDeleteData = false

    private static let symbols = ["briefcase.fill", "house.fill", "figure.2.and.child.holdinghands", "airplane", "graduationcap.fill",
                                  "heart.fill", "dumbbell.fill", "paintbrush.fill", "hammer.fill", "person.2.fill", "star.fill", "leaf.fill"]

    var body: some View {
        Form {
            Section {
                TextField("Name", text: $name)
                Picker("Icon", selection: $symbol) {
                    ForEach(Self.symbols, id: \.self) { Image(systemName: $0).tag($0) }
                }
            }

            Section {
                ForEach(Array(rules.enumerated()), id: \.offset) { _, rule in
                    Text(rule.summary)
                }
                .onDelete { rules.remove(atOffsets: $0) }
                Button { addingRule = true } label: { Label("Add Rule", systemImage: "plus.circle.fill") }
            } header: { Text("Rules") } footer: {
                Text("Accounts and calendars are strong signals. Places, people and words add up. Time of day only helps alongside another signal.")
            }

            Section {
                Picker("Matching", selection: $threshold) {
                    Text("Loose").tag(0.5)
                    Text("Balanced").tag(0.8)
                    Text("Strict").tag(1.2)
                }
                .pickerStyle(.segmented)
                Toggle("Require Face ID to Open", isOn: $requiresAuth)
            } header: { Text("Privacy") } footer: {
                Text("A locked profile is hidden from every list, search and trip until you unlock it, and is left out of notifications.")
            }

            if let model {
                Section("Data") {
                    Button("Export This Profile\u{2026}") {
                        exportURL = try? MemoryStore.shared.exportJSON(profile: model.profile)
                    }
                    Button("Delete This Profile\u{2019}s Data\u{2026}", role: .destructive) { confirmDeleteData = true }
                }
            }
        }
        .navigationTitle(model == nil ? "New Profile" : name)
        .navigationBarTitleDisplayMode(.inline)
        .toolbar {
            if model == nil {
                ToolbarItem(placement: .cancellationAction) { Button("Cancel") { dismiss() } }
            }
            ToolbarItem(placement: .confirmationAction) {
                Button(model == nil ? "Add" : "Save") { save() }.disabled(name.trimmingCharacters(in: .whitespaces).isEmpty)
            }
        }
        .onAppear(perform: load)
        .sheet(isPresented: $addingRule) { NavigationStack { AddRuleSheet { rules.append($0) } } }
        .sheet(item: Binding(get: { exportURL.map(IdentifiedURL.init) }, set: { exportURL = $0?.url })) { ShareSheet(items: [$0.url]) }
        .confirmationDialog("Delete everything in this profile?", isPresented: $confirmDeleteData, titleVisibility: .visible) {
            Button("Delete Data", role: .destructive) { if let model { MemoryStore.shared.deleteData(of: model.profile) } }
        } message: {
            Text("Items that also belong to another profile are kept.")
        }
    }

    private func load() {
        guard !loaded else { return }
        loaded = true
        guard let model else { return }
        name = model.name
        symbol = model.symbolName
        rules = model.profile.rules
        threshold = model.threshold
        requiresAuth = model.requiresAuth
    }

    private func save() {
        let trimmed = name.trimmingCharacters(in: .whitespaces)
        if let model {
            var p = model.profile
            p.name = trimmed; p.symbolName = symbol; p.rules = rules; p.threshold = threshold
            model.profile = p
            model.requiresAuth = requiresAuth
        } else {
            let p = Profile(name: trimmed, kind: .custom, symbolName: symbol, rules: rules, threshold: threshold)
            let context = Persistence.shared.context
            let count = (try? context.fetchCount(FetchDescriptor<ProfileModel>())) ?? 0
            context.insert(ProfileModel(profile: p, requiresAuth: requiresAuth, sortOrder: count))
        }
        MemoryStore.shared.reassignAll()
        dismiss()
    }
}

struct IdentifiedURL: Identifiable {
    let url: URL
    var id: URL { url }
}

struct AddRuleSheet: View {
    enum Kind: String, CaseIterable, Identifiable {
        case account = "Mail account", calendar = "Calendar", place = "Place", time = "Time of day", person = "Person", keyword = "Word", source = "Type of item"
        var id: String { rawValue }
    }

    let onAdd: (ProfileRule) -> Void
    @Environment(\.dismiss) private var dismiss
    @State private var kind: Kind = .keyword
    @State private var text = ""
    @State private var source: MemorySource = .mail
    @State private var radius = 200.0
    @State private var days: Set<Int> = [2, 3, 4, 5, 6]
    @State private var start = Calendar.current.date(bySettingHour: 9, minute: 0, second: 0, of: Date()) ?? Date()
    @State private var end = Calendar.current.date(bySettingHour: 17, minute: 0, second: 0, of: Date()) ?? Date()

    var body: some View {
        Form {
            Picker("Kind of Rule", selection: $kind) {
                ForEach(Kind.allCases) { Text($0.rawValue).tag($0) }
            }
            switch kind {
            case .account:
                Section {
                    TextField("For example @company.com", text: $text).textInputAutocapitalization(.never).autocorrectionDisabled()
                } footer: {
                    Text("Mail and calendar items from an account whose name contains this text.")
                }
            case .calendar:
                Section { TextField("Calendar name", text: $text) } footer: { Text("Events from the calendar with this name.") }
            case .place:
                Section {
                    TextField("Name, for example Office", text: $text)
                    Stepper("Within \(Int(radius)) m", value: $radius, in: 50...2000, step: 50)
                } footer: {
                    Text(LocationService.shared.lastKnown == nil ? "Needs your current location. Allow location access first."
                                                                  : "Uses where you are right now as the centre.")
                }
            case .time:
                Section("Days") {
                    ForEach(1...7, id: \.self) { d in
                        Toggle(Calendar.current.weekdaySymbols[d - 1], isOn: Binding(get: { days.contains(d) },
                                                                                       set: { on in if on { days.insert(d) } else { days.remove(d) } }))
                    }
                }
                Section("Hours") {
                    DatePicker("From", selection: $start, displayedComponents: .hourAndMinute)
                    DatePicker("To", selection: $end, displayedComponents: .hourAndMinute)
                }
            case .person:
                Section { TextField("Name", text: $text) } footer: { Text("Events, mail and to-dos that involve this person.") }
            case .keyword:
                Section { TextField("For example invoice", text: $text).textInputAutocapitalization(.never) } footer: { Text("Anything mentioning this word.") }
            case .source:
                Picker("Type", selection: $source) {
                    ForEach(MemorySource.allCases, id: \.self) { Label($0.displayName, systemImage: $0.symbolName).tag($0) }
                }
            }
        }
        .navigationTitle("Add Rule")
        .navigationBarTitleDisplayMode(.inline)
        .toolbar {
            ToolbarItem(placement: .cancellationAction) { Button("Cancel") { dismiss() } }
            ToolbarItem(placement: .confirmationAction) { Button("Add") { add() }.disabled(!canAdd) }
        }
    }

    private var canAdd: Bool {
        switch kind {
        case .source: return true
        case .time: return !days.isEmpty
        case .place: return !text.trimmingCharacters(in: .whitespaces).isEmpty && LocationService.shared.lastKnown != nil
        default: return !text.trimmingCharacters(in: .whitespaces).isEmpty
        }
    }

    private func add() {
        let t = text.trimmingCharacters(in: .whitespaces)
        let rule: ProfileRule
        switch kind {
        case .account: rule = .accountContains(t)
        case .calendar: rule = .calendarNamed(t)
        case .place:
            guard let here = LocationService.shared.lastKnown else { return }
            rule = .place(center: Coordinate(latitude: here.coordinate.latitude, longitude: here.coordinate.longitude), radiusMeters: radius, name: t)
        case .time:
            let cal = Calendar.current
            let s = cal.dateComponents([.hour, .minute], from: start), e = cal.dateComponents([.hour, .minute], from: end)
            rule = .timeWindow(weekdays: days.sorted(), startMinute: (s.hour ?? 9) * 60 + (s.minute ?? 0), endMinute: (e.hour ?? 17) * 60 + (e.minute ?? 0))
        case .person: rule = .person(t)
        case .keyword: rule = .keyword(t)
        case .source: rule = .source(source)
        }
        onAdd(rule)
        dismiss()
    }
}
