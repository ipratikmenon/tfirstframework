import SwiftUI
import SwiftData
import RecallCore

struct TimelineScreen: View {
    @State private var day = Calendar.current.startOfDay(for: Date())
    @State private var showSettings = false

    var body: some View {
        NavigationStack {
            DayList(day: day)
                .navigationTitle(title)
                .navigationBarTitleDisplayMode(.large)
                .toolbar {
                    ToolbarItem(placement: .topBarLeading) { ProfileMenu() }
                    ToolbarItem(placement: .topBarTrailing) {
                        Button { showSettings = true } label: { Label("Settings", systemImage: "gearshape") }
                    }
                    ToolbarItemGroup(placement: .bottomBar) {
                        Button { shift(-1) } label: { Label("Previous day", systemImage: "chevron.left") }
                        Spacer()
                        DatePicker("Day", selection: $day, in: ...Date(), displayedComponents: .date)
                            .labelsHidden()
                        Spacer()
                        Button { shift(1) } label: { Label("Next day", systemImage: "chevron.right") }
                            .disabled(Calendar.current.isDateInToday(day))
                    }
                }
                .sheet(isPresented: $showSettings) { SettingsScreen() }
        }
    }

    private var title: String {
        if Calendar.current.isDateInToday(day) { return "Today" }
        if Calendar.current.isDateInYesterday(day) { return "Yesterday" }
        return day.formatted(.dateTime.weekday(.wide).month().day())
    }

    private func shift(_ days: Int) {
        if let d = Calendar.current.date(byAdding: .day, value: days, to: day), d <= Date() {
            day = Calendar.current.startOfDay(for: d)
        }
    }
}

private struct DayList: View {
    let day: Date
    @Query private var entries: [MemoryEntry]
    @Query private var profiles: [ProfileModel]
    @Query private var todos: [TodoItem]
    @EnvironmentObject private var state: ProfileState
    @State private var confirmDelete = false

    init(day: Date) {
        self.day = day
        let start = Calendar.current.startOfDay(for: day)
        let end = Calendar.current.date(byAdding: .day, value: 1, to: start) ?? start
        _entries = Query(filter: #Predicate<MemoryEntry> { $0.date >= start && $0.date < end },
                         sort: \MemoryEntry.date, order: .reverse)
        _todos = Query(filter: #Predicate<TodoItem> { $0.createdAt >= start && $0.createdAt < end })
    }

    private var visible: [MemoryEntry] {
        entries.filter { state.canShow($0.profileIDs, lockedIDs: profiles.lockedIDs) }
    }

    private var groups: [(title: String, items: [MemoryEntry])] {
        let cal = Calendar.current
        func part(_ e: MemoryEntry) -> Int {
            let h = cal.component(.hour, from: e.date)
            return h < 12 ? 0 : h < 18 ? 1 : 2
        }
        let names = ["Morning", "Afternoon", "Evening"]
        let grouped = Dictionary(grouping: visible, by: part)
        return [2, 1, 0].compactMap { k -> (title: String, items: [MemoryEntry])? in
            guard let items = grouped[k] else { return nil }
            return (title: names[k], items: items)
        }
    }

    var body: some View {
        List {
            if visible.isEmpty {
                ContentUnavailableView {
                    Label("Nothing yet", systemImage: "tray")
                } description: {
                    Text(state.filter == .all ? "Places, photos, events and mail you allow will appear here as the day goes on."
                                              : "Nothing in this profile for this day.")
                }
                .listRowBackground(Color.clear)
            } else {
                Section {
                    Text(DailyDigest.summary(for: visible.map(\.record),
                                             newTodos: todos.filter { state.canShow($0.profileIDs, lockedIDs: profiles.lockedIDs) }.count))
                        .font(.body)
                } header: { Text("Summary") }

                ForEach(groups, id: \.title) { group in
                    Section(group.title) {
                        ForEach(group.items) { entry in
                            NavigationLink { EntryDetailView(entry: entry) } label: { EntryRow(entry: entry) }
                        }
                        .onDelete { offsets in
                            for i in offsets { Persistence.shared.context.delete(group.items[i]) }
                            MemoryStore.shared.save()
                        }
                    }
                }

                Section {
                    Button("Delete This Day…", role: .destructive) { confirmDelete = true }
                }
            }
        }
        .refreshable { await CaptureCoordinator.shared.syncAll() }
        .confirmationDialog("Delete everything captured on this day?", isPresented: $confirmDelete, titleVisibility: .visible) {
            Button("Delete Day", role: .destructive) { MemoryStore.shared.deleteDay(day) }
        } message: {
            Text("This removes the items from Total Recall only. Your photos, events and mail are not touched.")
        }
    }
}

struct EntryRow: View {
    let entry: MemoryEntry
    var body: some View {
        HStack(spacing: 12) {
            SourceBadge(source: entry.source)
            VStack(alignment: .leading, spacing: 2) {
                Text(entry.title).lineLimit(2)
                Text(subtitle).font(.subheadline).foregroundStyle(.secondary).lineLimit(1)
            }
        }
        .accessibilityElement(children: .combine)
    }

    private var subtitle: String {
        var parts = [entry.date.formatted(date: .omitted, time: .shortened)]
        if let end = entry.endDate, end > entry.date { parts[0] += "\u{2013}" + end.formatted(date: .omitted, time: .shortened) }
        if entry.source != .location, let place = entry.placeName, !place.isEmpty { parts.append(place) }
        if entry.source == .mail, let from = entry.people.first { parts.append(from) }
        if entry.source != .mail, entry.source != .location, !entry.body.isEmpty { parts.append(entry.body.replacingOccurrences(of: "\n", with: " ")) }
        return parts.joined(separator: " \u{00B7} ")
    }
}

struct EntryDetailView: View {
    @Bindable var entry: MemoryEntry
    @Query(sort: [SortDescriptor(\ProfileModel.sortOrder), SortDescriptor(\ProfileModel.createdAt)]) private var profiles: [ProfileModel]
    @Environment(\.dismiss) private var dismiss
    @State private var confirmDelete = false

    var body: some View {
        Form {
            if let c = entry.coordinate {
                Section { CoordinateMap(coordinate: c, title: entry.title).listRowInsets(EdgeInsets()) }
            }
            Section {
                LabeledContent("Type") { Label(entry.source.displayName, systemImage: entry.source.symbolName) }
                LabeledContent("When") { Text(entry.date.formatted(date: .abbreviated, time: .shortened)) }
                if let place = entry.placeName, !place.isEmpty { LabeledContent("Place", value: place) }
                if let account = entry.accountLabel, !account.isEmpty { LabeledContent("Account", value: account) }
            }
            if !entry.body.isEmpty {
                Section("Details") { Text(entry.body).textSelection(.enabled) }
            }
            if !entry.people.isEmpty {
                Section("People") { ForEach(entry.people, id: \.self) { Label($0, systemImage: "person") } }
            }
            Section {
                ForEach(profiles) { p in
                    Toggle(isOn: membership(p)) { Label(p.name, systemImage: p.symbolName) }
                }
                if entry.profilesLocked {
                    Button("Assign Automatically") {
                        entry.profilesLocked = false
                        MemoryStore.shared.assignProfiles(entry)
                        MemoryStore.shared.save()
                    }
                }
            } header: { Text("Profiles") } footer: {
                Text(entry.profilesLocked ? "You chose these profiles. Automatic assignment is off for this item."
                                          : "Assigned automatically from your profile rules. Change a switch to take control.")
            }
            Section {
                Button("Delete from Total Recall", role: .destructive) { confirmDelete = true }
            }
        }
        .navigationTitle(entry.title)
        .navigationBarTitleDisplayMode(.inline)
        .confirmationDialog("Delete this item?", isPresented: $confirmDelete, titleVisibility: .visible) {
            Button("Delete", role: .destructive) {
                Persistence.shared.context.delete(entry)
                MemoryStore.shared.save()
                dismiss()
            }
        }
    }

    private func membership(_ p: ProfileModel) -> Binding<Bool> {
        Binding(get: { entry.profileIDs.contains(p.id) },
                set: { on in
                    var ids = entry.profileIDs
                    if on { if !ids.contains(p.id) { ids.append(p.id) } } else { ids.removeAll { $0 == p.id } }
                    entry.profileIDs = ids
                    entry.profilesLocked = true
                    MemoryStore.shared.save()
                })
    }
}
