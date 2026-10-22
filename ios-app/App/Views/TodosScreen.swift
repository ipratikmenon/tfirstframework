import SwiftUI
import SwiftData
import RecallCore

struct TodosScreen: View {
    enum Segment: String, CaseIterable, Identifiable {
        case open = "Open", review = "To Review", done = "Done"
        var id: String { rawValue }
    }

    @Query(sort: \TodoItem.createdAt, order: .reverse) private var todos: [TodoItem]
    @Query private var profiles: [ProfileModel]
    @EnvironmentObject private var state: ProfileState
    @EnvironmentObject private var router: Router
    @State private var segment: Segment = .open
    @State private var editing: TodoItem?
    @State private var adding = false

    private var visible: [TodoItem] {
        let locked = profiles.lockedIDs
        return todos.filter { state.canShow($0.profileIDs, lockedIDs: locked) }
    }

    private var shown: [TodoItem] {
        switch segment {
        case .open:
            return visible.filter { $0.status == .open }.sorted { a, b in
                switch (a.due, b.due) {
                case let (x?, y?): return x < y
                case (_?, nil): return true
                case (nil, _?): return false
                default: return a.createdAt > b.createdAt
                }
            }
        case .review: return visible.filter { $0.status == .draft }
        case .done: return visible.filter { $0.status == .done }.sorted { ($0.completedAt ?? .distantPast) > ($1.completedAt ?? .distantPast) }
        }
    }

    private var draftCount: Int { visible.filter { $0.status == .draft }.count }

    var body: some View {
        NavigationStack {
            List {
                Section {
                    Picker("Show", selection: $segment) {
                        ForEach(Segment.allCases) { s in
                            Text(s == .review && draftCount > 0 ? "\(s.rawValue) (\(draftCount))" : s.rawValue).tag(s)
                        }
                    }
                    .pickerStyle(.segmented)
                }
                .listRowBackground(Color.clear)
                .listRowInsets(EdgeInsets())

                if shown.isEmpty {
                    emptyState.listRowBackground(Color.clear)
                } else {
                    Section {
                        ForEach(shown) { todo in
                            TodoRow(todo: todo, showConfidence: segment == .review)
                                .contentShape(Rectangle())
                                .onTapGesture { editing = todo }
                                .swipeActions(edge: .leading, allowsFullSwipe: true) { leadingActions(todo) }
                                .swipeActions(edge: .trailing, allowsFullSwipe: segment != .done) { trailingActions(todo) }
                        }
                    } footer: {
                        if segment == .review {
                            Text("Heard in your voice sessions. Swipe right to keep a to-do, left to dismiss it.")
                        }
                    }
                }
            }
            .navigationTitle("To-Dos")
            .toolbar {
                ToolbarItem(placement: .topBarLeading) { ProfileMenu() }
                ToolbarItemGroup(placement: .topBarTrailing) {
                    if segment == .review && draftCount > 1 {
                        Button("Keep All") { keepAll() }
                    }
                    Button { router.showListening = true } label: { Label("Capture by Voice", systemImage: "waveform") }
                    Button { adding = true } label: { Label("Add To-Do", systemImage: "plus") }
                }
            }
            .sheet(item: $editing) { todo in TodoEditor(todo: todo) }
            .sheet(isPresented: $adding) { TodoEditor(todo: nil) }
        }
    }

    @ViewBuilder
    private var emptyState: some View {
        switch segment {
        case .open:
            ContentUnavailableView {
                Label("No open to-dos", systemImage: "checklist")
            } description: {
                Text("Press the Action Button and say what you need to do, or tap the waveform.")
            } actions: {
                Button("Capture by Voice") { router.showListening = true }.buttonStyle(.bordered)
            }
        case .review:
            ContentUnavailableView("Nothing to review", systemImage: "checkmark.circle", description: Text("To-dos heard in voice sessions wait here for you to keep or dismiss."))
        case .done:
            ContentUnavailableView("Nothing done yet", systemImage: "circle.dashed")
        }
    }

    @ViewBuilder
    private func leadingActions(_ todo: TodoItem) -> some View {
        switch todo.status {
        case .draft:
            Button { keep(todo) } label: { Label("Keep", systemImage: "checkmark") }.tint(.green)
        case .open:
            Button { complete(todo) } label: { Label("Done", systemImage: "checkmark.circle") }.tint(.green)
        default:
            Button { reopen(todo) } label: { Label("Reopen", systemImage: "arrow.uturn.backward") }.tint(.blue)
        }
    }

    @ViewBuilder
    private func trailingActions(_ todo: TodoItem) -> some View {
        if todo.status == .draft {
            Button(role: .destructive) { dismiss(todo) } label: { Label("Dismiss", systemImage: "xmark") }
        } else {
            Button(role: .destructive) { remove(todo) } label: { Label("Delete", systemImage: "trash") }
        }
    }

    private func keep(_ todo: TodoItem) {
        todo.status = .open
        MemoryStore.shared.save()
        Task { await AlertEngine.shared.replan() }
    }

    private func keepAll() {
        for t in shown where t.status == .draft { t.status = .open }
        MemoryStore.shared.save()
        Task { await AlertEngine.shared.replan() }
    }

    private func complete(_ todo: TodoItem) {
        todo.status = .done
        todo.completedAt = Date()
        NotificationService.shared.cancel(todoID: todo.id)
        MemoryStore.shared.save()
    }

    private func reopen(_ todo: TodoItem) {
        todo.status = .open
        todo.completedAt = nil
        MemoryStore.shared.save()
    }

    private func dismiss(_ todo: TodoItem) {
        Persistence.shared.context.delete(todo)
        MemoryStore.shared.save()
    }

    private func remove(_ todo: TodoItem) {
        NotificationService.shared.cancel(todoID: todo.id)
        Persistence.shared.context.delete(todo)
        MemoryStore.shared.save()
    }
}

struct TodoRow: View {
    let todo: TodoItem
    var showConfidence = false

    var body: some View {
        HStack(alignment: .top, spacing: 12) {
            Image(systemName: todo.status == .done ? "checkmark.circle.fill" : "circle")
                .font(.title3)
                .foregroundStyle(todo.status == .done ? Color.green : Color.secondary)
                .accessibilityHidden(true)
            VStack(alignment: .leading, spacing: 3) {
                Text(todo.title).strikethrough(todo.status == .done).foregroundStyle(todo.status == .done ? .secondary : .primary)
                HStack(spacing: 8) {
                    if let due = todo.due {
                        Label(dueText(due), systemImage: "clock").foregroundStyle(due < Date() && todo.status == .open ? Color.red : Color.secondary)
                    }
                    if let cat = todo.category {
                        Label(todo.placeName ?? cat.displayName, systemImage: cat.symbolName)
                    } else if let place = todo.placeName {
                        Label(place, systemImage: "mappin")
                    }
                    if let person = todo.person { Label(person, systemImage: "person") }
                }
                .font(.caption)
                .foregroundStyle(.secondary)
                .labelStyle(.titleAndIcon)
                if showConfidence {
                    Text("Heard \(todo.createdAt.formatted(.relative(presentation: .named))) \u{00B7} \(confidenceText)")
                        .font(.caption2).foregroundStyle(.tertiary)
                }
            }
        }
        .padding(.vertical, 2)
        .accessibilityElement(children: .combine)
    }

    private var confidenceText: String {
        todo.confidence >= 0.85 ? "Very likely a to-do" : todo.confidence >= 0.65 ? "Likely a to-do" : "Maybe a to-do"
    }

    private func dueText(_ due: Date) -> String {
        todo.dueHasTime ? due.formatted(.dateTime.weekday(.abbreviated).hour().minute())
                        : due.formatted(.dateTime.weekday(.abbreviated).month(.abbreviated).day())
    }
}

struct TodoEditor: View {
    let todo: TodoItem?
    @Environment(\.dismiss) private var dismiss
    @Query(sort: [SortDescriptor(\ProfileModel.sortOrder), SortDescriptor(\ProfileModel.createdAt)]) private var profiles: [ProfileModel]

    @State private var title = ""
    @State private var hasDue = false
    @State private var due = Date().addingTimeInterval(3600)
    @State private var withTime = true
    @State private var category: PlaceCategory?
    @State private var placeName = ""
    @State private var person = ""
    @State private var chosenProfiles: Set<UUID> = []
    @State private var loaded = false

    var body: some View {
        NavigationStack {
            Form {
                Section { TextField("What needs doing?", text: $title, axis: .vertical) }
                Section("When") {
                    Toggle("Due", isOn: $hasDue.animation())
                    if hasDue {
                        Toggle("Time", isOn: $withTime)
                        DatePicker("Date", selection: $due, displayedComponents: withTime ? [.date, .hourAndMinute] : [.date])
                    }
                }
                Section {
                    Picker("Place type", selection: $category) {
                        Text("None").tag(PlaceCategory?.none)
                        ForEach(PlaceCategory.allCases, id: \.self) { c in
                            Label(c.displayName, systemImage: c.symbolName).tag(PlaceCategory?.some(c))
                        }
                    }
                    TextField("Place name (optional)", text: $placeName)
                    TextField("Person (optional)", text: $person)
                } header: { Text("Where and who") } footer: {
                    Text("A place type lets Total Recall remind you when your route passes one.")
                }
                Section {
                    ForEach(profiles) { p in
                        Toggle(isOn: Binding(get: { chosenProfiles.contains(p.id) },
                                             set: { on in if on { chosenProfiles.insert(p.id) } else { chosenProfiles.remove(p.id) } })) {
                            Label(p.name, systemImage: p.symbolName)
                        }
                    }
                } header: { Text("Profiles") } footer: {
                    Text(todo == nil ? "Leave all off to assign automatically." : "Changing these turns off automatic assignment.")
                }
            }
            .navigationTitle(todo == nil ? "New To-Do" : "Edit To-Do")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar {
                ToolbarItem(placement: .cancellationAction) { Button("Cancel") { dismiss() } }
                ToolbarItem(placement: .confirmationAction) {
                    Button(todo?.status == .draft ? "Keep" : "Save") { save() }
                        .disabled(title.trimmingCharacters(in: .whitespaces).isEmpty)
                }
            }
            .onAppear(perform: load)
        }
    }

    private func load() {
        guard !loaded else { return }
        loaded = true
        guard let todo else { return }
        title = todo.title
        if let d = todo.due { hasDue = true; due = d; withTime = todo.dueHasTime }
        category = todo.category
        placeName = todo.placeName ?? ""
        person = todo.person ?? ""
        chosenProfiles = Set(todo.profileIDs)
    }

    private func save() {
        let item = todo ?? TodoItem(title: title)
        if todo == nil { Persistence.shared.context.insert(item) }
        item.title = title.trimmingCharacters(in: .whitespaces)
        item.due = hasDue ? due : nil
        item.dueHasTime = hasDue && withTime
        item.category = category
        item.placeName = placeName.isEmpty ? nil : placeName
        item.person = person.isEmpty ? nil : person
        if item.status == .draft { item.status = .open }
        if !chosenProfiles.isEmpty || (todo != nil && Set(item.profileIDs) != chosenProfiles) {
            item.profileIDs = Array(chosenProfiles)
            item.profilesLocked = true
        } else {
            MemoryStore.shared.assignProfiles(item)
        }
        item.lastAlerted = nil
        MemoryStore.shared.save()
        Task {
            _ = await NotificationService.shared.requestAuthorization()
            await AlertEngine.shared.replan()
        }
        dismiss()
    }
}
