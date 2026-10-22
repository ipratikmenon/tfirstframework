import SwiftUI
import SwiftData
import RecallCore

/// Toolbar menu that switches which profile every screen shows.
struct ProfileMenu: View {
    @EnvironmentObject private var state: ProfileState
    @Query(sort: [SortDescriptor(\ProfileModel.sortOrder), SortDescriptor(\ProfileModel.createdAt)]) private var profiles: [ProfileModel]

    var body: some View {
        Menu {
            Picker("Profile", selection: binding) {
                Label("All", systemImage: "square.stack.3d.up").tag(ProfileFilter.all)
                ForEach(profiles) { p in
                    Label(p.name, systemImage: p.requiresAuth ? "lock.fill" : p.symbolName).tag(ProfileFilter.profile(p.id))
                }
                Label("Unassigned", systemImage: "questionmark.folder").tag(ProfileFilter.unassigned)
            }
        } label: {
            Label(currentName, systemImage: currentSymbol)
        }
        .accessibilityLabel("Profile: \(currentName)")
    }

    private var binding: Binding<ProfileFilter> {
        Binding(get: { state.filter },
                set: { new in Task { await state.select(new, lockedIDs: profiles.lockedIDs) } })
    }

    private var current: ProfileModel? {
        if case .profile(let id) = state.filter { return profiles.first { $0.id == id } }
        return nil
    }

    private var currentName: String {
        switch state.filter {
        case .all: return "All"
        case .unassigned: return "Unassigned"
        case .profile: return current?.name ?? "All"
        }
    }

    private var currentSymbol: String {
        switch state.filter {
        case .all: return "square.stack.3d.up"
        case .unassigned: return "questionmark.folder"
        case .profile: return current?.symbolName ?? "square.stack.3d.up"
        }
    }
}
