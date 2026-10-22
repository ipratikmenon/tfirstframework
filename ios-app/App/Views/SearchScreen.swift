import SwiftUI
import SwiftData
import NaturalLanguage
import RecallCore

struct SearchScreen: View {
    @Query(sort: \MemoryEntry.date, order: .reverse) private var entries: [MemoryEntry]
    @Query private var profiles: [ProfileModel]
    @EnvironmentObject private var state: ProfileState
    @EnvironmentObject private var router: Router
    @State private var query = ""
    @State private var hits: [SearchHit] = []

    private static let examples = [
        "What did I do yesterday?",
        "Emails this week",
        "Photos last Saturday",
        "Where was I on Monday?",
        "Meetings last week",
    ]

    var body: some View {
        NavigationStack {
            List {
                if query.trimmingCharacters(in: .whitespaces).isEmpty {
                    Section {
                        ForEach(Self.examples, id: \.self) { example in
                            Button { query = example } label: {
                                Label(example, systemImage: "text.magnifyingglass").foregroundStyle(.primary)
                            }
                        }
                    } header: { Text("Try asking") } footer: {
                        Text("Search understands days (\u{201C}last Tuesday\u{201D}), types (\u{201C}emails\u{201D}, \u{201C}photos\u{201D}) and names. Everything is matched on your device.")
                    }
                } else if hits.isEmpty {
                    ContentUnavailableView.search(text: query).listRowBackground(Color.clear)
                } else {
                    Section(hits.count == 1 ? "1 match" : "\(hits.count) matches") {
                        ForEach(hits) { hit in
                            if let entry = lookup[hit.id] {
                                NavigationLink { EntryDetailView(entry: entry) } label: { HitRow(entry: entry, snippet: hit.snippet) }
                            }
                        }
                    }
                }
            }
            .navigationTitle("Ask")
            .searchable(text: $query, placement: .navigationBarDrawer(displayMode: .always), prompt: "Ask your memory")
            .toolbar { ToolbarItem(placement: .topBarLeading) { ProfileMenu() } }
            .task(id: SearchKey(query: query, count: entries.count, filter: state.filter)) { await runSearch() }
            .onAppear {
                if let prefill = router.searchPrefill { query = prefill; router.searchPrefill = nil }
            }
        }
    }

    private struct SearchKey: Hashable {
        let query: String
        let count: Int
        let filter: ProfileFilter
    }

    private var lookup: [UUID: MemoryEntry] { Dictionary(uniqueKeysWithValues: entries.map { ($0.id, $0) }) }

    private func runSearch() async {
        let q = query.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !q.isEmpty else { hits = []; return }
        try? await Task.sleep(nanoseconds: 250_000_000)       // debounce typing
        if Task.isCancelled { return }
        let locked = profiles.lockedIDs
        let visible = entries.filter { state.canShow($0.profileIDs, lockedIDs: locked) }.map(\.record)
        let embedding = NLEmbedding.sentenceEmbedding(for: .english)
        let semantic: (@Sendable (String, MemoryRecord) -> Double)? = embedding.map { e in
            { q, r in
                let d = e.distance(between: q, and: r.title)
                return d < 0.8 ? (0.8 - d) : 0
            }
        }
        // Profile filtering already happened above through canShow; search the visible set as a whole.
        hits = SearchEngine().search(query: q, in: visible, filter: .all, now: Date(), limit: 60, semantic: semantic)
    }
}

private struct HitRow: View {
    let entry: MemoryEntry
    let snippet: String
    var body: some View {
        HStack(spacing: 12) {
            SourceBadge(source: entry.source)
            VStack(alignment: .leading, spacing: 2) {
                Text(entry.title).lineLimit(2)
                if !snippet.isEmpty, snippet != entry.title {
                    Text(snippet).font(.subheadline).foregroundStyle(.secondary).lineLimit(2)
                }
                Text(entry.date.formatted(date: .abbreviated, time: .shortened)).font(.caption).foregroundStyle(.tertiary)
            }
        }
        .accessibilityElement(children: .combine)
    }
}
