import SwiftUI
import RecallCore

/// Opened by the Action Button or the waveform button. Listens, shows what it understands, and saves to-dos for review.
struct ListeningView: View {
    @EnvironmentObject private var voice: VoiceSession
    @Environment(\.dismiss) private var dismiss
    @State private var savedMessage: String?

    var body: some View {
        NavigationStack {
            Group {
                switch voice.phase {
                case .failed(let message):
                    ContentUnavailableView {
                        Label("Can\u{2019}t Listen", systemImage: "mic.slash")
                    } description: {
                        Text(message)
                    } actions: {
                        Button("Open Settings") { openSystemSettings() }.buttonStyle(.bordered)
                    }
                default:
                    content
                }
            }
            .navigationTitle("Capture To-Dos")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar {
                ToolbarItem(placement: .cancellationAction) {
                    Button(voice.isActive ? "Cancel" : "Close") { voice.cancel(); dismiss() }
                }
            }
        }
        .task {
            if voice.phase == .idle { await voice.start() }
        }
        .onChange(of: voice.phase) { old, new in
            // A finished session closes itself and the to-do tab shows what was heard.
            if (old == .finishing || old == .listening), new == .idle {
                if voice.lastSavedCount > 0 {
                    savedMessage = "Saved \(voice.lastSavedCount) for review"
                    Task {
                        try? await Task.sleep(nanoseconds: 1_200_000_000)
                        Router.shared.selectedTab = .todos
                        dismiss()
                    }
                }
            }
        }
        .interactiveDismissDisabled(voice.phase == .listening)
    }

    private var content: some View {
        VStack(spacing: 20) {
            Spacer(minLength: 8)
            ZStack {
                Circle()
                    .fill(Color.accentColor.opacity(0.15))
                    .frame(width: 140, height: 140)
                    .scaleEffect(1 + CGFloat(voice.level) * 0.35)
                    .animation(.easeOut(duration: 0.1), value: voice.level)
                Circle().fill(Color.accentColor).frame(width: 96, height: 96)
                Image(systemName: icon)
                    .font(.system(size: 40, weight: .medium))
                    .foregroundStyle(.white)
                    .symbolEffect(.variableColor.iterative, isActive: voice.phase == .listening)
            }
            .accessibilityHidden(true)

            Text(statusText).font(.headline)

            if !voice.drafts.isEmpty {
                VStack(alignment: .leading, spacing: 8) {
                    ForEach(Array(voice.drafts.enumerated()), id: \.offset) { _, draft in
                        Label {
                            Text(draft.title)
                        } icon: {
                            Image(systemName: "checklist").foregroundStyle(.tint)
                        }
                        .frame(maxWidth: .infinity, alignment: .leading)
                    }
                }
                .padding()
                .background(.quaternary, in: RoundedRectangle(cornerRadius: 12, style: .continuous))
                .padding(.horizontal)
            } else if voice.phase == .listening {
                Text(voice.transcript.isEmpty ? "Say what you need to do, for example \u{201C}remind me to call Sam tomorrow at 5.\u{201D}" : voice.transcript)
                    .font(.callout)
                    .foregroundStyle(.secondary)
                    .multilineTextAlignment(.center)
                    .padding(.horizontal)
                    .lineLimit(4)
            }

            if let savedMessage {
                Label(savedMessage, systemImage: "checkmark.circle.fill").foregroundStyle(.green)
            }

            Spacer()

            if voice.phase == .idle && savedMessage == nil {
                Button { Task { await voice.start() } } label: { Text("Listen Again").frame(maxWidth: .infinity) }
                    .buttonStyle(.bordered)
                    .controlSize(.large)
                    .padding(.horizontal)
            }
            if voice.phase == .listening {
                Button { voice.stop() } label: {
                    Text("Done").frame(maxWidth: .infinity)
                }
                .buttonStyle(.borderedProminent)
                .controlSize(.large)
                .padding(.horizontal)
            }
            Text("Audio is understood on your iPhone and discarded. Only the to-dos are kept.")
                .font(.footnote).foregroundStyle(.secondary).multilineTextAlignment(.center).padding(.horizontal).padding(.bottom, 8)
        }
    }

    private var icon: String { voice.phase == .listening ? "waveform" : "mic.fill" }

    private var statusText: String {
        switch voice.phase {
        case .idle: return savedMessage == nil ? "No to-dos heard" : "Done"
        case .starting: return "Starting\u{2026}"
        case .listening: return voice.drafts.isEmpty ? "Listening\u{2026}" : "Found \(voice.drafts.count) to-do\(voice.drafts.count == 1 ? "" : "s")"
        case .finishing: return "Saving\u{2026}"
        case .failed: return ""
        }
    }
}
