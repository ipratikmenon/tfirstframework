import SwiftUI
import SwiftData
import UIKit
import MapKit
import RecallCore

/// Rounded-square icon tile in the style of the Settings app.
struct IconTile: View {
    let symbol: String
    var color: Color = .accentColor
    var body: some View {
        Image(systemName: symbol)
            .font(.system(size: 15, weight: .semibold))
            .foregroundStyle(.white)
            .frame(width: 30, height: 30)
            .background(color.gradient, in: RoundedRectangle(cornerRadius: 7, style: .continuous))
            .accessibilityHidden(true)
    }
}

extension MemorySource {
    var tint: Color {
        switch self {
        case .location: return .blue
        case .photo: return .orange
        case .calendar: return .red
        case .health: return .pink
        case .mail: return .cyan
        case .share: return .indigo
        case .voice: return .purple
        case .glasses: return .teal
        case .manual: return .gray
        }
    }
}

struct SourceBadge: View {
    let source: MemorySource
    var body: some View { IconTile(symbol: source.symbolName, color: source.tint) }
}

extension Array where Element == ProfileModel {
    var lockedIDs: Set<UUID> { Set(filter(\.requiresAuth).map(\.id)) }
}

/// Wraps UIActivityViewController for exporting files.
struct ShareSheet: UIViewControllerRepresentable {
    let items: [Any]
    func makeUIViewController(context: Context) -> UIActivityViewController {
        UIActivityViewController(activityItems: items, applicationActivities: nil)
    }
    func updateUIViewController(_ controller: UIActivityViewController, context: Context) {}
}

struct CoordinateMap: View {
    let coordinate: Coordinate
    let title: String
    var body: some View {
        let center = CLLocationCoordinate2D(latitude: coordinate.latitude, longitude: coordinate.longitude)
        Map(initialPosition: .region(MKCoordinateRegion(center: center, latitudinalMeters: 800, longitudinalMeters: 800)),
            interactionModes: []) {
            Marker(title, coordinate: center)
        }
        .frame(height: 180)
        .accessibilityLabel("Map showing \(title)")
    }
}

/// Shared UI state that several tabs need to drive.
@MainActor
final class Router: ObservableObject {
    static let shared = Router()
    @Published var showListening = false
    @Published var selectedTab: AppTab = .timeline
    @Published var searchPrefill: String?
}

enum AppTab: Hashable { case timeline, search, todos, trips }

func openSystemSettings() {
    if let url = URL(string: UIApplication.openSettingsURLString) { UIApplication.shared.open(url) }
}
