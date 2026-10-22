import SwiftUI
import SwiftData
import MapKit
import RecallCore

struct TripsScreen: View {
    @StateObject private var planner = DestinationPlanner()
    @Query(sort: \DestinationRecord.lastUsed, order: .reverse) private var recents: [DestinationRecord]
    @ObservedObject private var location = LocationService.shared
    @State private var position: MapCameraPosition = .automatic
    @State private var reminderSent = false

    var body: some View {
        NavigationStack {
            List {
                if !planner.suggestions.isEmpty && !planner.query.isEmpty {
                    Section {
                        ForEach(planner.suggestions.prefix(8), id: \.self) { s in
                            Button { Task { await planner.choose(s) } } label: {
                                VStack(alignment: .leading) {
                                    Text(s.title).foregroundStyle(.primary)
                                    if !s.subtitle.isEmpty { Text(s.subtitle).font(.subheadline).foregroundStyle(.secondary) }
                                }
                            }
                        }
                    }
                } else if let destination = planner.destination {
                    tripSections(destination)
                } else {
                    intro
                }
            }
            .navigationTitle("Trips")
            .searchable(text: $planner.query, prompt: "Where are you going?")
            .toolbar { ToolbarItem(placement: .topBarLeading) { ProfileMenu() } }
            .overlay { if planner.isPlanning { ProgressView().controlSize(.large) } }
        }
    }

    // MARK: sections

    @ViewBuilder
    private var intro: some View {
        Section {
            Label {
                Text("Tell Total Recall where you are going and it will bring back to-dos, places and plans that matter on the way.")
            } icon: { Image(systemName: "sparkles").foregroundStyle(.tint) }
            if location.authorization == .notDetermined || location.authorization == .denied || location.authorization == .restricted {
                Button("Allow Location for Routes") {
                    if location.authorization == .notDetermined { location.requestAuthorization() } else { openSystemSettings() }
                }
            }
        }
        if !recents.isEmpty {
            Section("Recent") {
                ForEach(recents.prefix(6)) { r in
                    Button { Task { await planner.choose(coordinate: r.coordinate, name: r.name) } } label: {
                        Label(r.name, systemImage: "clock.arrow.circlepath").foregroundStyle(.primary)
                    }
                }
            }
        }
        Section {
            Label("From Apple Maps or Google Maps, tap Share and choose Total Recall to get a reminder without opening the app.", systemImage: "square.and.arrow.up")
                .font(.footnote).foregroundStyle(.secondary)
        }
    }

    @ViewBuilder
    private func tripSections(_ destination: MKMapItem) -> some View {
        Section {
            Map(position: $position) {
                UserAnnotation()
                if let route = planner.route { MapPolyline(route.polyline).stroke(.blue, lineWidth: 5) }
                Marker(destination.name ?? "Destination", coordinate: destination.placemark.coordinate)
            }
            .mapControls { MapUserLocationButton() }
            .frame(height: 220)
            .listRowInsets(EdgeInsets())
            Picker("Travel", selection: $planner.transport) {
                ForEach(DestinationPlanner.Transport.allCases) { t in Label(t.rawValue, systemImage: t.symbol).tag(t) }
            }
            .pickerStyle(.segmented)
            if let route = planner.route {
                LabeledContent(destination.name ?? "Destination") {
                    Text("\(Duration.seconds(route.expectedTravelTime).formatted(.units(allowed: [.hours, .minutes], width: .abbreviated))) \u{00B7} \(Measurement(value: route.distance, unit: UnitLength.meters).formatted(.measurement(width: .abbreviated, usage: .road)))")
                }
            } else {
                Text(destination.name ?? "Destination").font(.headline)
            }
            if let message = planner.message { Text(message).font(.footnote).foregroundStyle(.secondary) }
        }

        Section {
            if planner.recalls.isEmpty {
                Label("Nothing to remember for this trip.", systemImage: "checkmark.circle").foregroundStyle(.secondary)
            } else {
                ForEach(planner.recalls) { r in
                    HStack(spacing: 12) {
                        IconTile(symbol: r.symbol, color: r.match.candidate.kind == .todo ? .orange : .blue)
                        VStack(alignment: .leading, spacing: 2) {
                            Text(r.title)
                            Text(r.detail).font(.subheadline).foregroundStyle(.secondary)
                        }
                    }
                    .accessibilityElement(children: .combine)
                }
            }
        } header: { Text("Remember on this trip") }

        Section {
            if !planner.recalls.isEmpty {
                Button(reminderSent ? "Reminder Sent" : "Remind Me Now") {
                    reminderSent = true
                    Task { _ = await NotificationService.shared.requestAuthorization(); await planner.remindAtDeparture() }
                }
                .disabled(reminderSent)
            }
            Button("Clear Trip", role: .destructive) { planner.clear(); reminderSent = false }
        }
    }
}
