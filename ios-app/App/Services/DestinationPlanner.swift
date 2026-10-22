import Foundation
import MapKit
import SwiftData
import Combine
import RecallCore

/// "I'm going to X": finds the route, then recalls everything relevant on the way or at the other end.
@MainActor
final class DestinationPlanner: NSObject, ObservableObject, MKLocalSearchCompleterDelegate {
    enum Transport: String, CaseIterable, Identifiable {
        case driving = "Drive", walking = "Walk", transit = "Transit"
        var id: String { rawValue }
        var mapType: MKDirectionsTransportType {
            switch self {
            case .driving: return .automobile
            case .walking: return .walking
            case .transit: return .transit
            }
        }
        var symbol: String {
            switch self {
            case .driving: return "car.fill"
            case .walking: return "figure.walk"
            case .transit: return "tram.fill"
            }
        }
    }

    struct Recall: Identifiable {
        let id: UUID
        let match: RouteMatch
        var title: String { match.candidate.title }
        var detail: String { match.reason.description }
        var symbol: String {
            switch match.candidate.kind {
            case .todo: return "checklist"
            case .event: return "calendar"
            case .memory: return "clock.arrow.circlepath"
            }
        }
    }

    @Published var query = "" { didSet { completer.queryFragment = query } }
    @Published private(set) var suggestions: [MKLocalSearchCompletion] = []
    @Published private(set) var destination: MKMapItem?
    @Published private(set) var route: MKRoute?
    @Published private(set) var recalls: [Recall] = []
    @Published private(set) var isPlanning = false
    @Published private(set) var message: String?
    @Published var transport: Transport = .driving { didSet { if destination != nil { Task { await plan() } } } }

    private let completer = MKLocalSearchCompleter()

    override init() {
        super.init()
        completer.delegate = self
        completer.resultTypes = [.address, .pointOfInterest, .query]
    }

    nonisolated func completerDidUpdateResults(_ completer: MKLocalSearchCompleter) {
        let results = completer.results
        Task { @MainActor in self.suggestions = results }
    }

    nonisolated func completer(_ completer: MKLocalSearchCompleter, didFailWithError error: Error) {}

    // MARK: choosing

    func choose(_ completion: MKLocalSearchCompletion) async {
        isPlanning = true
        defer { isPlanning = false }
        let request = MKLocalSearch.Request(completion: completion)
        guard let item = try? await MKLocalSearch(request: request).start().mapItems.first else {
            message = "Couldn't find that place."
            return
        }
        suggestions = []
        query = ""
        destination = item
        await plan()
    }

    func choose(coordinate: Coordinate, name: String) async {
        let item = MKMapItem(placemark: MKPlacemark(coordinate: CLLocationCoordinate2D(latitude: coordinate.latitude, longitude: coordinate.longitude)))
        item.name = name
        destination = item
        await plan()
    }

    func clear() {
        destination = nil; route = nil; recalls = []; message = nil; query = ""
    }

    // MARK: planning

    func plan() async {
        guard let destination else { return }
        isPlanning = true
        message = nil
        defer { isPlanning = false }

        let destCoord = Coordinate(latitude: destination.placemark.coordinate.latitude, longitude: destination.placemark.coordinate.longitude)

        // Route from the current position, when known.
        route = nil
        var path: [Coordinate] = [destCoord]
        if let here = LocationService.shared.lastKnown {
            let request = MKDirections.Request()
            request.source = MKMapItem(placemark: MKPlacemark(coordinate: here.coordinate))
            request.destination = destination
            request.transportType = transport.mapType
            if let response = try? await MKDirections(request: request).calculate(), let best = response.routes.first {
                route = best
                path = best.polyline.coordinatesArray.map { Coordinate(latitude: $0.latitude, longitude: $0.longitude) }
            } else {
                message = "No route found, showing what is relevant at the destination."
            }
        } else {
            message = "Turn on location to include your route. Showing what is relevant at the destination."
        }

        // Candidates: open to-dos, memories near the destination, upcoming events with a place.
        let store = MemoryStore.shared
        let locked = store.lockedProfileIDs()
        let filter = ProfileState.shared.filter
        var candidates: [RouteCandidate] = store.openTodos()
            .filter { filter.includes($0.profileIDs) && Set($0.profileIDs).isDisjoint(with: locked) }
            .map { RouteCandidate(id: $0.id, kind: .todo, title: $0.title, category: $0.category, coordinate: $0.coordinate, date: $0.due) }

        let since = Calendar.current.date(byAdding: .day, value: -365, to: Date()) ?? Date()
        let weekAhead = Calendar.current.date(byAdding: .day, value: 7, to: Date()) ?? Date()
        let memories = (try? Persistence.shared.context.fetch(
            FetchDescriptor<MemoryEntry>(predicate: #Predicate { $0.date >= since && $0.date <= weekAhead && $0.latitude != nil }))) ?? []
        for m in memories where filter.includes(m.profileIDs) && Set(m.profileIDs).isDisjoint(with: locked) {
            guard m.source != .location || m.placeName != nil else { continue }
            candidates.append(RouteCandidate(id: m.id, kind: m.source == .calendar ? .event : .memory,
                                             title: m.source == .location ? "Last time here: you were at \(m.placeName ?? m.title)" : m.title,
                                             coordinate: m.coordinate, date: m.date))
        }

        // Places for to-dos that name a kind of place ("pharmacy") but not a location.
        var nearby: [PlaceCategory: [NamedPlace]] = [:]
        let categories = Set(candidates.compactMap { $0.coordinate == nil ? $0.category : nil })
        if !categories.isEmpty {
            let region = Self.region(around: path)
            let center = path[path.count / 2]
            for c in categories { nearby[c] = await PlaceResolver.shared.search(c, in: region, center: center) }
        }

        let destCategory = PlaceCategory.detect(in: destination.name ?? "")
        let matches = RouteMatcher.match(candidates: candidates, route: path, destination: destCoord,
                                         destinationCategory: destCategory, nearbyPlaces: nearby)
        recalls = matches.prefix(12).map { Recall(id: $0.id, match: $0) }
        remember(destination: destination, coordinate: destCoord)
    }

    /// Posts the summary as a notification, for when the user is about to leave.
    func remindAtDeparture() async {
        guard !recalls.isEmpty else { return }
        let name = destination?.name ?? "your trip"
        let summary = recalls.prefix(3).map(\.title).joined(separator: ", ")
        await NotificationService.shared.postNow(title: "Before you go to \(name)", body: summary)
    }

    private func remember(destination: MKMapItem, coordinate: Coordinate) {
        let context = Persistence.shared.context
        let name = destination.name ?? "Destination"
        let all = (try? context.fetch(FetchDescriptor<DestinationRecord>())) ?? []
        if let existing = all.first(where: { Geo.distance($0.coordinate, coordinate) < 100 }) {
            existing.lastUsed = Date(); existing.useCount += 1
        } else {
            context.insert(DestinationRecord(name: name, coordinate: coordinate))
        }
        try? context.save()
    }

    static func region(around path: [Coordinate]) -> MKCoordinateRegion {
        let lats = path.map(\.latitude), lons = path.map(\.longitude)
        let minLat = lats.min() ?? 0, maxLat = lats.max() ?? 0, minLon = lons.min() ?? 0, maxLon = lons.max() ?? 0
        let center = CLLocationCoordinate2D(latitude: (minLat + maxLat) / 2, longitude: (minLon + maxLon) / 2)
        let span = MKCoordinateSpan(latitudeDelta: max(0.02, (maxLat - minLat) * 1.3), longitudeDelta: max(0.02, (maxLon - minLon) * 1.3))
        return MKCoordinateRegion(center: center, span: span)
    }
}

extension MKPolyline {
    var coordinatesArray: [CLLocationCoordinate2D] {
        var coords = [CLLocationCoordinate2D](repeating: CLLocationCoordinate2D(), count: pointCount)
        getCoordinates(&coords, range: NSRange(location: 0, length: pointCount))
        return coords
    }
}
