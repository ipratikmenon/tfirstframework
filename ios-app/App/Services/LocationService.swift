import Foundation
import CoreLocation
import MapKit
import SwiftData
import RecallCore

/// Visits and place triggers. Visit monitoring is the low-battery way to learn where
/// the user spends time: iOS wakes the app when a stay starts or ends.
@MainActor
final class LocationService: NSObject, ObservableObject, CLLocationManagerDelegate {
    static let shared = LocationService()

    @Published private(set) var authorization: CLAuthorizationStatus = .notDetermined
    private let manager = CLLocationManager()
    private let geocoder = CLGeocoder()
    private(set) var lastKnown: CLLocation?

    private override init() {
        super.init()
        authorization = manager.authorizationStatus
        manager.delegate = self
        manager.pausesLocationUpdatesAutomatically = true
        lastKnown = manager.location
    }

    var isEnabled: Bool { UserDefaults.standard.bool(forKey: Prefs.enabled("location")) }
    var hasAlways: Bool { authorization == .authorizedAlways }

    func requestAuthorization() {
        // iOS asks for While Using first; Always can be granted from the follow-up prompt or Settings.
        if authorization == .notDetermined { manager.requestWhenInUseAuthorization() }
        else if authorization == .authorizedWhenInUse { manager.requestAlwaysAuthorization() }
    }

    /// Starts or stops monitoring according to the user's setting. Call at every launch:
    /// when iOS relaunches the app for a visit this re-attaches the delegate.
    func applySettings() {
        if isEnabled && (authorization == .authorizedAlways || authorization == .authorizedWhenInUse) {
            manager.startMonitoringVisits()
            manager.startMonitoringSignificantLocationChanges()
        } else {
            manager.stopMonitoringVisits()
            manager.stopMonitoringSignificantLocationChanges()
            stopAllGeofences()
        }
    }

    func requestOneShot() { manager.requestLocation() }

    // MARK: delegate

    nonisolated func locationManagerDidChangeAuthorization(_ manager: CLLocationManager) {
        let status = manager.authorizationStatus
        Task { @MainActor in
            self.authorization = status
            self.applySettings()
        }
    }

    nonisolated func locationManager(_ manager: CLLocationManager, didUpdateLocations locations: [CLLocation]) {
        guard let last = locations.last else { return }
        Task { @MainActor in
            self.lastKnown = last
            GeofenceRefresher.refresh()
        }
    }

    nonisolated func locationManager(_ manager: CLLocationManager, didFailWithError error: Error) {}

    nonisolated func locationManager(_ manager: CLLocationManager, didVisit visit: CLVisit) {
        let arrival = visit.arrivalDate, departure = visit.departureDate, coord = visit.coordinate
        Task { @MainActor in await self.record(arrival: arrival, departure: departure, coordinate: coord) }
    }

    nonisolated func locationManager(_ manager: CLLocationManager, didEnterRegion region: CLRegion) {
        let id = region.identifier
        Task { @MainActor in await self.handleRegionEntry(identifier: id) }
    }

    // MARK: visits

    private func record(arrival: Date, departure: Date, coordinate: CLLocationCoordinate2D) async {
        guard isEnabled else { return }
        let start = arrival == .distantPast ? (departure == .distantFuture ? Date() : departure.addingTimeInterval(-1800)) : arrival
        let end: Date? = departure == .distantFuture ? nil : departure
        let name = await reverseGeocode(coordinate)
        let title = name ?? "Somewhere new"
        let entry = MemoryEntry(date: start, endDate: end, source: .location, title: title,
                                coordinate: Coordinate(latitude: coordinate.latitude, longitude: coordinate.longitude),
                                placeName: name,
                                externalID: String(format: "visit-%.4f-%.4f-%d", coordinate.latitude, coordinate.longitude, Int(start.timeIntervalSince1970 / 600)))
        MemoryStore.shared.upsert(entry)
        MemoryStore.shared.save()
        await AlertEngine.shared.replan()
    }

    func reverseGeocode(_ c: CLLocationCoordinate2D) async -> String? {
        let location = CLLocation(latitude: c.latitude, longitude: c.longitude)
        guard let placemark = try? await geocoder.reverseGeocodeLocation(location).first else { return nil }
        return placemark.name ?? placemark.thoroughfare ?? placemark.locality
    }

    // MARK: geofences (iOS allows 20 per app; we keep a few spare)

    static let maxRegions = 16

    func setGeofences(_ targets: [(id: UUID, coordinate: Coordinate)]) {
        guard CLLocationManager.isMonitoringAvailable(for: CLCircularRegion.self), hasAlways else { return }
        let wanted = Dictionary(uniqueKeysWithValues: targets.prefix(Self.maxRegions).map { ($0.id.uuidString, $0.coordinate) })
        for region in manager.monitoredRegions where wanted[region.identifier] == nil {
            manager.stopMonitoring(for: region)
        }
        let existing = Set(manager.monitoredRegions.map(\.identifier))
        for (id, c) in wanted where !existing.contains(id) {
            let region = CLCircularRegion(center: CLLocationCoordinate2D(latitude: c.latitude, longitude: c.longitude),
                                          radius: 150, identifier: id)
            region.notifyOnEntry = true
            region.notifyOnExit = false
            manager.startMonitoring(for: region)
        }
    }

    func stopAllGeofences() {
        for region in manager.monitoredRegions { manager.stopMonitoring(for: region) }
    }

    private func handleRegionEntry(identifier: String) async {
        guard let id = UUID(uuidString: identifier.components(separatedBy: "#").first ?? identifier) else { return }
        let context = Persistence.shared.context
        guard let todo = try? context.fetch(FetchDescriptor<TodoItem>(predicate: #Predicate { $0.id == id })).first,
              todo.status == .open else { return }
        if let last = todo.lastAlerted, Date().timeIntervalSince(last) < 4 * 3600 { return }
        todo.lastAlerted = Date()
        MemoryStore.shared.save()
        await NotificationService.shared.postNow(title: todo.title, body: todo.placeName.map { "You are near \($0)" } ?? "You are near a place for this", todoID: id)
    }
}

/// Chooses which open to-dos deserve a geofence right now.
@MainActor
enum GeofenceRefresher {
    static func refresh() {
        let service = LocationService.shared
        guard service.hasAlways, service.isEnabled, UserDefaults.standard.bool(forKey: Prefs.alertsEnabled),
              let here = service.lastKnown else { return }
        let hereCoord = Coordinate(latitude: here.coordinate.latitude, longitude: here.coordinate.longitude)
        let locked = MemoryStore.shared.lockedProfileIDs()
        var targets: [(UUID, Coordinate, Double)] = []
        for todo in MemoryStore.shared.openTodos() where Set(todo.profileIDs).isDisjoint(with: locked) {
            if let c = todo.coordinate {
                targets.append((todo.id, c, Geo.distance(c, hereCoord)))
            } else if let cat = todo.category, let place = PlaceResolver.shared.cachedNearest(cat, near: hereCoord) {
                targets.append((todo.id, place.coordinate, Geo.distance(place.coordinate, hereCoord)))
            }
        }
        // Nearest first: the limited slots go to what is most likely to matter soon.
        let nearest = targets.sorted { $0.2 < $1.2 }.map { ($0.0, $0.1) }
        service.setGeofences(nearest.map { (id: $0.0, coordinate: $0.1) })
    }
}

/// Finds places of a category near a coordinate (Apple Maps search) and remembers the answers.
@MainActor
final class PlaceResolver {
    static let shared = PlaceResolver()
    private var cache: [String: [NamedPlace]] = [:]

    private func key(_ c: PlaceCategory, _ near: Coordinate) -> String {
        String(format: "%@-%.2f-%.2f", c.rawValue, near.latitude, near.longitude)
    }

    func cachedNearest(_ category: PlaceCategory, near: Coordinate) -> NamedPlace? {
        cache[key(category, near)]?.min { Geo.distance($0.coordinate, near) < Geo.distance($1.coordinate, near) }
    }

    /// Places of `category` within a region, nearest first.
    func search(_ category: PlaceCategory, in region: MKCoordinateRegion, center: Coordinate) async -> [NamedPlace] {
        let k = key(category, center)
        if let hit = cache[k] { return hit }
        let request = MKLocalSearch.Request()
        request.naturalLanguageQuery = category.searchQuery
        request.region = region
        request.resultTypes = .pointOfInterest
        let places: [NamedPlace]
        if let response = try? await MKLocalSearch(request: request).start() {
            places = response.mapItems.compactMap { item in
                guard let name = item.name else { return nil }
                let c = item.placemark.coordinate
                return NamedPlace(name: name, coordinate: Coordinate(latitude: c.latitude, longitude: c.longitude))
            }
        } else {
            places = []
        }
        cache[k] = places
        return places
    }

    /// Warm the cache for every category that open to-dos need near `coordinate`.
    func warm(near coordinate: Coordinate) async {
        let categories = Set(MemoryStore.shared.openTodos().filter { $0.coordinate == nil }.compactMap(\.category))
        let region = MKCoordinateRegion(center: CLLocationCoordinate2D(latitude: coordinate.latitude, longitude: coordinate.longitude),
                                        latitudinalMeters: 6000, longitudinalMeters: 6000)
        for c in categories { _ = await search(c, in: region, center: coordinate) }
    }
}
