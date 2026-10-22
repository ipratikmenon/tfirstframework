import Foundation

public struct SharedDestination: Equatable, Sendable {
    public var name: String?
    public var address: String?
    public var coordinate: Coordinate?
}

/// Reads destinations out of links shared from Apple Maps or Google Maps.
/// Short links (maps.app.goo.gl) cannot be resolved offline and return nil.
public enum MapsLinkParser {
    public static func parse(_ string: String) -> SharedDestination? {
        guard let comps = URLComponents(string: string.trimmingCharacters(in: .whitespacesAndNewlines)),
              let host = comps.host?.lowercased() else { return nil }
        let isApple = host == "maps.apple.com" || host == "maps.apple"
        let isGoogle = host.hasSuffix("google.com") && (comps.path.contains("/maps") || host.hasPrefix("maps."))
        guard isApple || isGoogle else { return nil }

        func item(_ names: String...) -> String? {
            for n in names {
                if let v = comps.queryItems?.first(where: { $0.name == n })?.value, !v.isEmpty {
                    return v.replacingOccurrences(of: "+", with: " ")
                }
            }
            return nil
        }

        var result = SharedDestination()
        if isApple {
            if let ll = item("ll", "coordinate", "sll"), let c = coordinate(from: ll) { result.coordinate = c }
            if let d = item("daddr"), let c = coordinate(from: d) { result.coordinate = result.coordinate ?? c }
            else if let d = item("daddr") { result.address = d }
            result.name = item("q", "name")
            result.address = result.address ?? item("address")
        } else {
            if let q = item("q", "query", "destination", "daddr") {
                if let c = coordinate(from: q) { result.coordinate = c } else { result.name = q }
            }
            // /maps/place/Blue+Bottle/@37.77,-122.41,17z
            if let range = comps.path.range(of: "/place/") {
                let after = comps.path[range.upperBound...]
                let name = after.split(separator: "/").first.map(String.init)
                result.name = result.name ?? name?.replacingOccurrences(of: "+", with: " ").removingPercentEncoding
            }
            if let at = comps.path.range(of: "/@") {
                let tail = comps.path[at.upperBound...].split(separator: ",")
                if tail.count >= 2, let lat = Double(tail[0]), let lon = Double(tail[1]) {
                    result.coordinate = result.coordinate ?? Coordinate(latitude: lat, longitude: lon)
                }
            }
        }
        return result.coordinate == nil && result.name == nil && result.address == nil ? nil : result
    }

    static func coordinate(from text: String) -> Coordinate? {
        let parts = text.split(separator: ",").map { $0.trimmingCharacters(in: .whitespaces) }
        guard parts.count == 2, let lat = Double(parts[0]), let lon = Double(parts[1]),
              abs(lat) <= 90, abs(lon) <= 180 else { return nil }
        return Coordinate(latitude: lat, longitude: lon)
    }
}
