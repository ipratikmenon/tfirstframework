import Foundation

public enum Geo {
    static let earthRadius = 6_371_000.0

    /// Great-circle distance in metres.
    public static func distance(_ a: Coordinate, _ b: Coordinate) -> Double {
        let lat1 = a.latitude * .pi / 180, lat2 = b.latitude * .pi / 180
        let dLat = lat2 - lat1
        let dLon = (b.longitude - a.longitude) * .pi / 180
        let h = sin(dLat / 2) * sin(dLat / 2) + cos(lat1) * cos(lat2) * sin(dLon / 2) * sin(dLon / 2)
        return 2 * earthRadius * asin(min(1, sqrt(h)))
    }

    /// Distance in metres from a point to the segment a-b (local flat projection, fine at city scale).
    public static func distance(from p: Coordinate, toSegment a: Coordinate, _ b: Coordinate) -> Double {
        let lat0 = p.latitude * .pi / 180
        func xy(_ c: Coordinate) -> (Double, Double) {
            ((c.longitude - p.longitude) * .pi / 180 * cos(lat0) * earthRadius,
             (c.latitude - p.latitude) * .pi / 180 * earthRadius)
        }
        let (ax, ay) = xy(a)
        let (bx, by) = xy(b)
        let dx = bx - ax, dy = by - ay
        let len2 = dx * dx + dy * dy
        if len2 == 0 { return (ax * ax + ay * ay).squareRoot() }
        let t = min(1, max(0, -(ax * dx + ay * dy) / len2))
        let px = ax + t * dx, py = ay + t * dy
        return (px * px + py * py).squareRoot()
    }

    /// Distance in metres from a point to a polyline.
    public static func distance(from p: Coordinate, toPath path: [Coordinate]) -> Double {
        guard let first = path.first else { return .infinity }
        if path.count == 1 { return distance(p, first) }
        var best = Double.infinity
        for i in 0..<(path.count - 1) {
            best = min(best, distance(from: p, toSegment: path[i], path[i + 1]))
        }
        return best
    }
}
