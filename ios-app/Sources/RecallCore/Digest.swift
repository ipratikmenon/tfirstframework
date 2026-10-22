import Foundation

/// Plain-language summary of a day, shown at the top of the timeline.
public enum DailyDigest {
    public static func summary(for records: [MemoryRecord], newTodos: Int = 0) -> String {
        guard !records.isEmpty else { return "Nothing captured yet." }
        func count(_ s: MemorySource) -> Int { records.filter { $0.source == s }.count }
        var parts: [String] = []
        let places = Set(records.filter { $0.source == .location }.compactMap { $0.placeName }).count
        if places > 0 { parts.append(places == 1 ? "visited 1 place" : "visited \(places) places") }
        let events = count(.calendar)
        if events > 0 { parts.append(events == 1 ? "had 1 event" : "had \(events) events") }
        let photos = count(.photo) + count(.glasses)
        if photos > 0 { parts.append(photos == 1 ? "took 1 photo" : "took \(photos) photos") }
        let mails = count(.mail)
        if mails > 0 { parts.append(mails == 1 ? "received 1 email" : "received \(mails) emails") }
        if newTodos > 0 { parts.append(newTodos == 1 ? "noted 1 to-do" : "noted \(newTodos) to-dos") }
        if let steps = records.first(where: { $0.source == .health && $0.title.contains("steps") }) {
            parts.append("recorded " + steps.title.lowercased())
        }
        if parts.isEmpty { return "\(records.count) items captured." }
        let sentence = parts.count == 1 ? parts[0] : parts.dropLast().joined(separator: ", ") + " and " + parts.last!
        return "You " + sentence + "."
    }
}
