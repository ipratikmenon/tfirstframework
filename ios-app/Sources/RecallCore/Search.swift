import Foundation

public struct SearchHit: Hashable, Identifiable, Sendable {
    public var id: UUID
    public var score: Double
    public var snippet: String
}

/// Natural-language time expressions in a query: "last tuesday", "yesterday", "3 days ago".
public struct QueryTimeParser: Sendable {
    public var calendar: Calendar
    public init(calendar: Calendar = .current) { self.calendar = calendar }

    /// Returns the interval and the query with the time words removed.
    public func parse(_ query: String, now: Date) -> (interval: DateInterval?, remainder: String) {
        var q = query.lowercased()
        let today = calendar.startOfDay(for: now)
        func day(_ offset: Int) -> Date { calendar.date(byAdding: .day, value: offset, to: today) ?? today }
        func span(_ from: Date, _ to: Date) -> DateInterval { DateInterval(start: from, end: max(to, from.addingTimeInterval(1))) }
        func strip(_ phrase: String) { q = q.replacingOccurrences(of: phrase, with: " ") }

        if q.contains("yesterday") { strip("yesterday"); return (span(day(-1), today), q) }
        if q.contains("last night") { strip("last night"); return (span(day(-1).addingTimeInterval(18 * 3600), day(0).addingTimeInterval(6 * 3600)), q) }
        if q.contains("this morning") { strip("this morning"); return (span(today, today.addingTimeInterval(12 * 3600)), q) }
        if q.contains("earlier today") { strip("earlier today"); return (span(today, now), q) }
        if q.contains("today") { strip("today"); return (span(today, day(1)), q) }
        if q.contains("last week") {
            strip("last week")
            let startThis = startOfWeek(now)
            return (span(calendar.date(byAdding: .day, value: -7, to: startThis) ?? startThis, startThis), q)
        }
        if q.contains("this week") { strip("this week"); return (span(startOfWeek(now), day(1)), q) }
        if q.contains("last month") {
            strip("last month")
            let comps = calendar.dateComponents([.year, .month], from: now)
            let startThis = calendar.date(from: comps) ?? today
            return (span(calendar.date(byAdding: .month, value: -1, to: startThis) ?? startThis, startThis), q)
        }
        if q.contains("this month") {
            strip("this month")
            return (span(calendar.date(from: calendar.dateComponents([.year, .month], from: now)) ?? today, day(1)), q)
        }
        if let m = firstMatch(#"\b(\d+|a|one|two|three|four|five|six|seven) (day|week|month)s? ago\b"#, in: q) {
            let n = DateParser.number(m[1]) ?? 1
            let unit: Calendar.Component = m[2] == "day" ? .day : m[2] == "week" ? .weekOfYear : .month
            strip(m[0])
            let target = calendar.date(byAdding: unit, value: -n, to: today) ?? today
            let pad: Double = m[2] == "day" ? 1 : m[2] == "week" ? 3.5 : 15
            return (span(target.addingTimeInterval(-pad * 86400 / 2), target.addingTimeInterval(pad * 86400)), q)
        }
        if let m = firstMatch(#"\b(?:last|on|this past) "# + DateParser.weekdayPattern + #"\b"#, in: q) {
            let wd = DateParser.weekdays[m[1]] ?? 2
            var d = day(-1)
            for _ in 0..<7 where calendar.component(.weekday, from: d) != wd { d = calendar.date(byAdding: .day, value: -1, to: d) ?? d }
            strip(m[0])
            return (span(d, calendar.date(byAdding: .day, value: 1, to: d) ?? d), q)
        }
        let months = ["january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december"]
        for (i, name) in months.enumerated() {
            if let m = firstMatch(#"\bin "# + name + #"(?: (\d{4}))?\b"#, in: q) {
                var comps = calendar.dateComponents([.year], from: now)
                if let y = Int(m[1]) { comps.year = y }
                comps.month = i + 1; comps.day = 1
                guard let start = calendar.date(from: comps), let end = calendar.date(byAdding: .month, value: 1, to: start) else { break }
                strip(m[0])
                return (span(start, end), q)
            }
        }
        return (nil, q)
    }

    func startOfWeek(_ date: Date) -> Date {
        calendar.dateInterval(of: .weekOfYear, for: date)?.start ?? calendar.startOfDay(for: date)
    }
}

/// Ranked keyword search with light stemming, time filters and source hints (BM25).
/// An optional semantic scorer (for example Apple's NLEmbedding) can be blended in.
public struct SearchEngine: Sendable {
    public var calendar: Calendar
    public init(calendar: Calendar = .current) { self.calendar = calendar }

    static let stopWords: Set<String> = [
        "the", "a", "an", "and", "or", "of", "to", "in", "on", "at", "for", "with", "about", "was", "were", "is", "are",
        "did", "do", "does", "i", "me", "my", "we", "you", "what", "where", "when", "who", "which", "how", "that", "this",
        "it", "from", "by", "have", "had", "has", "say", "said", "tell", "find", "show", "any", "anything", "something"
    ]
    static let sourceHints: [(words: [String], source: MemorySource)] = [
        (["email", "emails", "mail", "inbox"], .mail),
        (["photo", "photos", "picture", "pictures", "screenshot", "screenshots"], .photo),
        (["meeting", "meetings", "event", "events", "appointment", "calendar"], .calendar),
        (["workout", "steps", "walked", "ran", "health"], .health),
        (["visited", "went", "been", "location", "parked"], .location),
        (["said", "voice", "conversation", "spoke"], .voice),
    ]

    public func search(query: String, in records: [MemoryRecord], filter: ProfileFilter = .all, now: Date = Date(),
                       limit: Int = 50, semantic: (@Sendable (String, MemoryRecord) -> Double)? = nil) -> [SearchHit] {
        let (interval, remainder) = QueryTimeParser(calendar: calendar).parse(query, now: now)
        var tokens = Self.tokenize(remainder)

        // Source hints narrow the pool but are also dropped from the scored terms.
        var hintedSources = Set<MemorySource>()
        for hint in Self.sourceHints where tokens.contains(where: { hint.words.contains($0.raw) }) {
            hintedSources.insert(hint.source)
        }
        tokens.removeAll { t in Self.sourceHints.contains { $0.words.contains(t.raw) } }

        var pool = records.filter { filter.includes($0.profileIDs) }
        if let interval { pool = pool.filter { interval.contains($0.date) } }
        if !hintedSources.isEmpty {
            let narrowed = pool.filter { hintedSources.contains($0.source) }
            if !narrowed.isEmpty { pool = narrowed }
        }
        guard !pool.isEmpty else { return [] }

        // No content words left ("yesterday", "emails from last week"): newest first.
        if tokens.isEmpty {
            return pool.sorted { $0.date > $1.date }.prefix(limit).map {
                SearchHit(id: $0.id, score: 0.1, snippet: Self.snippet(for: $0, terms: []))
            }
        }

        let docs = pool.map { (record: $0, terms: Self.tokenize($0.searchableText + " " + $0.title).map(\.stem)) }
        let avgLen = max(1, Double(docs.map { $0.terms.count }.reduce(0, +)) / Double(docs.count))
        var docFreq: [String: Int] = [:]
        for d in docs {
            for term in Set(tokens.map(\.stem)) where d.terms.contains(where: { Self.termMatches($0, term) > 0 }) {
                docFreq[term, default: 0] += 1
            }
        }

        var hits: [SearchHit] = []
        for d in docs {
            var score = 0.0
            var matched: [String] = []
            for t in tokens {
                var tf = 0.0
                for term in d.terms { tf += Self.termMatches(term, t.stem) }
                guard tf > 0 else { continue }
                matched.append(t.raw)
                let df = Double(docFreq[t.stem] ?? 0)
                let idf = log(1 + (Double(docs.count) - df + 0.5) / (df + 0.5))
                let k1 = 1.5, b = 0.75
                score += idf * (tf * (k1 + 1)) / (tf + k1 * (1 - b + b * Double(d.terms.count) / avgLen))
            }
            let sem = semantic?(query, d.record) ?? 0
            if score == 0 && sem < 0.55 { continue }
            // Coverage: matching more of the query beats matching one rare word.
            score *= 0.6 + 0.4 * Double(matched.count) / Double(max(1, tokens.count))
            score += sem * 2
            // Gentle recency: a month-old memory scores ~10% lower than one from today.
            let ageDays = max(0, now.timeIntervalSince(d.record.date) / 86400)
            score *= 1 / (1 + ageDays / 300)
            hits.append(SearchHit(id: d.record.id, score: score, snippet: Self.snippet(for: d.record, terms: matched)))
        }
        return hits.sorted { $0.score > $1.score }.prefix(limit).map { $0 }
    }

    // MARK: tokenising

    struct Token { var raw: String; var stem: String }

    static func tokenize(_ text: String) -> [Token] {
        var out: [Token] = []
        var word = ""
        func flush() {
            if word.count > 1 || Int(word) != nil, !stopWords.contains(word) { out.append(Token(raw: word, stem: stem(word))) }
            word = ""
        }
        for ch in text.lowercased() {
            if ch.isLetter || ch.isNumber || ch == "'" { if ch != "'" { word.append(ch) } } else { flush() }
        }
        flush()
        return out
    }

    static func stem(_ w: String) -> String {
        guard w.count > 4 else { return w }
        for suffix in ["ing", "ed", "es", "s"] where w.hasSuffix(suffix) && w.count - suffix.count >= 3 {
            return String(w.dropLast(suffix.count))
        }
        return w
    }

    /// 1 for an exact stem match, 0.5 for a prefix match on longer words.
    static func termMatches(_ docTerm: String, _ queryStem: String) -> Double {
        if docTerm == queryStem { return 1 }
        if queryStem.count >= 4, docTerm.hasPrefix(queryStem) { return 0.5 }
        return 0
    }

    static func snippet(for r: MemoryRecord, terms: [String]) -> String {
        let body = r.body.trimmingCharacters(in: .whitespacesAndNewlines)
        if body.isEmpty { return r.placeName ?? "" }
        let sentences = body.split(whereSeparator: { ".!?\n".contains($0) }).map { $0.trimmingCharacters(in: .whitespaces) }
        if let hit = sentences.first(where: { s in terms.contains { s.lowercased().contains($0) } }) { return String(hit.prefix(140)) }
        return String(body.prefix(140))
    }
}
