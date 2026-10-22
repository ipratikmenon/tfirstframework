import Foundation

/// A to-do found in spoken text, before the user confirms it.
public struct TodoDraft: Codable, Hashable, Sendable {
    public var title: String
    public var due: Date?
    public var dueHasTime: Bool
    public var category: PlaceCategory?
    public var placeName: String?
    public var person: String?
    /// 0...1. Drives the review queue and the optional auto-add threshold.
    public var confidence: Double

    public init(title: String, due: Date? = nil, dueHasTime: Bool = false, category: PlaceCategory? = nil,
                placeName: String? = nil, person: String? = nil, confidence: Double) {
        self.title = title
        self.due = due
        self.dueHasTime = dueHasTime
        self.category = category
        self.placeName = placeName
        self.person = person
        self.confidence = confidence
    }
}

/// Finds commitments and tasks in a transcript. Rule based and fully on device:
/// no model download, deterministic, easy to test.
public struct TodoExtractor: Sendable {
    public var dateParser: DateParser

    public init(calendar: Calendar = .current) {
        self.dateParser = DateParser(calendar: calendar)
    }

    struct Trigger {
        let pattern: String
        let confidence: Double
        var buildsNeed = false     // "we need milk" -> "Get milk"
        var outOf = false          // "we're out of milk" -> "Buy milk"
    }

    // Order matters: the first (strongest / most specific) trigger that matches wins.
    static let triggers: [Trigger] = [
        Trigger(pattern: #"\bremind me to\b"#, confidence: 0.95),
        Trigger(pattern: #"\b(?:don'?t|do not) forget to\b"#, confidence: 0.92),
        Trigger(pattern: #"\bremember to\b"#, confidence: 0.85),
        Trigger(pattern: #"\b(?:i|we) (?:really )?(?:need|have|has|got|gotta|must) to\b"#, confidence: 0.86),
        Trigger(pattern: #"\b(?:i've|i have|we've|we have) got to\b"#, confidence: 0.86),
        Trigger(pattern: #"\b(?:i|we)(?:'m| am| are|'re)? gotta\b"#, confidence: 0.86),
        Trigger(pattern: #"\b(?:we|i) (?:need|want) (?:some |more |a |an |to get )?(?=[a-z])"#, confidence: 0.7, buildsNeed: true),
        Trigger(pattern: #"\b(?:we|i)(?:'re| are|'m| am| have|'ve) (?:out of|run out of|ran out of|low on)\b"#, confidence: 0.7, outOf: true),
        Trigger(pattern: #"\bmake sure (?:to |that (?:i|we) )?"#, confidence: 0.72),
        Trigger(pattern: #"\bfollow up (?:with|on)\b"#, confidence: 0.7),
        Trigger(pattern: #"\b(?:i'll|i will|i'm going to|i am going to|i'm gonna|i am gonna|i'll go ahead and)\b"#, confidence: 0.72),
        Trigger(pattern: #"\b(?:we'll|we will|we're going to|we are going to)\b"#, confidence: 0.66),
        Trigger(pattern: #"\blet me\b"#, confidence: 0.6),
        Trigger(pattern: #"\b(?:can|could|would) you (?:please )?"#, confidence: 0.6),
        Trigger(pattern: #"\bplease\b"#, confidence: 0.55),
        Trigger(pattern: #"\b(?:i|we) (?:should|ought to)\b"#, confidence: 0.55),
    ]

    static let rejectPhrases = ["don't need to", "do not need to", "no need to", "won't", "will not",
                                "not going to", "never going to", "wish i", "wish we", "used to"]
    static let nameStopWords: Set<String> = ["The", "A", "An", "My", "Our", "Your", "His", "Her", "Their", "This", "That",
                                             "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday",
                                             "Today", "Tomorrow", "Tonight", "Next", "Some", "It", "Them", "Him", "Her", "Me", "Us"]
    static let personVerbs: Set<String> = ["call", "text", "email", "message", "ping", "ask", "tell", "meet", "remind",
                                           "thank", "visit", "invite", "congratulate", "ring"]

    public func extract(from transcript: String, now: Date = Date()) -> [TodoDraft] {
        var drafts: [TodoDraft] = []
        var seen = Set<String>()
        for sentence in Self.sentences(in: transcript) {
            guard var draft = analyse(sentence: sentence, now: now) else { continue }
            let key = draft.title.lowercased()
            if seen.contains(key) { continue }
            seen.insert(key)
            draft.title = Self.cleanTitle(draft.title)
            drafts.append(draft)
        }
        return drafts
    }

    // MARK: single sentence

    func analyse(sentence original: String, now: Date) -> TodoDraft? {
        let s = original.trimmingCharacters(in: .whitespacesAndNewlines)
        let lower = s.lowercased()
        guard s.split(separator: " ").count >= 3 else { return nil }
        if Self.rejectPhrases.contains(where: { lower.contains($0) }) { return nil }
        if lower.hasPrefix("if ") || lower.contains(", if ") || lower.hasPrefix("maybe ") || lower.hasPrefix("what if") { return nil }
        let isQuestion = s.hasSuffix("?")

        var matched: (Trigger, Range<String.Index>)?
        for t in Self.triggers {
            if let r = lower.range(of: t.pattern, options: .regularExpression) {
                matched = (t, r); break
            }
        }
        guard let (trigger, range) = matched else { return nil }
        // A question is only a task when it is a request ("could you send...").
        if isQuestion && trigger.confidence > 0.6 { return nil }

        // Work on the original-cased tail so names keep their capitals.
        let offset = lower.distance(from: lower.startIndex, to: range.upperBound)
        let tailStart = s.index(s.startIndex, offsetBy: min(offset, s.count))
        var tail = String(s[tailStart...]).trimmingCharacters(in: CharacterSet(charactersIn: " ,.-"))
        guard !tail.isEmpty else { return nil }

        if trigger.buildsNeed { tail = "Get " + tail }
        if trigger.outOf { tail = "Buy " + tail }

        var confidence = trigger.confidence
        var due: Date?
        var hasTime = false
        if let parsed = dateParser.parse(tail, relativeTo: now) {
            due = parsed.date
            hasTime = parsed.hasTime
            for p in parsed.phrases { tail = tail.replacingOccurrences(of: p, with: "", options: [.caseInsensitive]) }
            confidence += 0.05
        }
        tail = Self.tidy(tail)
        guard tail.count >= 3 else { return nil }

        let category = PlaceCategory.detect(in: tail)
        if category != nil { confidence += 0.04 }
        let placeName = Self.namedPlace(in: s)
        let person = Self.person(in: s)

        return TodoDraft(title: tail, due: due, dueHasTime: hasTime, category: category,
                         placeName: placeName, person: person, confidence: min(0.99, confidence))
    }

    // MARK: text helpers

    static func sentences(in text: String) -> [String] {
        var normalised = text.replacingOccurrences(of: "\n", with: ". ")
        for glue in [" and then ", " after that ", " also "] {
            normalised = normalised.replacingOccurrences(of: glue, with: ". ", options: .caseInsensitive)
        }
        var out: [String] = []
        var current = ""
        for ch in normalised {
            current.append(ch)
            if ch == "." || ch == "?" || ch == "!" || ch == ";" {
                out.append(current); current = ""
            }
        }
        if !current.trimmingCharacters(in: .whitespaces).isEmpty { out.append(current) }
        return out.map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }.filter { !$0.isEmpty }
    }

    static func tidy(_ s: String) -> String {
        var t = s
        for filler in [" or something", " or whatever", " I guess", " really quick", " real quick"] {
            t = t.replacingOccurrences(of: filler, with: "", options: .caseInsensitive)
        }
        t = t.replacingOccurrences(of: "  ", with: " ")
        let junk = CharacterSet(charactersIn: " ,.;:-!?").union(.whitespacesAndNewlines)
        t = t.trimmingCharacters(in: junk)
        // Trailing dangling prepositions left behind after removing a date ("send the deck by").
        let dangling: Set<String> = ["by", "on", "at", "in", "before", "until", "for", "this", "next", "tomorrow"]
        var words = t.split(separator: " ").map(String.init)
        while let last = words.last, dangling.contains(last.lowercased()) { words.removeLast() }
        return words.joined(separator: " ")
    }

    static func cleanTitle(_ s: String) -> String {
        guard let first = s.first else { return s }
        return first.uppercased() + s.dropFirst()
    }

    /// "pick up the parcel at Costco" -> "Costco"
    static func namedPlace(in sentence: String) -> String? {
        let pattern = #"\b(?:at|from|to) ((?:[A-Z][\w&'’-]*)(?: [A-Z][\w&'’-]*)*)"#
        guard let m = firstMatch(pattern, in: sentence) else { return nil }
        let name = m[1]
        let first = name.split(separator: " ").first.map(String.init) ?? name
        return nameStopWords.contains(first) ? nil : name
    }

    /// "call Sam about the invoice" -> "Sam"
    static func person(in sentence: String) -> String? {
        let words = sentence.split(separator: " ").map { $0.trimmingCharacters(in: CharacterSet(charactersIn: ",.;:!?")) }
        for (i, w) in words.enumerated() where personVerbs.contains(w.lowercased()) {
            var j = i + 1
            while j < words.count, ["to", "with", "up"].contains(words[j].lowercased()) { j += 1 }
            if j < words.count, let c = words[j].first, c.isUppercase, !nameStopWords.contains(words[j]) {
                return words[j]
            }
        }
        if let m = firstMatch(#"\b(?:follow up with|meet with|check in with) ([A-Z][a-z]+)"#, in: sentence) { return m[1] }
        return nil
    }
}
