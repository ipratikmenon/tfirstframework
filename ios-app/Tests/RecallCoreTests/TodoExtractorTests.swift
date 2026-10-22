import XCTest
@testable import RecallCore

final class TodoExtractorTests: XCTestCase {
    let extractor = TodoExtractor(calendar: Fixture.calendar)

    func testRemindMeWithTomorrowTime() {
        let d = extractor.extract(from: "Remind me to call Sam tomorrow at 5pm.", now: Fixture.now)
        XCTAssertEqual(d.count, 1)
        XCTAssertEqual(d[0].title, "Call Sam")
        XCTAssertEqual(d[0].person, "Sam")
        XCTAssertEqual(d[0].due, Fixture.date(2026, 10, 22, 17))
        XCTAssertTrue(d[0].dueHasTime)
        XCTAssertGreaterThan(d[0].confidence, 0.9)
    }

    func testNeedItemBecomesGetAndDetectsGrocery() {
        let d = extractor.extract(from: "We need milk and eggs.", now: Fixture.now)
        XCTAssertEqual(d.first?.title, "Get milk and eggs")
        XCTAssertEqual(d.first?.category, .grocery)
    }

    func testOutOf() {
        let d = extractor.extract(from: "We're out of coffee beans", now: Fixture.now)
        XCTAssertEqual(d.first?.title, "Buy coffee beans")
    }

    func testCommitmentWithWeekday() {
        let d = extractor.extract(from: "I'll send the deck to Priya by Friday.", now: Fixture.now)
        XCTAssertEqual(d.first?.title, "Send the deck to Priya")
        XCTAssertEqual(d.first?.due, Fixture.date(2026, 10, 23, 9))
        XCTAssertEqual(d.first?.dueHasTime, false)
    }

    func testMultipleSentencesAndNoise() {
        let text = "That was a great lunch. I need to pick up my prescription today. The weather is nice. Don't forget to book the dentist."
        let d = extractor.extract(from: text, now: Fixture.now)
        XCTAssertEqual(d.map(\.title), ["Pick up my prescription", "Book the dentist"])
        XCTAssertEqual(d[0].category, .pharmacy)
        XCTAssertEqual(d[1].category, .hospital)
    }

    func testRejectsNegationQuestionsAndHypotheticals() {
        let text = "I won't go to the gym. Should I call Sam? If I have time I'll clean the garage. We don't need to buy bread."
        XCTAssertTrue(extractor.extract(from: text, now: Fixture.now).isEmpty)
    }

    func testRequestKeepsLowerConfidence() {
        let d = extractor.extract(from: "Could you send me the invoice?", now: Fixture.now)
        XCTAssertEqual(d.count, 1)
        XCTAssertLessThan(d[0].confidence, 0.7)
    }

    func testNamedPlace() {
        let d = extractor.extract(from: "I need to return the drill at Home Depot.", now: Fixture.now)
        XCTAssertEqual(d.first?.placeName, "Home Depot")
    }

    func testDedupesSameTask() {
        let d = extractor.extract(from: "I need to call Sam. I need to call Sam.", now: Fixture.now)
        XCTAssertEqual(d.count, 1)
    }
}

final class DateParserTests: XCTestCase {
    let p = DateParser(calendar: Fixture.calendar)

    func testInHours() {
        let r = p.parse("in 2 hours", relativeTo: Fixture.now)
        XCTAssertEqual(r?.date, Fixture.date(2026, 10, 21, 12))
        XCTAssertEqual(r?.hasTime, true)
    }

    func testTonightAndTomorrowMorning() {
        XCTAssertEqual(p.parse("tonight", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 21, 20))
        XCTAssertEqual(p.parse("tomorrow morning", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 22, 9))
    }

    func testTimeOnlyRollsToTomorrowWhenPassed() {
        XCTAssertEqual(p.parse("at 9am", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 22, 9))
        XCTAssertEqual(p.parse("at 3pm", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 21, 15))
    }

    func testBareHourAssumesDaytime() {
        XCTAssertEqual(p.parse("tomorrow at 4", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 22, 16))
    }

    func testNextWeekdayAndNextWeek() {
        XCTAssertEqual(p.parse("next monday", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 26, 9))
        XCTAssertEqual(p.parse("next week", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 26, 9))
        // Friday is two days away, so "next friday" means the one after.
        XCTAssertEqual(p.parse("next friday", relativeTo: Fixture.now)?.date, Fixture.date(2026, 10, 30, 9))
    }

    func testNoDate() {
        XCTAssertNil(p.parse("buy some milk", relativeTo: Fixture.now))
    }
}
