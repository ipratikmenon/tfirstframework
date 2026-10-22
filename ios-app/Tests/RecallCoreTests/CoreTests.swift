import XCTest
@testable import RecallCore

final class ProfileTests: XCTestCase {
    let work = Profile.defaultWork()
    let assigner = ProfileAssigner(calendar: Fixture.calendar)

    func testWorkCalendarEventAssignedByCalendarName() {
        let r = MemoryRecord(date: Fixture.date(2026, 10, 24, 20), source: .calendar, title: "Dinner", calendarName: "Work")
        XCTAssertEqual(assigner.assign(r, profiles: [work]), [work.id])
    }

    func testTimeWindowAloneIsNotEnough() {
        // Wednesday 11:00, but nothing else says "work".
        let r = MemoryRecord(date: Fixture.date(2026, 10, 21, 11), source: .photo, title: "Lunch photo")
        XCTAssertTrue(assigner.assign(r, profiles: [work]).isEmpty)
    }

    func testTimeWindowPlusKeywordIsEnough() {
        let r = MemoryRecord(date: Fixture.date(2026, 10, 21, 11), source: .mail, title: "Invoice 42 overdue")
        XCTAssertEqual(assigner.assign(r, profiles: [work]), [work.id])
        let weekend = MemoryRecord(date: Fixture.date(2026, 10, 24, 11), source: .mail, title: "Invoice 42 overdue")
        XCTAssertTrue(assigner.assign(weekend, profiles: [work]).isEmpty, "keyword alone is below the threshold")
    }

    func testCustomProfileByPlaceAndAccount() {
        let home = Coordinate(latitude: 37.77, longitude: -122.41)
        let family = Profile(name: "Family", rules: [.place(center: home, radiusMeters: 150, name: "Home"),
                                                     .accountContains("gmail.com")], threshold: 0.8)
        let near = MemoryRecord(date: Fixture.now, source: .location, title: "Home", coordinate: Coordinate(latitude: 37.7705, longitude: -122.4101))
        let far = MemoryRecord(date: Fixture.now, source: .location, title: "Cafe", coordinate: Coordinate(latitude: 37.8, longitude: -122.4))
        XCTAssertEqual(assigner.assign(near, profiles: [family, work]), [family.id])
        XCTAssertTrue(assigner.assign(far, profiles: [family]).isEmpty)
        let mail = MemoryRecord(date: Fixture.now, source: .mail, title: "Hi", accountLabel: "me@gmail.com")
        XCTAssertEqual(assigner.assign(mail, profiles: [family]), [family.id])
    }

    func testScheduledProfileAndFilter() {
        XCTAssertEqual(assigner.scheduledProfile(at: Fixture.date(2026, 10, 21, 10), profiles: [work])?.id, work.id)
        XCTAssertNil(assigner.scheduledProfile(at: Fixture.date(2026, 10, 24, 10), profiles: [work]))
        XCTAssertTrue(ProfileFilter.unassigned.includes([]))
        XCTAssertFalse(ProfileFilter.profile(work.id).includes([]))
        XCTAssertTrue(ProfileFilter.all.includes([]))
    }
}

final class SearchTests: XCTestCase {
    let engine = SearchEngine(calendar: Fixture.calendar)
    let work = UUID()

    func records() -> [MemoryRecord] {
        [
            MemoryRecord(date: Fixture.date(2026, 10, 20, 14), source: .mail, title: "Re: Deposit for the flat",
                         body: "Hi, the landlord confirms the deposit of 1200 will be returned by Friday.", people: ["Marta"], profileIDs: []),
            MemoryRecord(date: Fixture.date(2026, 10, 14, 9), source: .calendar, title: "Meeting with Priya",
                         body: "Quarterly planning at the office", placeName: "Office", people: ["Priya"], profileIDs: [work]),
            MemoryRecord(date: Fixture.date(2026, 10, 20, 18), source: .location, title: "Blue Bottle Coffee",
                         placeName: "Blue Bottle Coffee"),
            MemoryRecord(date: Fixture.date(2026, 9, 2, 12), source: .photo, title: "Photo", body: "Receipt from the pharmacy"),
        ]
    }

    func testFindsEmailAboutDeposit() {
        let r = records()
        let hits = engine.search(query: "what did the landlord email about the deposit", in: r, now: Fixture.now)
        XCTAssertEqual(hits.first?.id, r[0].id)
        XCTAssertTrue(hits.first?.snippet.contains("deposit") == true)
    }

    func testPersonAndTimeRange() {
        let r = records()
        XCTAssertEqual(engine.search(query: "Priya", in: r, now: Fixture.now).first?.id, r[1].id)
        // "yesterday" is Tuesday 20 Oct: two records, newest first.
        let y = engine.search(query: "yesterday", in: r, now: Fixture.now)
        XCTAssertEqual(y.map(\.id), [r[2].id, r[0].id])
    }

    func testSourceHintNarrowsPool() {
        let r = records()
        let hits = engine.search(query: "emails this week", in: r, now: Fixture.now)
        XCTAssertEqual(hits.map(\.id), [r[0].id])
    }

    func testProfileFilterHidesOtherProfiles() {
        let r = records()
        XCTAssertEqual(engine.search(query: "planning", in: r, filter: .profile(work), now: Fixture.now).map(\.id), [r[1].id])
        XCTAssertTrue(engine.search(query: "planning", in: r, filter: .unassigned, now: Fixture.now).isEmpty)
    }

    func testStemmingAndPrefix() {
        let r = records()
        XCTAssertEqual(engine.search(query: "pharmacies receipts", in: r, now: Fixture.now).first?.id, r[3].id)
        XCTAssertEqual(engine.search(query: "plan", in: r, now: Fixture.now).first?.id, r[1].id)
    }

    func testTimeParser() {
        let p = QueryTimeParser(calendar: Fixture.calendar)
        let lastTue = p.parse("coffee last tuesday", now: Fixture.now)
        XCTAssertEqual(lastTue.interval?.start, Fixture.date(2026, 10, 20))
        XCTAssertEqual(lastTue.remainder.trimmingCharacters(in: .whitespaces), "coffee")
        XCTAssertEqual(p.parse("in september 2025", now: Fixture.now).interval?.start, Fixture.date(2025, 9, 1))
        XCTAssertEqual(p.parse("3 days ago", now: Fixture.now).interval != nil, true)
        XCTAssertNil(p.parse("hello", now: Fixture.now).interval)
    }
}

final class RouteTests: XCTestCase {
    // A straight east-west route about 2.2 km long at latitude 37.77.
    let route = (0...10).map { Coordinate(latitude: 37.77, longitude: -122.42 + Double($0) * 0.0025) }
    var destination: Coordinate { route.last! }

    func testDistances() {
        let a = Coordinate(latitude: 37.77, longitude: -122.42)
        let b = Coordinate(latitude: 37.78, longitude: -122.42)
        XCTAssertEqual(Geo.distance(a, b), 1111, accuracy: 15)
        let p = Coordinate(latitude: 37.7718, longitude: -122.41)   // ~200 m north of the route
        XCTAssertEqual(Geo.distance(from: p, toPath: route), 200, accuracy: 10)
    }

    func testMatchesOnRouteDestinationAndCategory() {
        let near = RouteCandidate(id: UUID(), kind: .todo, title: "Drop parcel", coordinate: Coordinate(latitude: 37.7718, longitude: -122.41))
        let far = RouteCandidate(id: UUID(), kind: .todo, title: "Far away", coordinate: Coordinate(latitude: 37.80, longitude: -122.41))
        let atDest = RouteCandidate(id: UUID(), kind: .memory, title: "Last time here", coordinate: Coordinate(latitude: 37.7701, longitude: destination.longitude))
        let pharmacy = RouteCandidate(id: UUID(), kind: .todo, title: "Pick up prescription", category: .pharmacy)
        let stamps = RouteCandidate(id: UUID(), kind: .todo, title: "Buy stamps", category: .postOffice)
        let places: [PlaceCategory: [NamedPlace]] = [
            .pharmacy: [NamedPlace(name: "Corner Pharmacy", coordinate: Coordinate(latitude: 37.7712, longitude: -122.405))],
            .postOffice: [NamedPlace(name: "Main PO", coordinate: Coordinate(latitude: 37.80, longitude: -122.405))],
        ]
        let m = RouteMatcher.match(candidates: [near, far, atDest, pharmacy, stamps], route: route, destination: destination,
                                   nearbyPlaces: places, now: Fixture.now)
        XCTAssertEqual(Set(m.map(\.candidate.title)), ["Drop parcel", "Last time here", "Pick up prescription"])
        XCTAssertEqual(m.first?.candidate.title, "Last time here")
        XCTAssertEqual(m.first?.reason, .atDestination)
        let rx = m.first { $0.candidate.title == "Pick up prescription" }
        if case .categoryAlongRoute(let name, _)? = rx?.reason { XCTAssertEqual(name, "Corner Pharmacy") } else { XCTFail("expected category match") }
    }

    func testDestinationCategoryMatchesTodoWithoutCoordinate() {
        let pharmacy = RouteCandidate(id: UUID(), kind: .todo, title: "Pick up prescription", category: .pharmacy)
        let m = RouteMatcher.match(candidates: [pharmacy], route: route, destination: destination, destinationCategory: .pharmacy, now: Fixture.now)
        XCTAssertEqual(m.first?.reason, .categoryAtDestination(.pharmacy))
    }
}

final class MapsLinkTests: XCTestCase {
    func testApple() {
        let d = MapsLinkParser.parse("https://maps.apple.com/?ll=37.7749,-122.4194&q=Blue%20Bottle")
        XCTAssertEqual(d?.name, "Blue Bottle")
        XCTAssertEqual(d?.coordinate, Coordinate(latitude: 37.7749, longitude: -122.4194))
        XCTAssertEqual(MapsLinkParser.parse("https://maps.apple.com/?daddr=1+Infinite+Loop,+Cupertino")?.address, "1 Infinite Loop, Cupertino")
    }

    func testGoogle() {
        let d = MapsLinkParser.parse("https://www.google.com/maps/place/Ferry+Building/@37.7955,-122.3937,17z/data=abc")
        XCTAssertEqual(d?.name, "Ferry Building")
        XCTAssertEqual(d?.coordinate, Coordinate(latitude: 37.7955, longitude: -122.3937))
        XCTAssertEqual(MapsLinkParser.parse("https://maps.google.com/?q=37.1,-122.2")?.coordinate?.latitude, 37.1)
    }

    func testRejectsOthers() {
        XCTAssertNil(MapsLinkParser.parse("https://example.com/maps"))
        XCTAssertNil(MapsLinkParser.parse("not a url"))
    }
}

final class RoutineAndAlertTests: XCTestCase {
    let gym = Coordinate(latitude: 37.7800, longitude: -122.4100)

    /// The user visited the gym every Thursday 18:00-19:00 for six weeks.
    func thursdayVisits() -> [Visit] {
        (1...6).map { week in
            // Thursdays before Wed 21 Oct 2026: 15 Oct, 8 Oct, 1 Oct, 24 Sep, 17 Sep, 10 Sep
            Visit(coordinate: gym, arrival: Fixture.date(2026, 10, 15, 18).addingTimeInterval(Double(-(week - 1) * 7) * 86400),
                  departure: Fixture.date(2026, 10, 15, 19).addingTimeInterval(Double(-(week - 1) * 7) * 86400))
        }
    }

    func testRoutinePredictsThursdayGym() {
        let model = RoutineModel(visits: thursdayVisits(), calendar: Fixture.calendar)
        let thursdayNoon = Fixture.date(2026, 10, 22, 12)
        let p = model.predictions(from: thursdayNoon, hours: 8)
        XCTAssertEqual(p.count, 1)
        XCTAssertEqual(Geo.distance(p[0].coordinate, gym), 0, accuracy: 200)
        XCTAssertEqual(p[0].probability, 1.0, accuracy: 0.01)
        // Not a gym day.
        XCTAssertTrue(model.predictions(from: Fixture.date(2026, 10, 21, 12), hours: 8).isEmpty)
    }

    func testPlannerAlertsBeforeUsualGymTimeAndForDueTodo() {
        let model = RoutineModel(visits: thursdayVisits(), calendar: Fixture.calendar)
        let now = Fixture.date(2026, 10, 22, 12)
        let shaker = TodoRecord(title: "Bring protein shaker", coordinate: gym, placeName: "Gym")
        let due = TodoRecord(title: "Send invoice", due: Fixture.date(2026, 10, 22, 14), dueHasTime: true)
        let unrelated = TodoRecord(title: "Water plants")
        let alerts = AlertPlanner.plan(todos: [shaker, due, unrelated], routine: model, knownPlaces: [:], now: now, horizonHours: 8)
        XCTAssertEqual(Set(alerts.map(\.title)), ["Bring protein shaker", "Send invoice"])
        let gymAlert = alerts.first { $0.title == "Bring protein shaker" }
        XCTAssertEqual(gymAlert?.fireDate, Fixture.date(2026, 10, 22, 17, 45))
        XCTAssertEqual(alerts.first { $0.title == "Send invoice" }?.fireDate, Fixture.date(2026, 10, 22, 13, 30))
    }

    func testPlannerHonoursCapCooldownAndProfile() {
        let now = Fixture.date(2026, 10, 22, 12)
        let work = UUID()
        let t1 = TodoRecord(title: "A", due: now.addingTimeInterval(3600), profileIDs: [work])
        let t2 = TodoRecord(title: "B", due: now.addingTimeInterval(3600), profileIDs: [])
        let t3 = TodoRecord(title: "C", due: now.addingTimeInterval(3600), profileIDs: [work], lastAlerted: now.addingTimeInterval(-600))
        let workOnly = AlertPlanner.plan(todos: [t1, t2, t3], routine: nil, knownPlaces: [:], now: now, filter: .profile(work))
        XCTAssertEqual(workOnly.map(\.title), ["A"])
        XCTAssertTrue(AlertPlanner.plan(todos: [t1, t2], routine: nil, knownPlaces: [:], now: now, dailyCap: 2, alreadyAlertedToday: 2).isEmpty)
    }

    func testDigest() {
        let r = [
            MemoryRecord(date: Fixture.now, source: .location, title: "Cafe", placeName: "Cafe"),
            MemoryRecord(date: Fixture.now, source: .location, title: "Office", placeName: "Office"),
            MemoryRecord(date: Fixture.now, source: .photo, title: "Photo"),
            MemoryRecord(date: Fixture.now, source: .health, title: "8,412 steps"),
        ]
        XCTAssertEqual(DailyDigest.summary(for: r, newTodos: 2),
                       "You visited 2 places, took 1 photo, noted 2 to-dos and recorded 8,412 steps.")
        XCTAssertEqual(DailyDigest.summary(for: []), "Nothing captured yet.")
    }
}

final class SnapshotTests: XCTestCase {
    func testQuickMatchAtDestinationAndByCategory() throws {
        let cafe = Coordinate(latitude: 37.7749, longitude: -122.4194)
        let snapshot = ExtensionSnapshot(items: [
            .init(id: UUID(), kind: .todo, title: "Ask barista about beans", coordinate: Coordinate(latitude: 37.7751, longitude: -122.4195)),
            .init(id: UUID(), kind: .todo, title: "Pick up prescription", category: .pharmacy),
            .init(id: UUID(), kind: .todo, title: "Water plants"),
            .init(id: UUID(), kind: .memory, title: "Left umbrella here", coordinate: cafe),
        ])
        let blue = SharedDestination(name: "Blue Bottle Coffee", address: nil, coordinate: cafe)
        let m = QuickDestinationMatcher.match(destination: blue, snapshot: snapshot, now: Fixture.now)
        XCTAssertEqual(Set(m.map(\.candidate.title)), ["Ask barista about beans", "Left umbrella here"])
        XCTAssertEqual(QuickDestinationMatcher.summary(for: m, destinationName: "Blue Bottle"), m.map(\.candidate.title).prefix(3).joined(separator: ", "))

        let pharm = SharedDestination(name: "Walgreens Pharmacy", address: nil, coordinate: Coordinate(latitude: 37.9, longitude: -122.5))
        let p = QuickDestinationMatcher.match(destination: pharm, snapshot: snapshot, now: Fixture.now)
        XCTAssertEqual(p.map(\.candidate.title), ["Pick up prescription"])

        // JSON round trip, as the extension reads it from disk.
        let data = try JSONEncoder().encode(snapshot)
        XCTAssertEqual(try JSONDecoder().decode(ExtensionSnapshot.self, from: data).items.count, 4)
    }
}
