import Foundation
import HealthKit
import RecallCore

/// Daily step totals and workouts. Read-only; HealthKit does not reveal whether read access was granted.
@MainActor
final class HealthService: ObservableObject {
    static let shared = HealthService()
    private let store = HKHealthStore()

    var isAvailable: Bool { HKHealthStore.isHealthDataAvailable() }
    var isEnabled: Bool { UserDefaults.standard.bool(forKey: Prefs.enabled("health")) }

    func requestAuthorization() async {
        guard isAvailable, let steps = HKQuantityType.quantityType(forIdentifier: .stepCount) else { return }
        try? await store.requestAuthorization(toShare: [], read: [steps, HKObjectType.workoutType()])
    }

    func sync(days: Int = 7) async {
        guard isEnabled, isAvailable else { return }
        let cal = Calendar.current
        for offset in 0..<days {
            guard let day = cal.date(byAdding: .day, value: -offset, to: cal.startOfDay(for: Date())),
                  let end = cal.date(byAdding: .day, value: 1, to: day),
                  let count = await steps(from: day, to: end), count > 0 else { continue }
            let formatted = NumberFormatter.localizedString(from: NSNumber(value: Int(count)), number: .decimal)
            let key = ISO8601DateFormatter().string(from: day)
            let entry = MemoryEntry(date: min(end, Date()).addingTimeInterval(-1), source: .health, title: "\(formatted) steps",
                                    externalID: "steps-\(key)")
            MemoryStore.shared.upsert(entry)
        }
        for w in await workouts(since: cal.date(byAdding: .day, value: -days, to: Date()) ?? Date()) {
            let km = w.totalDistance?.doubleValue(for: .meterUnit(with: .kilo))
            var title = Self.name(of: w.workoutActivityType)
            if let km, km > 0 { title += String(format: ", %.1f km", km) }
            let minutes = Int((w.duration / 60).rounded())
            let entry = MemoryEntry(date: w.startDate, endDate: w.endDate, source: .health, title: title,
                                    body: "\(minutes) minutes", externalID: "workout-\(w.uuid.uuidString)")
            MemoryStore.shared.upsert(entry)
        }
        MemoryStore.shared.save()
    }

    private func steps(from start: Date, to end: Date) async -> Double? {
        guard let type = HKQuantityType.quantityType(forIdentifier: .stepCount) else { return nil }
        return await withCheckedContinuation { continuation in
            let predicate = HKQuery.predicateForSamples(withStart: start, end: end, options: .strictStartDate)
            let query = HKStatisticsQuery(quantityType: type, quantitySamplePredicate: predicate, options: .cumulativeSum) { _, stats, _ in
                continuation.resume(returning: stats?.sumQuantity()?.doubleValue(for: .count()))
            }
            store.execute(query)
        }
    }

    private func workouts(since: Date) async -> [HKWorkout] {
        await withCheckedContinuation { continuation in
            let predicate = HKQuery.predicateForSamples(withStart: since, end: nil, options: .strictStartDate)
            let query = HKSampleQuery(sampleType: HKObjectType.workoutType(), predicate: predicate, limit: 50,
                                      sortDescriptors: [NSSortDescriptor(key: HKSampleSortIdentifierStartDate, ascending: false)]) { _, samples, _ in
                continuation.resume(returning: (samples as? [HKWorkout]) ?? [])
            }
            store.execute(query)
        }
    }

    static func name(of type: HKWorkoutActivityType) -> String {
        switch type {
        case .running: return "Run"
        case .walking: return "Walk"
        case .cycling: return "Ride"
        case .swimming: return "Swim"
        case .hiking: return "Hike"
        case .yoga: return "Yoga"
        case .functionalStrengthTraining, .traditionalStrengthTraining: return "Strength workout"
        default: return "Workout"
        }
    }
}
