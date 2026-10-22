import Foundation

/// Kinds of place a to-do can be tied to ("buy stamps" -> post office).
public enum PlaceCategory: String, Codable, CaseIterable, Sendable {
    case pharmacy, grocery, gym, bank, postOffice, gasStation, hardware
    case coffee, restaurant, airport, office, school, hospital, laundry, library, florist

    public var displayName: String {
        switch self {
        case .pharmacy: return "Pharmacy"
        case .grocery: return "Grocery store"
        case .gym: return "Gym"
        case .bank: return "Bank"
        case .postOffice: return "Post office"
        case .gasStation: return "Gas station"
        case .hardware: return "Hardware store"
        case .coffee: return "Coffee shop"
        case .restaurant: return "Restaurant"
        case .airport: return "Airport"
        case .office: return "Office"
        case .school: return "School"
        case .hospital: return "Hospital"
        case .laundry: return "Dry cleaner"
        case .library: return "Library"
        case .florist: return "Florist"
        }
    }

    public var symbolName: String {
        switch self {
        case .pharmacy: return "cross.case.fill"
        case .grocery: return "cart.fill"
        case .gym: return "dumbbell.fill"
        case .bank: return "building.columns.fill"
        case .postOffice: return "envelope.fill"
        case .gasStation: return "fuelpump.fill"
        case .hardware: return "hammer.fill"
        case .coffee: return "cup.and.saucer.fill"
        case .restaurant: return "fork.knife"
        case .airport: return "airplane"
        case .office: return "building.2.fill"
        case .school: return "graduationcap.fill"
        case .hospital: return "cross.fill"
        case .laundry: return "tshirt.fill"
        case .library: return "books.vertical.fill"
        case .florist: return "leaf.fill"
        }
    }

    /// Query used with a map search to find places of this kind.
    public var searchQuery: String { displayName.lowercased() }

    /// Words in a sentence that imply this kind of place.
    var keywords: [String] {
        switch self {
        case .pharmacy: return ["pharmacy", "prescription", "medicine", "medication", "vitamins", "drugstore", "ibuprofen"]
        case .grocery: return ["groceries", "grocery", "supermarket", "milk", "eggs", "bread", "butter", "vegetables", "fruit", "coffee beans", "cereal", "rice", "chicken", "cheese"]
        case .gym: return ["gym", "workout", "yoga class"]
        case .bank: return ["bank", "deposit", "withdraw cash", "atm", "cheque", "check deposit"]
        case .postOffice: return ["post office", "stamps", "parcel", "package", "mail the", "ship the"]
        case .gasStation: return ["gas", "petrol", "fuel", "fill up", "refuel", "gas station"]
        case .hardware: return ["hardware", "screws", "nails", "paint", "batteries", "light bulbs", "lightbulbs", "drill"]
        case .coffee: return ["coffee shop", "cafe", "café"]
        case .restaurant: return ["restaurant", "reservation", "book a table"]
        case .airport: return ["airport"]
        case .office: return ["office", "the office"]
        case .school: return ["school", "daycare", "pick up the kids"]
        case .hospital: return ["hospital", "clinic", "doctor", "dentist"]
        case .laundry: return ["dry cleaning", "dry cleaner", "laundry", "tailor"]
        case .library: return ["library", "return the book", "return books"]
        case .florist: return ["flowers", "florist", "bouquet"]
        }
    }

    /// First category whose keywords appear in the text (word-boundary match).
    public static func detect(in text: String) -> PlaceCategory? {
        let padded = " " + text.lowercased().filter { $0.isLetter || $0.isNumber || $0 == " " || $0 == "é" } + " "
        var best: (PlaceCategory, Int)?
        for category in PlaceCategory.allCases {
            for keyword in category.keywords where padded.contains(" \(keyword) ") || padded.contains(" \(keyword)s ") {
                // Prefer the longest (most specific) keyword.
                if best == nil || keyword.count > best!.1 { best = (category, keyword.count) }
            }
        }
        return best?.0
    }
}
