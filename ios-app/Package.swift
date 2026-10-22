// swift-tools-version: 5.10
import PackageDescription

// RecallCore holds all platform-independent logic (to-do extraction, profile
// assignment, search, route matching, routines). It only imports Foundation so
// it builds and tests on Linux CI as well as on macOS.
let package = Package(
    name: "RecallCore",
    platforms: [.iOS(.v17), .macOS(.v14)],
    products: [.library(name: "RecallCore", targets: ["RecallCore"])],
    targets: [
        .target(name: "RecallCore"),
        .testTarget(name: "RecallCoreTests", dependencies: ["RecallCore"]),
    ]
)
