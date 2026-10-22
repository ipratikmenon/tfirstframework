# Total Recall (iOS)

A private life log for iPhone. It quietly collects what you saw, said and did (places, photos,
events, mail, to-dos spoken by voice), keeps it on the phone, lets you ask for it in plain
language, and brings the right item back when you head somewhere relevant.

Product requirements: see the PRD doc. This folder is the implementation.

## Layout

| Path | What it is |
| --- | --- |
| `Sources/RecallCore` | Platform-independent logic: to-do extraction, date parsing, profile assignment, search, routine learning, alert planning, route matching, Maps link parsing. Foundation only. |
| `Tests/RecallCoreTests` | 37 unit tests for the above. Run anywhere with `swift test`. |
| `App/` | The SwiftUI app: SwiftData models, capture services, views, App Intents. |
| `ShareExtension/` | Share-sheet target. Saves shared items and recalls what is noted near a shared Maps destination. |
| `project.yml` | XcodeGen description of the Xcode project. |

## Build

```sh
brew install xcodegen
cd ios-app
xcodegen generate          # writes TotalRecall.xcodeproj, Info.plists and entitlements
open TotalRecall.xcodeproj
```

Before running on a device, in `project.yml` set `DEVELOPMENT_TEAM`, and change the
`com.example.*` bundle ids and the `group.com.example.totalrecall` app group (also in
`App/Support/AppGroup.swift`) to ones you own. Requires Xcode 16 or later, iOS 17+.

Run the logic tests without Xcode (Linux or macOS): `swift test`.

## Gmail (optional)

Apple Mail cannot be read by any app. To read Gmail, create a Google Cloud OAuth client of type
iOS, put its client id in `GMAIL_CLIENT_ID` in `project.yml`, and sign in from Settings > Mail.
Only sender, subject and Google's short preview are stored. Google may require app verification
for the `gmail.readonly` scope outside your own test users.

## Action Button

Settings (iOS) > Action Button > Shortcut > Total Recall > Capture To-Dos.

## Status

`RecallCore` is compiled and tested on Linux (37 tests). The SwiftUI app target has not been compiled
in the environment that wrote it (no Xcode there), so expect to fix a few compiler errors on first build.
The CI workflow in `.github/workflows/ios-build.yml` builds it on a macOS runner.
