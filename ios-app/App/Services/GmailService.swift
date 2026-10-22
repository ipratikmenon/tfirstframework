import Foundation
import AuthenticationServices
import CryptoKit
import UIKit
import RecallCore

/// Reads recent Gmail subjects and snippets with the user's own Google sign-in (OAuth with PKCE).
/// Apple Mail has no API for third-party apps, so connecting an account is the only route to mail.
/// Needs a Google Cloud iOS OAuth client: put its id in the GMAIL_CLIENT_ID build setting.
@MainActor
final class GmailService: NSObject, ObservableObject, ASWebAuthenticationPresentationContextProviding {
    static let shared = GmailService()

    @Published private(set) var accountEmail: String? = Keychain.get("gmail.email")
    @Published private(set) var lastError: String?

    static var clientID: String { (Bundle.main.object(forInfoDictionaryKey: "GmailClientID") as? String) ?? "" }
    var isConfigured: Bool { Self.clientID.hasSuffix(".apps.googleusercontent.com") }
    var isConnected: Bool { Keychain.get("gmail.refresh") != nil }
    var isEnabled: Bool { UserDefaults.standard.bool(forKey: Prefs.enabled("mail")) }

    private var redirectScheme: String {
        "com.googleusercontent.apps." + Self.clientID.replacingOccurrences(of: ".apps.googleusercontent.com", with: "")
    }
    private var redirectURI: String { redirectScheme + ":/oauth2redirect" }

    // MARK: sign in

    func connect() async {
        lastError = nil
        guard isConfigured else { lastError = "Add a Google OAuth client id (GMAIL_CLIENT_ID) to enable mail."; return }
        let verifier = Self.randomString(64)
        let challenge = Self.base64URL(Data(SHA256.hash(data: Data(verifier.utf8))))
        var comps = URLComponents(string: "https://accounts.google.com/o/oauth2/v2/auth")!
        comps.queryItems = [
            .init(name: "client_id", value: Self.clientID),
            .init(name: "redirect_uri", value: redirectURI),
            .init(name: "response_type", value: "code"),
            .init(name: "scope", value: "https://www.googleapis.com/auth/gmail.readonly email"),
            .init(name: "code_challenge", value: challenge),
            .init(name: "code_challenge_method", value: "S256"),
            .init(name: "access_type", value: "offline"),
            .init(name: "prompt", value: "consent"),
        ]
        do {
            let callback = try await authenticate(url: comps.url!)
            guard let code = URLComponents(url: callback, resolvingAgainstBaseURL: false)?.queryItems?.first(where: { $0.name == "code" })?.value
            else { lastError = "Sign-in was cancelled."; return }
            try await exchange(["code": code, "code_verifier": verifier, "grant_type": "authorization_code", "redirect_uri": redirectURI])
            try await loadProfile()
        } catch {
            lastError = error.localizedDescription
        }
    }

    func disconnect() {
        for key in ["gmail.access", "gmail.refresh", "gmail.expiry", "gmail.email"] { Keychain.remove(key) }
        accountEmail = nil
    }

    private func authenticate(url: URL) async throws -> URL {
        try await withCheckedThrowingContinuation { continuation in
            let session = ASWebAuthenticationSession(url: url, callbackURLScheme: redirectScheme) { callback, error in
                if let callback { continuation.resume(returning: callback) }
                else { continuation.resume(throwing: error ?? URLError(.cancelled)) }
            }
            session.presentationContextProvider = self
            session.prefersEphemeralWebBrowserSession = false
            session.start()
        }
    }

    nonisolated func presentationAnchor(for session: ASWebAuthenticationSession) -> ASPresentationAnchor {
        MainActor.assumeIsolated {
            UIApplication.shared.connectedScenes.compactMap { $0 as? UIWindowScene }.flatMap(\.windows).first { $0.isKeyWindow } ?? ASPresentationAnchor()
        }
    }

    // MARK: tokens

    private struct TokenResponse: Decodable {
        let access_token: String
        let expires_in: Double
        let refresh_token: String?
    }

    private func exchange(_ params: [String: String]) async throws {
        var request = URLRequest(url: URL(string: "https://oauth2.googleapis.com/token")!)
        request.httpMethod = "POST"
        request.setValue("application/x-www-form-urlencoded", forHTTPHeaderField: "Content-Type")
        var body = params
        body["client_id"] = Self.clientID
        request.httpBody = body.map { "\($0.key)=\(Self.formEncode($0.value))" }.joined(separator: "&").data(using: .utf8)
        let (data, response) = try await URLSession.shared.data(for: request)
        guard (response as? HTTPURLResponse)?.statusCode == 200 else { throw URLError(.userAuthenticationRequired) }
        let token = try JSONDecoder().decode(TokenResponse.self, from: data)
        Keychain.set(token.access_token, for: "gmail.access")
        Keychain.set(String(Date().timeIntervalSince1970 + token.expires_in - 60), for: "gmail.expiry")
        if let refresh = token.refresh_token { Keychain.set(refresh, for: "gmail.refresh") }
    }

    private func accessToken() async throws -> String {
        if let expiry = Keychain.get("gmail.expiry").flatMap(Double.init), expiry > Date().timeIntervalSince1970,
           let token = Keychain.get("gmail.access") { return token }
        guard let refresh = Keychain.get("gmail.refresh") else { throw URLError(.userAuthenticationRequired) }
        try await exchange(["refresh_token": refresh, "grant_type": "refresh_token"])
        guard let token = Keychain.get("gmail.access") else { throw URLError(.userAuthenticationRequired) }
        return token
    }

    // MARK: Gmail API

    private func get<T: Decodable>(_ path: String, query: [URLQueryItem] = []) async throws -> T {
        var comps = URLComponents(string: "https://gmail.googleapis.com/gmail/v1/users/me/" + path)!
        comps.queryItems = query.isEmpty ? nil : query
        var request = URLRequest(url: comps.url!)
        request.setValue("Bearer \(try await accessToken())", forHTTPHeaderField: "Authorization")
        let (data, response) = try await URLSession.shared.data(for: request)
        guard (response as? HTTPURLResponse)?.statusCode == 200 else { throw URLError(.badServerResponse) }
        return try JSONDecoder().decode(T.self, from: data)
    }

    private struct GmailProfile: Decodable { let emailAddress: String }
    private struct MessageList: Decodable { struct Ref: Decodable { let id: String }; let messages: [Ref]? }
    private struct Message: Decodable {
        struct Header: Decodable { let name: String; let value: String }
        struct Payload: Decodable { let headers: [Header]? }
        let id: String
        let snippet: String?
        let internalDate: String?
        let payload: Payload?
    }

    private func loadProfile() async throws {
        let profile: GmailProfile = try await get("profile")
        Keychain.set(profile.emailAddress, for: "gmail.email")
        accountEmail = profile.emailAddress
    }

    /// Stores subject, sender name and Gmail's own ~100 character snippet. Bodies are never downloaded.
    func sync(days: Int = 3) async {
        guard isEnabled, isConnected else { return }
        do {
            let list: MessageList = try await get("messages", query: [
                .init(name: "q", value: "newer_than:\(days)d -category:promotions -category:social -category:forums"),
                .init(name: "maxResults", value: "50")])
            for ref in list.messages ?? [] {
                if MemoryStore.shared.hasEntry(externalID: "mail-\(ref.id)") { continue }
                let m: Message = try await get("messages/\(ref.id)", query: [
                    .init(name: "format", value: "metadata"),
                    .init(name: "metadataHeaders", value: "From"),
                    .init(name: "metadataHeaders", value: "Subject")])
                let headers = m.payload?.headers ?? []
                func header(_ n: String) -> String? { headers.first { $0.name.caseInsensitiveCompare(n) == .orderedSame }?.value }
                let date = m.internalDate.flatMap(Double.init).map { Date(timeIntervalSince1970: $0 / 1000) } ?? Date()
                let from = header("From") ?? ""
                let sender = from.components(separatedBy: "<").first?.trimmingCharacters(in: CharacterSet(charactersIn: " \"")) ?? from
                let snippet = (m.snippet ?? "").replacingOccurrences(of: "&#39;", with: "'").replacingOccurrences(of: "&amp;", with: "&")
                let entry = MemoryEntry(date: date, source: .mail, title: header("Subject") ?? "(no subject)", body: snippet,
                                        people: sender.isEmpty ? [] : [sender], accountLabel: accountEmail, externalID: "mail-\(ref.id)")
                MemoryStore.shared.upsert(entry)
            }
            MemoryStore.shared.save()
            lastError = nil
        } catch {
            lastError = "Mail sync failed: \(error.localizedDescription)"
        }
    }

    // MARK: helpers

    static func randomString(_ length: Int) -> String {
        let chars = Array("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-._~")
        var generator = SystemRandomNumberGenerator()
        return String((0..<length).map { _ in chars.randomElement(using: &generator)! })
    }

    static func base64URL(_ data: Data) -> String {
        data.base64EncodedString().replacingOccurrences(of: "+", with: "-").replacingOccurrences(of: "/", with: "_").replacingOccurrences(of: "=", with: "")
    }

    static func formEncode(_ s: String) -> String {
        s.addingPercentEncoding(withAllowedCharacters: CharacterSet(charactersIn: "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-._~")) ?? s
    }
}
