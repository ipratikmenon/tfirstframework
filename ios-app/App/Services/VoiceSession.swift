import Foundation
import AVFoundation
import Speech
import RecallCore

/// Listens, understands, forgets. Audio goes straight from the microphone into on-device speech
/// recognition and is never written to disk; only the to-dos found in it are kept.
@MainActor
final class VoiceSession: NSObject, ObservableObject {
    static let shared = VoiceSession()

    enum Phase: Equatable { case idle, starting, listening, finishing, failed(String) }

    @Published private(set) var phase: Phase = .idle
    @Published private(set) var transcript = ""
    @Published private(set) var level: Float = 0
    @Published private(set) var drafts: [TodoDraft] = []
    /// Set after a session ends so the UI can confirm what was saved.
    @Published private(set) var lastSavedCount = 0

    private let engine = AVAudioEngine()
    private var request: SFSpeechAudioBufferRecognitionRequest?
    private var task: SFSpeechRecognitionTask?
    private var timeout: Task<Void, Never>?
    private var lastChange = Date()
    private let extractor = TodoExtractor()

    static let maxDuration: TimeInterval = 180
    static let silenceLimit: TimeInterval = 30

    var isActive: Bool { phase == .starting || phase == .listening || phase == .finishing }

    // MARK: permissions

    private func requestPermissions() async -> String? {
        let speech = await withCheckedContinuation { c in
            SFSpeechRecognizer.requestAuthorization { c.resume(returning: $0) }
        }
        guard speech == .authorized else { return "Allow Speech Recognition in Settings to capture to-dos by voice." }
        let mic = await withCheckedContinuation { c in
            AVAudioApplication.requestRecordPermission { c.resume(returning: $0) }
        }
        return mic ? nil : "Allow Microphone access in Settings to capture to-dos by voice."
    }

    // MARK: session

    func start() async {
        guard !isActive else { return }
        phase = .starting
        transcript = ""
        drafts = []
        lastSavedCount = 0
        if let problem = await requestPermissions() { phase = .failed(problem); return }

        guard let recognizer = SFSpeechRecognizer(locale: Locale.current) ?? SFSpeechRecognizer(locale: Locale(identifier: "en-US")),
              recognizer.isAvailable else {
            phase = .failed("Speech recognition isn't available right now.")
            return
        }
        // Privacy rule: audio is never sent to a server.
        guard recognizer.supportsOnDeviceRecognition else {
            phase = .failed("On-device speech recognition isn't available for \(Locale.current.identifier). Total Recall never sends audio off your device.")
            return
        }

        do {
            let session = AVAudioSession.sharedInstance()
            try session.setCategory(.record, mode: .measurement, options: [.duckOthers])
            try session.setActive(true, options: .notifyOthersOnDeactivation)

            let request = SFSpeechAudioBufferRecognitionRequest()
            request.shouldReportPartialResults = true
            request.requiresOnDeviceRecognition = true
            request.addsPunctuation = true
            request.taskHint = .dictation
            self.request = request

            let input = engine.inputNode
            let format = input.outputFormat(forBus: 0)
            input.removeTap(onBus: 0)
            input.installTap(onBus: 0, bufferSize: 1024, format: format,
                             block: Self.tapBlock(request: request) { [weak self] level in
                                 DispatchQueue.main.async { self?.level = level }
                             })
            engine.prepare()
            try engine.start()

            task = recognizer.recognitionTask(with: request, resultHandler: Self.resultHandler { [weak self] text, isFinal, failed in
                Task { @MainActor in self?.handle(text: text, isFinal: isFinal, failed: failed) }
            })
            phase = .listening
            lastChange = Date()
            startWatchdog()
        } catch {
            teardown()
            phase = .failed("Couldn't start the microphone: \(error.localizedDescription)")
        }
    }

    private func handle(text: String?, isFinal: Bool, failed: Bool) {
        guard phase == .listening || phase == .finishing else { return }
        if let text, text != transcript {
            transcript = text
            lastChange = Date()
            drafts = extractor.extract(from: text)
        }
        if isFinal || (failed && phase == .finishing) { finalise() }
        else if failed { stop() }
    }

    private func startWatchdog() {
        let started = Date()
        timeout?.cancel()
        timeout = Task { [weak self] in
            while !Task.isCancelled {
                try? await Task.sleep(nanoseconds: 2_000_000_000)
                guard let self, self.phase == .listening else { return }
                if Date().timeIntervalSince(started) > Self.maxDuration || Date().timeIntervalSince(self.lastChange) > Self.silenceLimit {
                    self.stop()
                    return
                }
            }
        }
    }

    /// Ends listening; the recogniser then delivers a final result and `finalise` saves the to-dos.
    func stop() {
        guard phase == .listening else { return }
        phase = .finishing
        timeout?.cancel()
        engine.stop()
        engine.inputNode.removeTap(onBus: 0)
        request?.endAudio()
        // If the recogniser never answers, finish anyway.
        Task { [weak self] in
            try? await Task.sleep(nanoseconds: 3_000_000_000)
            if self?.phase == .finishing { self?.finalise() }
        }
    }

    func cancel() {
        teardown()
        transcript = ""
        drafts = []
        phase = .idle
    }

    private func finalise() {
        guard phase == .finishing else { return }
        let found = extractor.extract(from: transcript)
        teardown()
        let store = MemoryStore.shared
        let defaults = UserDefaults.standard
        let auto = defaults.bool(forKey: Prefs.autoAddEnabled)
        let threshold = defaults.double(forKey: Prefs.autoAddThreshold)
        var titles: [String] = []
        for draft in found {
            let status: TodoStatus = auto && draft.confidence >= threshold ? .open : .draft
            store.addTodo(from: draft, status: status)
            titles.append(draft.title)
        }
        if !titles.isEmpty {
            // The timeline records that a session happened and what came of it, never the words spoken.
            let entry = MemoryEntry(date: Date(), source: .voice,
                                    title: titles.count == 1 ? "Captured 1 to-do" : "Captured \(titles.count) to-dos",
                                    body: titles.joined(separator: "\n"))
            store.upsert(entry)
        }
        store.save()
        lastSavedCount = titles.count
        // Forget the transcript immediately.
        transcript = ""
        drafts = []
        phase = .idle
        Task { await AlertEngine.shared.replan() }
    }

    private func teardown() {
        timeout?.cancel()
        if engine.isRunning { engine.stop() }
        engine.inputNode.removeTap(onBus: 0)
        request?.endAudio()
        task?.cancel()
        request = nil
        task = nil
        try? AVAudioSession.sharedInstance().setActive(false, options: .notifyOthersOnDeactivation)
        level = 0
    }

    // These closures are built outside the main actor on purpose: AVAudioEngine and the speech
    // framework call them on their own queues, and a main-actor closure would trap there.
    nonisolated static func tapBlock(request: SFSpeechAudioBufferRecognitionRequest,
                                     onLevel: @escaping @Sendable (Float) -> Void) -> AVAudioNodeTapBlock {
        { buffer, _ in
            request.append(buffer)
            onLevel(rms(buffer))
        }
    }

    nonisolated static func resultHandler(_ deliver: @escaping @Sendable (String?, Bool, Bool) -> Void)
        -> (SFSpeechRecognitionResult?, Error?) -> Void {
        { result, error in
            deliver(result?.bestTranscription.formattedString, result?.isFinal ?? false, error != nil)
        }
    }

    nonisolated static func rms(_ buffer: AVAudioPCMBuffer) -> Float {
        guard let data = buffer.floatChannelData?[0] else { return 0 }
        let n = Int(buffer.frameLength)
        guard n > 0 else { return 0 }
        var sum: Float = 0
        for i in 0..<n { sum += data[i] * data[i] }
        let rms = (sum / Float(n)).squareRoot()
        return min(1, rms * 8)      // roughly 0...1 for speech
    }
}
