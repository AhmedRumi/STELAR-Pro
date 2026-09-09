package stelarx.pro;

import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Runs ASTRAL-Pro3's rooting/tagging pass ({@code -T}) without involving a shell.
 *
 * <p>The requested output is committed only after ASTRAL-Pro3 exits successfully
 * and produces one non-empty Newick line per non-empty input line. A failed run
 * therefore cannot replace a previous valid tagged-tree file.</p>
 *
 * <p>The bundled executable is machine-specific (glibc symbol versions,
 * {@code -march=native}) and is built per machine by {@code ensure_backends.sh};
 * {@link #preflight} lets callers detect an unusable binary before any work is
 * done, and every failure keeps the backend's own output so a loader error such
 * as {@code version `GLIBC_2.38' not found} is reported verbatim.</p>
 */
public final class GeneTreeRooterTagger {
    private static final String EXECUTABLE_RELATIVE = "ASTER-Linux/bin/astral-pro3";
    /** How much of the backend's combined stdout/stderr tail is kept for error reports. */
    private static final int OUTPUT_TAIL_BYTES = 4096;
    private static final long PREFLIGHT_TIMEOUT_SECONDS = 60;

    /** Appended to every backend failure; the usual cause is a binary from another machine. */
    public static final String REBUILD_HINT =
        "The bundled ASTRAL-Pro3 binary is machine-specific and not tracked by git: build it "
        + "for this machine with ./ensure_backends.sh (run.sh does this automatically), or point "
        + "--astral-pro-executable / STELAR_PRO_EXECUTABLE at a working ASTRAL-Pro3.";

    private static Preflight cachedPreflight;
    private static String cachedPreflightKey;

    private GeneTreeRooterTagger() {}

    /** Result of one successful rooting/tagging run. */
    public record Result(Path executable, Path output, int treeCount) {}

    /** Outcome of the cheap startup check; {@code detail} is a short human-readable reason. */
    public record Preflight(Path executable, boolean usable, String detail) {}

    public static Result run(String inputFile, String outputFile,
                             String executableOverride, String mappingFile)
            throws IOException, InterruptedException {
        Path input = requireRegularFile(inputFile, "input gene-tree file");
        Path output = Path.of(outputFile).toAbsolutePath().normalize();
        if (input.equals(output)) {
            throw new IllegalArgumentException(
                "Rooted/tagged output file must differ from the unrooted input file");
        }

        Path mapping = mappingFile == null
            ? null : requireRegularFile(mappingFile, "gene-to-species mapping file");
        Path executable = resolveExecutable(executableOverride);
        if (!Files.isRegularFile(executable) || !Files.isExecutable(executable)) {
            throw new IllegalArgumentException(
                "ASTRAL-Pro3 executable is missing or not executable: " + executable
                + ". " + REBUILD_HINT);
        }

        Path parent = output.getParent();
        if (parent != null) Files.createDirectories(parent);
        String prefix = output.getFileName() == null ? "stelar-pro-tagged" : output.getFileName().toString();
        if (prefix.length() < 3) prefix = "stelar-pro-" + prefix;
        Path temporary = Files.createTempFile(parent, prefix + ".", ".tmp");

        List<String> command = new ArrayList<>();
        command.add(executable.toString());
        command.add("-T");
        if (mapping != null) {
            command.add("-a");
            command.add(mapping.toString());
        }
        command.add("-o");
        command.add(temporary.toString());
        command.add(input.toString());

        Process process = null;
        try {
            process = new ProcessBuilder(command)
                .redirectErrorStream(true)
                .start();
            String backendOutput = readTail(process.getInputStream(), OUTPUT_TAIL_BYTES);
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new IOException("STELAR-Pro rooting/tagging failed (backend exit code "
                    + exitCode + ")" + describeExit(exitCode) + "\n  executable: " + executable
                    + formatBackendOutput(backendOutput) + "\n  " + REBUILD_HINT);
            }
            if (!Files.isRegularFile(temporary) || Files.size(temporary) == 0L) {
                throw new IOException(
                    "STELAR-Pro rooting/tagging produced no tagged gene trees"
                    + formatBackendOutput(backendOutput));
            }

            int inputTrees = countNonEmptyLines(input);
            int outputTrees = countNonEmptyLines(temporary);
            if (inputTrees == 0) {
                throw new IllegalArgumentException("Input gene-tree file is empty: " + input);
            }
            if (outputTrees != inputTrees) {
                throw new IOException("STELAR-Pro rooting/tagging produced " + outputTrees
                    + " tree(s) for " + inputTrees + " input tree(s)"
                    + formatBackendOutput(backendOutput));
            }

            moveIntoPlace(temporary, output);
            return new Result(executable, output, outputTrees);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw e;
        } finally {
            if (process != null && process.isAlive()) process.destroyForcibly();
            Files.deleteIfExists(temporary);
        }
    }

    /**
     * Cheap startup check: resolve the executable and ask it for its help text.
     * ASTRAL-Pro3 answers {@code -h} with exit status 0 in a few milliseconds; a
     * binary built for another machine fails before {@code main} (dynamic-loader
     * exit status 1, or a signal for an unsupported instruction set). The result
     * is cached per executable so the banner, diagnostics, and the run itself
     * share one probe.
     */
    public static synchronized Preflight preflight(String executableOverride) {
        String key = executableOverride == null ? "" : executableOverride;
        if (cachedPreflight != null && key.equals(cachedPreflightKey)) return cachedPreflight;
        cachedPreflight = probe(resolveExecutable(executableOverride));
        cachedPreflightKey = key;
        return cachedPreflight;
    }

    private static Preflight probe(Path executable) {
        if (!Files.isRegularFile(executable)) {
            return new Preflight(executable, false, "executable not found");
        }
        if (!Files.isExecutable(executable)) {
            return new Preflight(executable, false, "file is not executable");
        }
        Process process = null;
        try {
            process = new ProcessBuilder(executable.toString(), "-h")
                .redirectErrorStream(true)
                .start();
            // Read on a helper thread so a hung executable cannot block the probe.
            final Process started = process;
            final String[] captured = {""};
            Thread reader = new Thread(() -> {
                try { captured[0] = readTail(started.getInputStream(), OUTPUT_TAIL_BYTES); }
                catch (IOException ignored) { /* partial output is fine for a diagnostic */ }
            }, "astral-pro3-preflight");
            reader.setDaemon(true);
            reader.start();
            if (!process.waitFor(PREFLIGHT_TIMEOUT_SECONDS, TimeUnit.SECONDS)) {
                process.destroyForcibly();
                return new Preflight(executable, false,
                    "did not answer -h within " + PREFLIGHT_TIMEOUT_SECONDS + " s");
            }
            reader.join(TimeUnit.SECONDS.toMillis(5));
            int exitCode = process.exitValue();
            if (exitCode != 0) {
                String firstLine = firstNonBlankLine(captured[0]);
                return new Preflight(executable, false,
                    "exit code " + exitCode + describeExit(exitCode) + " when asked for -h"
                    + (firstLine.isEmpty() ? "" : ": " + firstLine));
            }
            return new Preflight(executable, true, "answers -h (built for this machine)");
        } catch (IOException e) {
            return new Preflight(executable, false, "cannot start: " + e.getMessage());
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return new Preflight(executable, false, "interrupted while probing");
        } finally {
            if (process != null && process.isAlive()) process.destroyForcibly();
        }
    }

    /**
     * Resolution order: explicit CLI path, STELAR_PRO_EXECUTABLE, launcher home,
     * classpath checkout, then the current working directory.
     */
    public static Path resolveExecutable(String override) {
        if (override != null && !override.isBlank()) {
            return Path.of(override).toAbsolutePath().normalize();
        }
        String environment = System.getenv("STELAR_PRO_EXECUTABLE");
        if (environment != null && !environment.isBlank()) {
            return Path.of(environment).toAbsolutePath().normalize();
        }
        String home = System.getProperty("stelarpro.home");
        if (home != null && !home.isBlank()) {
            Path candidate = Path.of(home).resolve(EXECUTABLE_RELATIVE).normalize();
            if (Files.isRegularFile(candidate)) return candidate.toAbsolutePath();
        }
        try {
            Path classes = Path.of(GeneTreeRooterTagger.class.getProtectionDomain()
                .getCodeSource().getLocation().toURI()).toAbsolutePath().normalize();
            Path base = classes.getParent();
            if (base != null) {
                Path candidate = base.resolve(EXECUTABLE_RELATIVE).normalize();
                if (Files.isRegularFile(candidate)) return candidate;
            }
        } catch (URISyntaxException | SecurityException ignored) {
            // Fall through to the working-directory diagnostic path.
        }
        return Path.of(EXECUTABLE_RELATIVE).toAbsolutePath().normalize();
    }

    private static Path requireRegularFile(String value, String description) {
        Path path = Path.of(value).toAbsolutePath().normalize();
        if (!Files.isRegularFile(path)) {
            throw new IllegalArgumentException(description + " does not exist: " + path);
        }
        return path;
    }

    private static int countNonEmptyLines(Path path) throws IOException {
        int count = 0;
        try (var lines = Files.lines(path, StandardCharsets.UTF_8)) {
            var iterator = lines.iterator();
            while (iterator.hasNext()) {
                if (!iterator.next().isBlank()) count++;
            }
        }
        return count;
    }

    /** Drains the stream to EOF, keeping only its last {@code limit} bytes. */
    static String readTail(InputStream input, int limit) throws IOException {
        byte[] ring = new byte[limit];
        byte[] buffer = new byte[8192];
        int size = 0, next = 0;
        boolean truncated = false;
        int n;
        while ((n = input.read(buffer)) >= 0) {
            for (int i = 0; i < n; i++) {
                ring[next] = buffer[i];
                next = (next + 1) % limit;
                if (size < limit) size++; else truncated = true;
            }
        }
        byte[] tail = new byte[size];
        int start = size < limit ? 0 : next;
        for (int i = 0; i < size; i++) tail[i] = ring[(start + i) % limit];
        String text = new String(tail, StandardCharsets.UTF_8).strip();
        return truncated ? "…" + text : text;
    }

    private static String formatBackendOutput(String output) {
        if (output == null || output.isBlank()) return "\n  backend output: (none)";
        StringBuilder sb = new StringBuilder("\n  backend output:");
        for (String line : output.split("\\R")) {
            if (!line.isBlank()) sb.append("\n    ").append(line);
        }
        return sb.toString();
    }

    private static String firstNonBlankLine(String output) {
        if (output == null) return "";
        for (String line : output.split("\\R")) {
            if (!line.isBlank()) return line.strip();
        }
        return "";
    }

    /** Names the usual meaning of an exit status that did not come from ASTRAL-Pro3 itself. */
    private static String describeExit(int exitCode) {
        if (exitCode == 1) return " — ASTRAL-Pro3 never exits 1 itself; this is the dynamic loader"
            + " refusing the binary (typically a glibc symbol version missing on this machine)";
        if (exitCode == 127) return " — a required shared library was not found";
        if (exitCode == 132) return " — SIGILL: compiled with -march=native for another CPU";
        if (exitCode > 128) return " — terminated by signal " + (exitCode - 128);
        return "";
    }

    private static void moveIntoPlace(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }

}
