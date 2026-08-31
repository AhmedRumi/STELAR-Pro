package stelarx.pro;

import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.List;

/** Uniquifies copies, resolves polytomies, then restores species labels. */
public final class GeneTreePolytomyResolver {
    private static final String SCRIPT_RELATIVE = "scripts/arb_resolve_polytomies.py";

    private GeneTreePolytomyResolver() {}

    public record Result(Path script, Path output, int treeCount) {}

    public static Result run(String inputFile, String outputFile)
            throws IOException, InterruptedException {
        Path input = requireRegularFile(inputFile, "input gene-tree file");
        Path output = Path.of(outputFile).toAbsolutePath().normalize();
        if (input.equals(output)) {
            throw new IllegalArgumentException(
                "Resolved output file must differ from the input gene-tree file");
        }

        Path script = resolveScript();
        if (!Files.isRegularFile(script)) {
            throw new IllegalArgumentException(
                "Polytomy-resolution script is missing: " + script);
        }

        Path parent = output.getParent();
        if (parent != null) Files.createDirectories(parent);
        String prefix = output.getFileName() == null
            ? "stelar-pro-resolved" : output.getFileName().toString();
        if (prefix.length() < 3) prefix = "stelar-pro-" + prefix;
        Path temporary = Files.createTempFile(parent, prefix + ".", ".tmp");

        String python = System.getenv().getOrDefault("STELAR_PRO_PYTHON", "python3");
        Process process = null;
        try {
            process = new ProcessBuilder(List.of(
                    python, script.toString(), input.toString(), temporary.toString()))
                .redirectErrorStream(true)
                .start();
            drain(process.getInputStream());
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new IOException("STELAR-Pro polytomy resolution failed (exit code "
                    + exitCode + ")");
            }

            int inputTrees = countNonEmptyLines(input);
            int outputTrees = countNonEmptyLines(temporary);
            if (inputTrees == 0) {
                throw new IllegalArgumentException("Input gene-tree file is empty: " + input);
            }
            if (outputTrees != inputTrees) {
                throw new IOException("STELAR-Pro polytomy resolution produced "
                    + outputTrees + " tree(s) for " + inputTrees + " input tree(s)");
            }

            moveIntoPlace(temporary, output);
            return new Result(script, output, outputTrees);
        } catch (InterruptedException error) {
            Thread.currentThread().interrupt();
            throw error;
        } finally {
            if (process != null && process.isAlive()) process.destroyForcibly();
            Files.deleteIfExists(temporary);
        }
    }

    /** Locate the script relative to the launcher, classpath, or working tree. */
    public static Path resolveScript() {
        String home = System.getProperty("stelarx.home");
        if (home != null && !home.isBlank()) {
            Path candidate = Path.of(home).resolve(SCRIPT_RELATIVE).normalize();
            if (Files.isRegularFile(candidate)) return candidate.toAbsolutePath();
        }
        try {
            Path classes = Path.of(GeneTreePolytomyResolver.class.getProtectionDomain()
                .getCodeSource().getLocation().toURI()).toAbsolutePath().normalize();
            Path base = classes.getParent();
            if (base != null) {
                Path candidate = base.resolve(SCRIPT_RELATIVE).normalize();
                if (Files.isRegularFile(candidate)) return candidate;
            }
        } catch (URISyntaxException | SecurityException ignored) {
            // Fall through to the working-directory diagnostic path.
        }
        return Path.of(SCRIPT_RELATIVE).toAbsolutePath().normalize();
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

    private static void drain(InputStream input) throws IOException {
        byte[] buffer = new byte[8192];
        while (input.read(buffer) >= 0) { /* suppress script output */ }
    }

    private static void moveIntoPlace(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException error) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
