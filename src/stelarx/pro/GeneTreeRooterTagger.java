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

/**
 * Runs ASTRAL-Pro3's rooting/tagging pass ({@code -T}) without involving a shell.
 *
 * <p>The requested output is committed only after ASTRAL-Pro3 exits successfully
 * and produces one non-empty Newick line per non-empty input line. A failed run
 * therefore cannot replace a previous valid tagged-tree file.</p>
 */
public final class GeneTreeRooterTagger {
    private static final String EXECUTABLE_RELATIVE = "ASTER-Linux/bin/astral-pro3";

    private GeneTreeRooterTagger() {}

    /** Result of one successful rooting/tagging run. */
    public record Result(Path executable, Path output, int treeCount) {}

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
                "ASTRAL-Pro3 executable is missing or not executable: " + executable);
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
            drain(process.getInputStream());
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new IOException("STELAR-Pro rooting/tagging failed (backend exit code "
                    + exitCode + ")");
            }
            if (!Files.isRegularFile(temporary) || Files.size(temporary) == 0L) {
                throw new IOException(
                    "STELAR-Pro rooting/tagging produced no tagged gene trees");
            }

            int inputTrees = countNonEmptyLines(input);
            int outputTrees = countNonEmptyLines(temporary);
            if (inputTrees == 0) {
                throw new IllegalArgumentException("Input gene-tree file is empty: " + input);
            }
            if (outputTrees != inputTrees) {
                throw new IOException("STELAR-Pro rooting/tagging produced " + outputTrees
                    + " tree(s) for " + inputTrees + " input tree(s)");
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

    private static void drain(InputStream input) throws IOException {
        byte[] buffer = new byte[8192];
        while (input.read(buffer) >= 0) { /* suppress backend output */ }
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
