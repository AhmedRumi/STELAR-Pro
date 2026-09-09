package stelarx.pro;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/** Focused process-boundary and atomic-output checks for STELAR-Pro Phase 0. */
public final class GeneTreeRooterTaggerTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);

        Path input = work.resolve("unrooted genes.tre");
        Path mapping = work.resolve("gene species.map");
        Path output = work.resolve("rooted tagged.tre");
        Files.writeString(input, "(a1,a2,(b1,c1));\n((a1,b1),(a2,c1),d1);\n",
            StandardCharsets.UTF_8);
        Files.writeString(mapping, "a1 A\na2 A\nb1 B\nc1 C\nd1 D\n",
            StandardCharsets.UTF_8);

        Path success = executable(work.resolve("fake astral-pro3"), """
            #!/bin/sh
            set -eu
            [ "$1" = "-T" ] || exit 11
            [ "$2" = "-a" ] || exit 12
            [ "$3" = "%s" ] || exit 13
            [ "$4" = "-o" ] || exit 14
            output=$5
            [ "$6" = "%s" ] || exit 15
            printf '((A,A)D,(B,C));\n(A,(B,(A,C)))D;\n' > "$output"
            """.formatted(mapping, input));

        GeneTreeRooterTagger.Result result = GeneTreeRooterTagger.run(
            input.toString(), output.toString(), success.toString(), mapping.toString());
        check(result.treeCount() == 2, "tree count");
        check(Files.readString(output).contains(")D"), "duplication tag preserved");

        Files.writeString(output, "previous-valid-output\n", StandardCharsets.UTF_8);
        Path failure = executable(work.resolve("failing-astral-pro3"), """
            #!/bin/sh
            set -eu
            while [ "$1" != "-o" ]; do shift; done
            printf 'partial-output\n' > "$2"
            exit 17
            """);
        expectFailure(() -> GeneTreeRooterTagger.run(
            input.toString(), output.toString(), failure.toString(), null),
            "exit code 17");
        check(Files.readString(output).equals("previous-valid-output\n"),
            "failed run replaced existing output");

        // A binary built on another machine dies in the dynamic loader with exit
        // status 1 and a one-line stderr message; both must survive into the error.
        Path foreign = executable(work.resolve("foreign-astral-pro3"), """
            #!/bin/sh
            echo 'foreign-astral-pro3: /lib/x86_64-linux-gnu/libc.so.6: version GLIBC_9.9 not found' >&2
            exit 1
            """);
        expectFailure(() -> GeneTreeRooterTagger.run(
            input.toString(), output.toString(), foreign.toString(), null),
            "GLIBC_9.9 not found");
        expectFailure(() -> GeneTreeRooterTagger.run(
            input.toString(), output.toString(), foreign.toString(), null),
            "ensure_backends.sh");

        Path helpful = executable(work.resolve("helpful astral-pro3"), """
            #!/bin/sh
            [ "$1" = "-h" ] && { echo 'usage: astral-pro3 [ --ARG_NAME ARG_VALUE ... ] input_file' >&2; exit 0; }
            exit 3
            """);
        GeneTreeRooterTagger.Preflight usable = GeneTreeRooterTagger.preflight(helpful.toString());
        check(usable.usable(), "preflight accepts a backend that answers -h: " + usable.detail());
        GeneTreeRooterTagger.Preflight broken = GeneTreeRooterTagger.preflight(foreign.toString());
        check(!broken.usable() && broken.detail().contains("GLIBC_9.9"),
            "preflight reports the loader message: " + broken.detail());
        check(broken.detail().contains("exit code 1"),
            "preflight reports the exit status: " + broken.detail());
        GeneTreeRooterTagger.Preflight absent =
            GeneTreeRooterTagger.preflight(work.resolve("missing").toString());
        check(!absent.usable() && absent.detail().contains("not found"),
            "preflight reports a missing executable: " + absent.detail());

        expectFailure(() -> GeneTreeRooterTagger.run(
            input.toString(), input.toString(), success.toString(), mapping.toString()),
            "must differ");
        expectFailure(() -> GeneTreeRooterTagger.run(
            input.toString(), output.toString(), work.resolve("missing").toString(), null),
            "missing or not executable");

        System.out.println("STELAR-Pro rooting/tagging process boundary: PASS");
    }

    private static Path executable(Path path, String content) throws IOException {
        Files.writeString(path, content, StandardCharsets.UTF_8);
        if (!path.toFile().setExecutable(true, true)) {
            throw new IOException("could not make test executable: " + path);
        }
        return path;
    }

    private static void expectFailure(ThrowingRunnable action, String expected) throws Exception {
        try {
            action.run();
            throw new AssertionError("expected failure containing: " + expected);
        } catch (IllegalArgumentException | IOException e) {
            check(e.getMessage().contains(expected), "unexpected failure: " + e.getMessage());
        }
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }

    @FunctionalInterface
    private interface ThrowingRunnable { void run() throws Exception; }
}
