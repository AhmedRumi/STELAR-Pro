package stelarx;

/** Focused, dependency-free tests for the user-facing CLI mappings. */
public final class CliPresetsTest {
    public static void main(String[] args) {
        Config cfg = Config.getInstance();
        check(cfg.getSearchSpace() == Config.SearchSpace.S1,
            "default search space is S1");
        check(cfg.getWeightIntersectionMethod()
                == Config.WeightIntersectionMethod.SMALLER_SIDE_TRAVERSAL,
            "built-in duplicate-aware intersection method");
        assertS1Defaults("default", cfg);

        check(!cfg.isKeepPolytomyDuringInference(),
            "default inference polytomy refinement");
        cfg.setKeepPolytomyDuringInference(true);
        check(cfg.isKeepPolytomyDuringInference(),
            "keep-polytomy-during-inference config");
        cfg.setKeepPolytomyDuringInference(false);
        // S2/S3 are recognized names, but remain inert reserved markers. Main's
        // current-scope validation rejects them before any analysis starts.
        for (int level = 1; level <= 3; level++) {
            cfg.setCompletionMethod(Config.CompletionMethod.DISTANCE);
            cfg.setStepBFastRestriction(false);
            CliPresets.applySearchSpace("S" + level, cfg);
            check(cfg.getSearchSpace() == Config.SearchSpace.valueOf("S" + level),
                "S" + level + " recognized");
            assertS1Defaults("S" + level, cfg);

            CliPresets.applySearchSpace(Integer.toString(level), cfg);
            check(cfg.getSearchSpace() == Config.SearchSpace.valueOf("S" + level),
                level + " recognized");
            assertS1Defaults(Integer.toString(level), cfg);
        }

        expectInvalid(() -> CliPresets.applySearchSpace("S4", cfg));
        expectInvalid(() -> CliPresets.applySearchSpace("S8", cfg));
        System.out.println("CLI preset mappings: PASS");
    }

    private static void assertS1Defaults(String label, Config cfg) {
        check(cfg.getSearchMode() == Config.SearchMode.LOCAL, label + " search mode");
        check(!cfg.isAutoCompleteIncompleteTrees(), label + " autocomplete");
        check(!cfg.isConsensusExperimental(), label + " consensus");
        check(!cfg.isStepBQuadraticNnBalls(), label + " quadratic");
        check(!cfg.isStepBRandomLeftoverResolution(), label + " random resolution");
        check(!cfg.isStepBProcessLargePolytomies(), label + " large polytomies");
        check(!cfg.isResolveInputGeneTreePolytomies(), label + " input polytomies");
        check(cfg.getCompletionMethod() == Config.CompletionMethod.SIMILARITY,
            label + " completion method");
        check(cfg.isStepBFastRestriction(), label + " fast restriction");
    }

    private static void expectInvalid(Runnable action) {
        try {
            action.run();
            throw new AssertionError("expected IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
            // Expected.
        }
    }

    private static void check(boolean condition, String label) {
        if (!condition) throw new AssertionError("failed: " + label);
    }
}
