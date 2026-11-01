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
        check(cfg.getMaxSearchPasses() == 1, "one search pass by default");
        assertS1Defaults("default", cfg);

        check(!cfg.isKeepPolytomyDuringInference(),
            "default inference polytomy refinement");
        cfg.setKeepPolytomyDuringInference(true);
        check(cfg.isKeepPolytomyDuringInference(),
            "keep-polytomy-during-inference config");
        cfg.setKeepPolytomyDuringInference(false);
        for (int level = 1; level <= 6; level++) {
            cfg.setCompletionMethod(Config.CompletionMethod.DISTANCE);
            cfg.setStepBFastRestriction(false);
            CliPresets.applySearchSpace("S" + level, cfg);
            check(cfg.getSearchSpace() == Config.SearchSpace.valueOf("S" + level),
                "S" + level + " recognized");
            assertPresetDefaults("S" + level, cfg, level == 3 || level == 4);

            CliPresets.applySearchSpace(Integer.toString(level), cfg);
            check(cfg.getSearchSpace() == Config.SearchSpace.valueOf("S" + level),
                level + " recognized");
            assertPresetDefaults(Integer.toString(level), cfg, level == 3 || level == 4);
        }

        cfg.setMaxSearchPasses(3);
        check(cfg.getS5AdaptivePasses() == 2, "max passes includes initial pass");
        cfg.setS5AdaptivePasses(0);
        check(cfg.getMaxSearchPasses() == 1, "legacy zero extra passes");
        expectInvalid(() -> cfg.setMaxSearchPasses(0));
        expectInvalid(() -> CliPresets.applySearchSpace("S8", cfg));
        System.out.println("CLI preset mappings: PASS");
    }

    private static void assertS1Defaults(String label, Config cfg) {
        assertPresetDefaults(label, cfg, false);
    }

    private static void assertPresetDefaults(String label, Config cfg, boolean full) {
        check(cfg.getSearchMode() == (full
                ? Config.SearchMode.FULL : Config.SearchMode.LOCAL),
            label + " search mode");
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
