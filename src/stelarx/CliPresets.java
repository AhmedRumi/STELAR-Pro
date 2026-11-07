package stelarx;

import java.util.Locale;

/**
 * User-facing CLI presets.  These methods only translate friendly names into
 * the existing Config switches; they do not change inference algorithms.
 */
final class CliPresets {
    private CliPresets() {}

    static void applySearchSpace(String value, Config cfg) {
        String preset = normalize(value);

        // A preset chooses the candidate source and its default DP search mode.
        // Explicit cross-tree CLI flags are applied after all presets are parsed.
        cfg.setSearchMode(Config.SearchMode.LOCAL);
        cfg.setAutoCompleteIncompleteTrees(false);
        cfg.setCompletionMethod(Config.CompletionMethod.SIMILARITY);
        cfg.setConsensusExperimental(false);
        cfg.setStepBFastRestriction(true);
        cfg.setStepBQuadraticNnBalls(false);
        cfg.setStepBRandomLeftoverResolution(false);
        cfg.setStepBProcessLargePolytomies(false);
        cfg.setResolveInputGeneTreePolytomies(false);

        switch (preset) {
            case "1", "s1" -> cfg.setSearchSpace(Config.SearchSpace.S1);
            case "2", "s2" -> {
                cfg.setSearchSpace(Config.SearchSpace.S2);
                // DISCO/completion/UPGMA local transitions are included by default.
            }
            case "3", "s3" -> {
                cfg.setSearchSpace(Config.SearchSpace.S3);
                // S3 adds greedy-consensus candidates and enables cross-tree
                // recombination so synthesized clades can receive DP splits.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            case "4", "s4" -> {
                cfg.setSearchSpace(Config.SearchSpace.S4);
                // S4 uses S2's local candidate pipeline with a triplet guide.
                // Cross-tree transitions remain available by explicit flag.
            }
            case "5", "s5" -> {
                cfg.setSearchSpace(Config.SearchSpace.S5);
                // S5 uses A-Pro greedy and one UPGMA guide, excluding S1 candidates.
            }
            case "6", "s6" -> {
                cfg.setSearchSpace(Config.SearchSpace.S6);
                // S6 independently adds rooted triplet-greedy trees to S1.
            }
            case "7", "s7" -> {
                cfg.setSearchSpace(Config.SearchSpace.S7);
                // S7 combines A-Pro greedy and streamed DISCO/UPGMA candidates.
            }
            default -> throw new IllegalArgumentException(
                "unknown search space '" + value + "' (expected S1-S7 or 1-7)");
        }
    }

    private static String normalize(String value) {
        return value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
    }
}
