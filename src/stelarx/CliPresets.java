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
                // S3 is S2 plus greedy-consensus candidates and therefore uses
                // the same cross-tree recombination stage.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            case "4", "s4" -> {
                cfg.setSearchSpace(Config.SearchSpace.S4);
                // S4 keeps S3's candidate pipeline, replacing its guide
                // similarity for completion, UPGMA, and consensus with triplets.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            case "5", "s5" -> {
                cfg.setSearchSpace(Config.SearchSpace.S5);
                // S5 adds A-Pro greedy-tree candidates to S1; recombination is optional.
            }
            case "6", "s6" -> {
                cfg.setSearchSpace(Config.SearchSpace.S6);
                // S6 independently adds rooted triplet-greedy trees to S1.
            }
            default -> throw new IllegalArgumentException(
                "unknown search space '" + value + "' (expected S1-S6 or 1-6)");
        }
    }

    private static String normalize(String value) {
        return value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
    }
}
