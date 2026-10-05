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

        // A preset is a complete, deterministic configuration.  This also makes
        // command-line ordering unsurprising: options are applied left-to-right,
        // so a later legacy flag can fine-tune a preset and a later preset resets
        // all search-space flags.
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
                // S2 contains every local transition plus compatible cross-tree
                // recombinations over the enriched candidate cluster set.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            case "3", "s3" -> {
                cfg.setSearchSpace(Config.SearchSpace.S3);
                // S3 is S2 plus greedy-consensus candidates and therefore uses
                // the same cross-tree recombination stage.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            case "4", "s4" -> {
                cfg.setSearchSpace(Config.SearchSpace.S4);
                // S4 keeps S3's UPGMA/consensus guide and changes only the
                // incomplete-tree insertion matrix to rooted-triplet similarity.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            case "5", "s5" -> {
                cfg.setSearchSpace(Config.SearchSpace.S5);
                // S5 adds A-Pro greedy-tree candidates to native S1 candidates;
                // cross-tree transitions let STELAR's triplet DP recombine them.
                cfg.setSearchMode(Config.SearchMode.FULL);
            }
            default -> throw new IllegalArgumentException(
                "unknown search space '" + value + "' (expected S1-S5 or 1-5)");
        }
    }

    private static String normalize(String value) {
        return value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
    }
}
