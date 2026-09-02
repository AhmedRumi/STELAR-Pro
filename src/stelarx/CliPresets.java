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
            case "2", "s2" -> cfg.setSearchSpace(Config.SearchSpace.S2);
            case "3", "s3" -> cfg.setSearchSpace(Config.SearchSpace.S3);
            default -> throw new IllegalArgumentException(
                "unknown search space '" + value + "' (expected S1-S3 or 1-3)");
        }
    }

    private static String normalize(String value) {
        return value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
    }
}
