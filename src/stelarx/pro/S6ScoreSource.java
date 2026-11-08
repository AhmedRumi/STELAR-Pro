package stelarx.pro;

/** Scores S6 greedy proposals without prescribing how triplet support is stored. */
public interface S6ScoreSource extends AutoCloseable {
    @FunctionalInterface
    interface Focus {
        /** The focused taxon must be one of the three arguments. */
        long support(int a, int b, int outgroup);
    }

    int taxa();

    /** Borrow an independent mutable scoring state for a guide or proposal. */
    S6ScoreSource fork();

    /** Allow bounded reuse of focused slices shared by one guide batch. */
    default void prepareGuideFocus(int[] selected) {}
    default void clearGuideFocus() {}

    /** Support involving one taxon currently being inserted. */
    Focus focus(int taxon);

    /** Restrict an insertion query to taxa in its fixed backbone. */
    default Focus focusOn(int taxon, int[] backboneTaxa) {
        return focus(taxon);
    }

    /** Scores AB|C, AC|B, and BC|A for three disjoint species sets. */
    long[] groupScores(int[] a, int[] b, int[] c);

    /** Exact split contribution used only to report a completed greedy tree's score. */
    long splitWeight(int[] a, int[] b);

    @Override void close();
}
