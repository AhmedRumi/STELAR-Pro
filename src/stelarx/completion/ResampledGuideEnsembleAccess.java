package stelarx.completion;

import stelarx.tree.Tree;

/** Test-only bridge to the package-private shared clade-set definition. */
final class ResampledGuideEnsembleAccess {
    private ResampledGuideEnsembleAccess() {}

    static boolean sameClades(Tree a, Tree b, int n) {
        return ResampledGuideEnsemble.cladeSet(a, n)
            .equals(ResampledGuideEnsemble.cladeSet(b, n));
    }
}
