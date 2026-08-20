# Rooted objective and data flow

For a candidate rooted split `A | B` and a rooted gene-tree node with child leaf
sets `M0 | M1`, define:

```text
a0 = |A intersect M0|    a1 = |A intersect M1|
b0 = |B intersect M0|    b1 = |B intersect M1|
```

The doubled number of agreeing rooted triplets is:

```text
2 w = a0*b1*(a0+b1-2) + a1*b0*(a1+b0-2)
```

This is exactly twice
`C(a0,2)b1 + a0C(b1,2) + C(a1,2)b0 + a1C(b0,2)`.
Dividing by two is exact. Summing this weight over gene nodes and then summing
chosen split weights over the species tree counts every agreeing rooted triplet
once, at its species-tree and gene-tree LCAs.

For a gene-tree polytomy with rooted child groups `Mi`, the linear-time form is:

```text
2 w = sum_i ai(ai-1)(sumB-bi) + bi(bi-1)(sumA-ai)
```

Only rooted child groups participate; taxa outside the node are not a child
group. Legacy third-part storage remains in compact structures solely to avoid
expanding the Java/CUDA ABI, and is ignored by the objective.

## Pipeline

```text
rooted Newick -> descendant clusters X -> rooted child transitions
              -> rooted gene child partitions -> split triplet weights
              -> memoized DP -> rooted Newick species tree
```

Incomplete gene trees contribute their native rooted triplets to scoring.
Completion is used only to enrich the candidate search space. Missing leaves
are attached next to their nearest present similarity anchor without reversing
any existing edge, preserving the supplied root and relationships among
observed taxa.
