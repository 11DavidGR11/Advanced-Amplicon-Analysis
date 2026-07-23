# Results and interpretation

## Read tables before plots

Figures summarise results but do not replace numerical output. Check the number of retained samples and taxa, transformations, distance metric, model formula, effect size, raw P-value and adjusted P-value.

## Diversity and community structure

Alpha diversity describes within-sample diversity and depends on sequencing depth and filtering. PCoA and NMDS visualise dissimilarities but do not themselves test group differences. NMDS stress should be reported. PCA assumes Euclidean geometry and is normally applied after an appropriate transformation.

## Group comparisons

PERMANOVA evaluates whether group centroids differ in multivariate space. Interpret it together with beta dispersion, because unequal within-group dispersion can influence the result. Pairwise tests require multiple-testing correction. ANOSIM is a rank-based complementary test and should not be treated as interchangeable with effect-size estimates from PERMANOVA.

## Supervised discrimination

PLS-DA and sPLS-DA construct the axes that best separate the declared groups, so they separate groups even on random data. The score plot alone is therefore not evidence. Report the cross-validated performance and the permutation P-value.

Q2 is the value to cite as predictive ability, ahead of R2Y: R2Y measures fit on the training data and Q2 measures prediction on held-out samples, so a high R2Y with a low Q2 is the signature of overfitting. sPLS-DA is judged by its balanced error rate; the confusion matrix it reports is built on the training data and is labelled descriptive for that reason.

## Differential abundance

ANCOM-BC2 and MaAsLin2 answer different modelling questions. Confirm that the chosen model matches the experimental design and that covariates, repeated measures and reference levels are correctly defined. Prioritise adjusted P-values and effect sizes rather than raw significance alone.

ANCOM-BC2 models sequencing depth, so declaring the abundance table as raw counts gives the most reliable result. With proportions or percentages the depth information does not exist and counts per million are used instead.

The direction of each pairwise comparison follows the group order, which comes from the metadata file. Check which group acts as the reference before interpreting the sign of a fold change.

## Environmental analyses

RDA, dbRDA, envfit and variance partitioning require matched metadata and an adequate number of complete samples. Avoid models with too many predictors for the available sample size. Partial analyses condition on declared variables; they do not automatically remove all confounding.

For RDA, report the adjusted R2. The per-axis percentages in the constrained-variance figure are relative to the constrained part only and always sum to 100, so they do not state how much of the community variation the environment explains; the summary row `Constrained variance (% of total inertia)` does, and equals 100 x R2. Inspect the VIF sheet as well: it is reported but not enforced, and values above 10 make individual predictor coefficients uninterpretable.

## Functional potential

Functional-potential inference links classified taxa to curated genes in representative reference annotations. A positive result means that the function may be represented by organisms compatible with the taxonomic assignment. It does not prove that the sampled strain contains the gene or that the function was active.

## Minimum reporting information

Record the Triple_A version, input-file versions, filtering choices, abundance scale, transformation, distance metric, statistical formula, permutations, multiple-testing method, software-package versions and any samples removed during validation.
