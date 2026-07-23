# Biological analyses in Triple_A

## Diversity

**Question:** How diverse is each sample and how different are samples from one another? Alpha diversity reports richness and evenness within samples. Beta diversity calculates pairwise ecological dissimilarities.

## Community structure

**Question:** What are the main compositional patterns? PCA analyses transformed abundances. PCoA represents a selected distance matrix. NMDS preserves ranked dissimilarities and reports stress. Hierarchical clustering groups similar samples.

## Community comparison

**Question:** Are entire communities different between groups? PERMANOVA tests group or model effects and reports R². Pairwise PERMANOVA identifies which levels differ after a significant global result and adjusts P-values. ANOSIM is a rank-based alternative. Beta dispersion tests whether within-group variability differs and must be interpreted alongside PERMANOVA.

## Environmental analysis

**Question:** Which measured variables explain community composition? envfit overlays environmental associations on unconstrained ordinations. RDA is appropriate for linear models of transformed abundance. dbRDA constrains ecological distance matrices. Partial variants condition on confounders. Variance partitioning separates unique and shared contributions of predictor sets.

## Taxon associations

**Question:** Which taxa differ between groups or track environmental gradients? ANCOM-BC2 is a compositional differential-abundance model with bias correction. MaAsLin fits multivariable feature–metadata associations. Both require matched sample metadata and should report adjusted P-values.

## Functional analysis

**Question:** Which biological capabilities may be represented? Functional potential uses curated taxon-to-function evidence. Functional abundance weights functions by contributing taxa. Differential and enrichment analyses compare functional profiles or test over-representation. These are inferred potentials and are not direct measurements of gene expression or metabolic activity.

Each functional call is now accompanied by the **diagnostic-module completeness** (the fraction of the essential/diagnostic gene set detected) and a **confidence tier** — *High* (module essentially complete), *Medium*, *Low*, *No evidence detected*, or *Insufficient evidence* (too little of the module could be evaluated to judge the call). Diagnostic genes are treated as the essential module and are distinguished from supporting and accessory genes.

Version 4.0 adds curated functions for further biogeochemical cycles: phosphonate degradation (phosphorus), (per)chlorate reduction and organohalide respiration (chlorine), dissimilatory iron reduction and iron oxidation (iron), arsenate respiration, arsenite oxidation and arsenic detoxification (arsenic), and periplasmic nitrate reduction. Manganese, selenium and anaerobic methane oxidation were deliberately not added because their canonical markers are poorly conserved or share gene symbols with unrelated pathways, which would produce false positives under a gene-name search.

## Graphical summaries

Whenever an analysis produces a statistical table, Triple_A can also generate an optional publication-ready figure from the same result, without introducing new statistics: a pairwise-PERMANOVA heatmap, a beta-dispersion boxplot, an ANOSIM R gauge and a PERMANOVA explained-variance plot for community comparison, an ANCOM-BC2 volcano plot, and a MaAsLin coefficient forest plot.

## Choosing between overlapping methods

Triple_A deliberately offers several methods for some questions so you can match the analysis to your data and to the conventions of your field. When methods overlap, this is how to choose — you rarely need to run and report all of them.

**Which taxa differ between groups? (differential taxa)**

- **Differential abundance** (pairwise Wilcoxon or t-test on relative abundances): fast and exploratory, easy to read as volcano/MA plots. It does not account for compositionality, so treat borderline hits with caution and confirm them with one of the models below.
- **ANCOM-BC2**: bias-corrected and compositionality-aware. Preferred when you need a statistically defensible list of differentially abundant taxa. Requires biological replication.
- **MaAsLin2**: multivariable models (several metadata variables, optional random effects). Use it when the effect of interest must be adjusted for covariates, or with continuous/complex metadata.
- *Rule of thumb:* explore with Differential abundance, then confirm and report with ANCOM-BC2 (or MaAsLin2 when you need covariate adjustment).

**Which samples are similar? (ordination — PCA / PCoA / NMDS)**

All three place samples in two dimensions by overall similarity; choose one as your main figure rather than reporting all three.

- **PCA**: linear, on transformed abundances; emphasises dominant variance gradients. Suitable when the transformation makes Euclidean distance meaningful.
- **PCoA**: works on the ecological distance you select (e.g. Bray-Curtis) and stays faithful to it. The usual default in microbial ecology.
- **NMDS**: rank-based on the same distance; robust to non-linearity and reports a stress value for fit quality. Preferred when the distance structure is strongly non-linear.

**Do communities differ between groups? (PERMANOVA vs ANOSIM)**

- **PERMANOVA** is the recommended test: it reports an effect size (R²), handles multivariable designs and is the current standard. Read it together with the **beta-dispersion** diagnostic, because unequal within-group dispersion can by itself drive a significant PERMANOVA.
- **ANOSIM** answers the same question with a rank-based statistic (R). It is more sensitive to dispersion and provides no effect size, so it is largely superseded and is kept mainly for comparison or continuity with older studies. If you run only one, choose PERMANOVA.

**Environmental gradients (RDA vs dbRDA and partial variants)**

- **RDA** assumes approximately linear responses to the standardised environmental variables; **dbRDA** constrains an ecological *distance* matrix instead, so it suits non-Euclidean dissimilarities such as Bray-Curtis. Use the **partial** variants to condition out confounding factors, and **variance partitioning** to split explained variation between two predictor sets. **envfit** is a lighter option that only overlays environmental vectors on an unconstrained ordination.

**Supervised discrimination (PLS-DA vs sPLS-DA)**

- **PLS-DA** discriminates predefined groups using all taxa; **sPLS-DA** does the same but also selects a compact taxon signature, so prefer it when you want a short, interpretable list of discriminating taxa. Judge both by cross-validated performance, not by visual separation alone.
