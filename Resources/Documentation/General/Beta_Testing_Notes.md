# Triple_A — Beta testing notes

Thank you for testing Triple_A. This build is a **beta**: it is stable and fully
functional, but the goal of this round is to exercise it on real, varied datasets
and surface anything that needs fixing before a final release.

## What this build is

Triple_A (Advanced Amplicon Analysis) is a reproducible platform for amplicon
sequencing analysis: taxonomic profiling, community ecology statistics,
supervised multivariate methods, differential/compositional abundance, and
reference-genome-based functional inference — all from a single Shiny launcher.

**Recently added**: graphical summaries for statistics that
previously produced only tables (pairwise-PERMANOVA heatmap, beta-dispersion
boxplot, ANOSIM gauge, PERMANOVA-variance plot, ANCOM-BC2 volcano, MaAsLin forest);
a functional-evidence **confidence system** (High / Medium / Low / Insufficient)
with diagnostic-module completeness; new curated biogeochemical functions
(phosphorus, chlorine, iron, arsenic cycles and periplasmic nitrate reduction);
and an automatic self-contained HTML report plus PDF export.

## Changes in this build

**Experimental design.** The design can now be read from a metadata column
instead of the column order. Group membership is then matched by sample
identifier, so the abundance columns may be in any order, groups may be
interleaved, and **groups may have different numbers of replicates**. The
previous consecutive-block mode remains available and unchanged as the default.
Group order follows the metadata row order and sets the reference level of every
pairwise comparison.

**Analysis locking.** An analysis whose requirements are not met is disabled and
explains, under *How to unlock*, the steps that enable it. The same rules are
re-checked when the run is launched.

**Four statistical corrections.** If you have results from an earlier build,
these values need recomputing:

- `R2Y` and `Q2` of PLS-DA were computed without centring the response matrix and
  came out at or below zero even for a model classifying correctly;
- ANCOM-BC2 received rounded percentages instead of real counts when the table
  was declared as counts, losing the sequencing depth it models;
- PLS-DA and sPLS-DA were capped at one component fewer than the number of
  groups, which flattened every two-group score plot into a line;
- the RDA summary row *Constrained variance (%)* summed per-axis percentages and
  therefore always reported 100. It is now *Constrained variance (% of total
  inertia)* and equals 100 x R2.

The estimated-execution-time indicator was removed: it was not reliable.

## Getting started

1. Unzip the package. A folder `Triple_A` appears.
2. Make sure **R** is installed (R 4.4.x recommended).
3. Open and run **`Run_Triple_A.R`**. It opens the Launcher in your browser.
4. The Launcher checks required R packages and offers to install anything missing.
   Most analyses only need common CRAN packages; a few advanced ones need
   Bioconductor packages (see below).

This beta ships with a **pre-filled reference cache** (`Cache/GenomeCache.sqlite`
and `Cache/GFF`), so functional-potential analyses reuse reference lookups that
have already been resolved instead of re-querying NCBI for every taxon.

## Known limitations and things to watch

- **Functional potential (reference-genome / NCBI path)** is the most complex,
  network-dependent part and the least covered by automated tests. If a taxon
  cannot be resolved to a representative RefSeq genome, or NCBI is unreachable
  and the taxon is not already cached, that taxon is reported as **Insufficient
  evidence** — this is expected behaviour, not a crash. Please report actual
  errors (with the message), not the Insufficient-evidence label itself.
- **First-time genome/GFF lookups are serial** and can be slow on large datasets.
  Later runs reuse the cache. (Parallelisation is planned for a later version.)
- **FASTQ module** requires the Bioconductor package **dada2** (large download).
  The Launcher offers to install it. It loads only when you press *Process FASTQ*,
  so the module itself opens instantly. ANCOM-BC2 and MaAsLin similarly need their
  Bioconductor packages, offered on demand.
- **Automatic PDF report** needs a Chrome/Chromium browser (via the `chromote`
  package). If it is unavailable, the **self-contained `Triple_A_report.html`** is
  still generated — open it and use your browser's *Print → Save as PDF*.
- **Windows long paths**: keep the project folder in a reasonably short path. Very
  deep output paths can exceed the Windows 260-character limit and make a figure
  fail to save.
- **Interpretation**: functional calls, completeness and confidence tiers describe
  *inferred genomic potential* from a representative reference genome, not measured
  gene expression or in-situ activity. Manganese, selenium and anaerobic methane
  oxidation are intentionally **not** included yet (their marker genes are not
  reliably identifiable by name and would produce false positives).

## How to report a bug

Please include as much of the following as you can — it makes issues far faster to
reproduce and fix:

1. **What you did**: which module, which analysis, and the parameters you set.
2. **What you expected vs. what happened** (a screenshot helps for UI issues).
3. **The exact error text** (copy it, don't paraphrase).
4. **Logs**: attach or paste the relevant log from `Results/Logs/` and, for a
   specific run, the files under that run's `…/metadata` folder
   (`Run_metadata.json`, `Session_information.txt`).
5. **Your data shape** (not the data itself unless you can share it): number of
   samples, number of treatments and replicates, and the taxonomy format
   (single lineage column vs. separate ranks).
6. **Environment**: operating system and R version
   (run `R.version.string` and `sessionInfo()`).

Reproducibility metadata is written automatically for every run (parameters,
package versions, methods and session info), so pointing us at the run folder is
usually enough.

Thank you — every report makes the final release better.
