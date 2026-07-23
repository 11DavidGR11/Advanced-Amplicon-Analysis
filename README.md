# Triple A — Advanced Amplicon Analysis

An integrated, reproducible framework for amplicon sequencing analysis: from raw reads or an
existing abundance table through to statistical evidence and biological interpretation, inside a
single environment.

Triple A is organised around the biological questions a study asks rather than around the
statistical methods it uses. An analysis is offered when the experimental design can support it,
and is locked — with an explanation of how to unlock it — when it cannot.

## What it does

- **17 analyses** producing **44 catalogued outputs** (32 figures and 12 tables)
- **55 curated biological functions** for reference-genome-based functional inference
- Community ecology (PCA, PCoA, NMDS, PERMANOVA, ANOSIM, beta dispersion, alpha diversity)
- Supervised multivariate methods (PLS-DA, sPLS-DA) with cross-validated performance
- Constrained ordination (RDA, dbRDA, partial RDA, envfit, variance partitioning)
- Compositional differential abundance (ANCOM-BC2, MaAsLin2) and pairwise testing
- Optional DADA2 preprocessing from raw FASTQ files

Every figure is backed by an Excel workbook containing the numbers behind it, and every run writes
its own reproducibility metadata: parameters, methods, package versions and a file manifest.

## Requirements

R 4.4.x. Nothing else needs to be installed by hand — the Launcher checks the required packages and
offers to install anything missing. Most analyses need only common CRAN packages; ANCOM-BC2,
MaAsLin2 and the FASTQ module need Bioconductor packages, which are offered on demand.

Optional: a Chrome or Chromium installation enables PDF export of the analysis report. Without it
the self-contained HTML report is still produced.

## Quick start

```r
source("Run_Triple_A.R")
```

The Launcher opens in your browser. From there you choose the workflow that matches your starting
data. Each application runs in its own local R process.

Keep the project folder in a reasonably short path: deep output paths can exceed the Windows
260-character limit and make a figure fail to save.

## The applications

| Application | Purpose |
|---|---|
| **Launcher** | Entry point. Checks dependencies and opens each application. |
| **Biological Analysis** | The main interface: import, validation, experimental design, analysis selection, execution and result browsing. |
| **FASTQ Pipeline** | DADA2 preprocessing: quality assessment, trimming, ASV inference, taxonomy. |
| **Amplicon Integrator** | Combines compatible abundance tables from several sequencing runs. |
| **Reference Annotation Cache** | Inspects, imports, merges, verifies and exports the local SQLite/GFF cache. |
| **Biological Function Builder** | Adds new biological functions and pathways declaratively, without editing R code. |

## Repository layout

```
Run_Triple_A.R        Single entry point
AAApp/                Launcher, the five applications and the analysis engine
  Common/Engine/      All the science, independent of any interface
Plugins/              Analysis plugins, discovered at runtime
Resources/            Contextual documentation and the curated functional database
Cache/                Reference-genome cache (not versioned; see .gitignore)
Results/              Run outputs, one directory per run
QualityAssurance/     Automated test suite (development only)
Tools/                Diagnostics and release packaging (development only)
```

`QualityAssurance/` and `Tools/` are development trees and are deliberately excluded from the
distribution an end user receives.

## Design principles

**The engine does not depend on the interface.** Every analysis in `AAApp/Common/Engine/` can be
called directly from an R script through the public API in `Triple_A.R`, which is what makes results
reproducible outside the application and testable without a browser.

**The distribution is self-contained and portable.** Code, documentation, reference cache and
results all live inside the project folder, so moving or copying that single directory moves the
whole working environment with it.

## Documentation

User guides live in `Resources/Documentation/` and are also shown by the Help button inside each
application. Start with `General/Installation_and_Quick_Start.md`, then
`General/Input_Data_Formats.md`.

`Triple_A_User_Manual.docx` is the full user manual: it covers every analysis and how to read every
figure and table the platform produces.

## Quality assurance

```bash
Rscript QualityAssurance/Run_All_Tests.R
```

97 tests covering scientific correctness, interface behaviour, caching and regression, including
numerical checks against reference results. The suite writes a report to
`QualityAssurance/reports/`.

Set `TRIPLE_A_TEST_LEVEL` to `smoke`, `functional`, `regression` or `release` to control the depth
of the run (default: `release`, which runs everything).

## Building a release

```bash
Rscript Tools/Release/Build_Clean_Release.R
```

Produces a distribution archive containing only what an end user needs. Set
`TRIPLE_A_KEEP_CACHE=TRUE` to ship the pre-resolved reference cache with it, which is what the
beta-testing build does.

## Citing

See `CITATION.cff`, or use the "Cite this repository" button on GitHub.

## License

MIT — see [LICENSE](LICENSE).

## Author

**David Garrido Rodríguez**
Universidad de Valladolid, Institute of Sustainable Processes (ISP)
[ORCID 0000-0001-5180-0006](https://orcid.org/0000-0001-5180-0006) · david.garrido23@uva.es
