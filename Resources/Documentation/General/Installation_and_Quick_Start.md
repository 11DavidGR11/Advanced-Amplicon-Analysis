# Installation and quick start

## Requirements

Use a current 64-bit installation of R (R 4.3 or later is recommended) and RStudio Desktop. Internet access is required for first-time package installation and for downloading new reference annotations. Triple_A itself runs locally and opens in your web browser.

## Install the launcher packages

Open RStudio and run:

```r
install.packages(c("shiny", "bslib", "callr", "httpuv"))
```

Additional packages are requested only when a workflow needs them. Bioconductor analyses such as DADA2, ANCOM-BC2 or MaAsLin2 may require Bioconductor installation. Follow the package message shown by Triple_A when an optional dependency is missing.

## Validate the installation

From the Triple_A folder, run:

```r
source("Tools/Diagnostics/Verify_Installation.R")
```

The validator checks the required folder structure, parses the R source files and reports missing packages. A missing optional package does not prevent unrelated analyses from running.

## Start Triple_A

```r
source("Run_Triple_A.R")
```

The launcher opens in the default browser. Each workflow starts in an independent local R process. Keep the launcher and R session open while using a child application.

## Recommended first analysis

1. Open **Biological Analysis**.
2. Import a taxonomic abundance table and declare its abundance scale.
3. Confirm taxonomy and sample columns.
4. Import sample metadata when group comparisons or environmental analyses are required.
5. Choose the **Sample ID column** and, if the design comes from the metadata, the **grouping column**. Assign a role to the remaining columns.
6. Check the design summary, for example `3 group(s): Control (n=5), Treated (n=7)`. It states how the experiment was interpreted and is the last chance to catch a mis-declared design.
7. Review input validation.
8. Select a small set of analyses first, such as alpha diversity, PCoA and PERMANOVA. Analyses whose requirements are not met stay locked and explain how to unlock them.
9. Run the analyses and inspect both tables and figures.
10. Download the complete result package.

## Safe file handling

Do not rename or separate internal folders. Store original input data outside the program folder or in a clearly named project folder. Do not manually edit `Cache/GenomeCache.sqlite`. Use **Cache Manager** for cache operations.

## Closing the application

Close child-application browser tabs and then stop the R process from RStudio. Results already written to `Results/` remain available.
