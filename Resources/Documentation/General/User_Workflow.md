# Triple_A — User workflow

## 1. Import and describe the data

Open **Data and metadata** and import the amplicon abundance table. Confirm the abundance scale, taxonomy columns, sample columns and biological replicate definition. Taxonomy may be supplied as a single lineage column or as separate ranks, depending on the imported table.

Sample metadata are optional. When supplied, use one row per sample and one column containing identifiers that match the abundance-table sample names. Triple_A does not restrict variable names. Assign every column as **identifier**, **experimental factor**, **environmental variable** or **ignore**. The selected roles determine which analyses are available.

## 2. Select analyses

The **Analyses and parameters** tab is organised by biological question. Enable only the analyses required for the study. Disabled controls indicate which input or package is missing.

- Diversity: alpha and beta diversity.
- Community structure: PCA, PCoA, NMDS and hierarchical clustering.
- Community comparison: PERMANOVA, pairwise PERMANOVA, ANOSIM and beta dispersion.
- Environmental analysis: envfit, RDA, partial RDA, dbRDA, partial dbRDA and variance partitioning.
- Taxon associations: ANCOM-BC2 and MaAsLin.
- Functional analysis: functional potential, abundance, differential functions and enrichment.

## 3. Validate and run

Review **Input validation**, correct any sample-name mismatch, choose the required figures and press **Run selected analyses**. A metadata-dependent method cannot run until compatible sample metadata have been imported and assigned.

## 4. Inspect and export

Use **Results** to browse tables and figures. Statistical tables that previously had no figure now also produce an optional graphical summary derived from the same result. Functional-potential runs report, per taxon, the diagnostic-module completeness and a confidence tier (High / Medium / Low / Insufficient evidence) alongside the classification. **Downloads** creates a ZIP of the current run and a report. Every run records parameters, package versions and session information, and writes a self-contained `Triple_A_report.html` plus, when a compatible browser engine is available, an automatic `Triple_A_report.pdf`.
