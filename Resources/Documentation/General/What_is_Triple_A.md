# What Is Triple_A?

## Purpose

Triple_A (Advanced Amplicon Analysis) is a local R/Shiny platform for reproducible analysis of amplicon sequencing data. It connects four stages that are often handled with separate tools: FASTQ preprocessing, integration of amplicon abundance tables, biological and statistical analysis, and management of reusable reference annotations.

Triple_A is designed for 16S and 18S projects in which the user has sequencing reads or taxonomic abundance tables together with sample metadata. The application produces analysis-specific tables, figures, logs and run metadata inside the project `Results/` directory.

## Main applications

- **Biological Analysis** imports an abundance table and sample design, validates them, prepares taxonomy and runs selected analyses.
- **MultiAmplicon Integrator** combines compatible abundance tables and exports a tab-separated table for downstream analysis.
- **Reference Annotation Cache** inspects, imports, merges, verifies, cleans and exports the local SQLite/GFF cache.
- **Contextual help** is available directly inside each application and is loaded from that application's documentation folder.
- **FASTQ Pipeline** processes raw reads with DADA2. This is the most recently incorporated module and should be considered an upstream route into the established table-based workflow.

## End-to-end workflow

1. Start with raw FASTQ files or an existing taxonomic abundance table.
2. When necessary, process FASTQ files to obtain an ASV abundance table and taxonomy.
3. Optionally combine multiple amplicon tables in MultiAmplicon Integrator.
4. Prepare a sample-design table with treatment and experimental metadata.
5. Import the abundance and design tables into Biological Analysis.
6. Review validation, taxonomy parsing and sample matching before running analyses.
7. Select analyses and scientifically appropriate parameters.
8. Inspect exported tables before interpreting plots.
9. Record the software release identifier, input data, parameters and limitations in the study report.

## Project architecture

`Run_Triple_A.R` is the single entry point. It opens a launcher that starts each Shiny application in an isolated R process. Shared scientific code is located in `AAApp/Common/Engine/`; applications are in `AAApp/`; biological extensions are discovered from `Plugins/`; immutable resources are in `Resources/`; runtime cache content is in `Cache/`; and generated outputs are in `Results/`.

The runtime reference cache contains:

```text
Cache/
├── GenomeCache.sqlite
├── GFF/
└── Backups/
```

Despite the historical database name, Triple_A does **not** download genome FASTA files. `GenomeCache.sqlite` stores NCBI assembly accessions and computed gene results. When functional potential requires a new reference, Triple_A downloads only the corresponding GFF annotation into `Cache/GFF/`.

## Scientific scope and limitations

Triple_A supports exploratory and inferential analyses, but no software can replace experimental design or biological validation. Functional-potential results are predictions based on curated genes found in a representative reference annotation. They do not demonstrate gene expression, pathway activity or the exact genotype of the sampled strain. Ordination plots are exploratory unless accompanied by an appropriate statistical test. Statistical significance should always be interpreted together with effect size, sample size, multiple-testing correction and study design.

## Portability and reproducibility

The project is portable when its relative folder structure is retained. Do not hard-code personal paths inside analysis code. The cache can be exported as a ZIP containing `GenomeCache.sqlite`, `GFF/` and a manifest. Results should be archived together with the original input tables, parameter choices and the Triple_A release identifier.
