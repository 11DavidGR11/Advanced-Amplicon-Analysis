# Complete functional verification

## Launcher

- Start Triple_A through `Run_Triple_A.R`.
- Open Biological Analysis, FASTQ Pipeline, MultiAmplicon Integrator, Cache Manager and Function Builder.
- Confirm each child application opens and writes a readable startup log.

## Biological Analysis

- Import a valid abundance table.
- Verify sample selection, original sample names and grouping metadata.
- Confirm exploratory analyses run.
- Confirm PLS-DA is unavailable without biological replicates.
- Save and reopen a project from History.

## FASTQ

- Test valid paired-end and single-end datasets.
- Confirm orphan, malformed and truncated files produce clear errors.
- Review read-retention summaries and exported abundance tables.

## Reference Genome Cache

- Export a portable cache ZIP.
- Import it into a temporary empty cache.
- Confirm SQLite, genome and GFF content are transferred together.
- Merge the same package twice and confirm file counts do not increase on the second merge.
- Verify integrity and preview orphan cleanup.

## Documentation

- Open contextual help in every application and review each discovered topic.
- Confirm guides render inside the interface and all referenced paths exist.

## Distribution hygiene

- Confirm runtime results and personal cache content are not bundled.
- Confirm no duplicate cache directories or obsolete documentation files exist.
