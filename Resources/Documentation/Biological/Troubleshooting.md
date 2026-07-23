# Troubleshooting

## The launcher does not open

Run `Run_Triple_A.R` from the extracted project and review the package message shown in the R console. Do not change the internal folder structure. Use `Tools/Diagnostics/Verify_Installation.R` to identify missing directories or dependencies.

## An application opens and closes immediately

Check its console output and `Results/Logs/`. A package may be missing, a port may be unavailable or the project root may not have been detected. Start from the launcher rather than the internal `app.R`.

## No usable taxa remain after filtering

Check that sample columns are numeric, sample names match the design table, taxonomy is in a supported format, abundance values are not all zero and filters are not too strict. Inspect duplicated taxonomic rows and unintended text columns. Test with filtering relaxed, then reintroduce justified thresholds.

## Sample names do not match

Differences in letter case and in the separators `-`, `.`, `_` and spaces are tolerated, so `S-1` and `s_1` are the same sample; anything else must match. Each selected sample column needs a metadata row when the design is read from the metadata. Do not silently rename samples after an analysis has started.

## MultiAmplicon export downloads HTML instead of TSV

Use the module's TSV download control after a successful merge. If a browser error page is saved, inspect the R console for a download-handler exception and verify that the merged reactive table is available. Do not rename an HTML error response to `.tsv`.

## Functional-potential results are empty

Possible causes include unresolved taxonomy, no suitable species/genus reference, unavailable NCBI service, missing GFF, gene names not present in the annotation, or a biological rule not met. Review the taxon-level reference and gene-evidence tables before interpreting the summary.

## Does Triple_A download genomes?

No. It stores an NCBI assembly accession in SQLite and downloads only the assembly GFF annotation. There should be no `Cache/Genomes/` directory in current releases.

## GFF download failed

Confirm internet access and NCBI availability. Check whether a stale `.gff.lock` remains after an interrupted run. Close all Triple_A processes before removing a stale lock. Do not delete a lock while another process may still be downloading the same annotation.

## Cache verification reports missing files

The SQLite index contains an accession without a corresponding local GFF. Re-run functional potential with network access or import a portable cache containing the required GFF. A SQLite-only import cannot restore annotation files.

## Cache verification reports orphan files

An orphan GFF is not linked to an indexed accession. Preview cleanup first. A backup is created before deletion, but retain it until the application has been tested.

## Cache import reports conflicts

The same GFF filename exists with different content, or the same taxonomy has different timestamps. Preserve existing records unless the imported cache is known to be authoritative. Review the operation history and backup before choosing replacement.

## Cache appears not to be reused

Run the same taxon/function twice and inspect the operation messages, GFF modification time and SQLite records. A cache hit should avoid downloading a new GFF and should not change the existing file time. Reuse can fail when taxonomy strings differ, the accession changed, the GFF is missing/corrupt, or requested genes were not previously evaluated.

## FASTQ paired files do not match

Ensure each R1 has exactly one R2, names follow a consistent convention and upload order is corresponding. Do not pair files only by alphabetical position when names are ambiguous.

## Too many reads are lost during filtering

Inspect quality profiles, reduce truncation only when necessary, verify primer/adaptor removal and relax expected-error thresholds cautiously. For paired reads, retain sufficient overlap. Record all modified parameters.

## Paired reads fail to merge

The truncated forward and reverse reads may not overlap sufficiently, or the amplicon may be longer than expected. Recalculate expected overlap from amplicon length and truncation settings. Do not force merging with implausibly small overlap.

## Taxonomy assignment is poor

Verify that the training FASTA matches the marker region and taxonomy convention. Poor read quality, short retained sequences and an unsuitable reference database all reduce assignment quality.

## An analysis is slow

Large numbers of taxa, permutations, cross-validation and network retrieval can be expensive. First verify that the run is progressing. Reusing the cache speeds repeated functional analyses. Reduce parameters only when scientifically acceptable, not solely to obtain a faster plot.

## Functional potential is slow on a first (cold-cache) run

NCBI limits unauthenticated requests to 3 per second, which functional potential respects with a short delay before each taxonomy/assembly lookup. A free NCBI API key (from the user's NCBI account "Settings" page) raises that limit to 10 per second and roughly triples throughput on new taxa. Set it before launching Triple_A with `options(triple_a_ncbi_api_key = "your-key-here")` or by setting the `NCBI_API_KEY` environment variable; already-cached taxa are unaffected either way.

## A result folder already contains older outputs

Use a new run name. Old files can be mistaken for current outputs even when the current analysis failed. Archive or remove obsolete result directories outside the active run; do not overwrite evidence needed for reproducibility.

## A QA test still expects Cache/Genomes

The test belongs to an older cache design. Current architecture requires `Cache/GenomeCache.sqlite` and `Cache/GFF/`. Update the test rather than recreating an unused genome directory.
