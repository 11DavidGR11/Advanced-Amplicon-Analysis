# Input data formats

## Abundance table

Accepted tabular formats include TSV, CSV, TXT, XLS and XLSX. TSV is recommended because it preserves plain text and avoids spreadsheet auto-formatting.

A valid table contains:

- one row per taxon, ASV or OTU;
- one or more taxonomy columns or a lineage column;
- one numeric abundance column per sample;
- unique and non-empty sample names;
- non-negative abundance values.

Example:

| FeatureID | Kingdom | Phylum | Genus | Sample_01 | Sample_02 |
|---|---|---|---|---:|---:|
| ASV_001 | Bacteria | Bacillota | Clostridium | 120 | 85 |
| ASV_002 | Bacteria | Bacteroidota | Prevotella | 42 | 73 |

Counts are required by methods that model sequencing counts. Relative abundances should not be presented as counts. Select the abundance scale accurately during import.

## Sample metadata

Use one row per sample. One column must contain sample identifiers matching the abundance-table sample names. Other columns may contain treatment, time, batch, subject, environmental measurements or other study variables.

Matching tolerates differences in letter case and in the separators `-`, `.`, `_` and spaces, so `S-1`, `S.1`, `S_1` and `s1` are treated as the same sample. Anything else must match exactly.

Every selected sample needs a row. A sample present in the abundance table but absent from the metadata blocks the design rather than being dropped silently.

Example:

| SampleID | Group | Time_day | Batch | pH |
|---|---|---:|---|---:|
| Sample_01 | Control | 0 | B1 | 7.2 |
| Sample_02 | Additive_A | 0 | B1 | 6.9 |

Assign each metadata column one role:

- **Sample ID:** the identifier used for matching. Chosen in the *Sample ID column* selector, not in the per-column role list;
- **Experimental factor:** categorical predictor such as treatment, time, batch or subject;
- **Environmental variable:** numeric covariate. RDA accepts only numeric predictors; for categorical ones use envfit, partial RDA or variance partitioning;
- **Ignore:** retained in the file but excluded from modelling.

Any column may additionally be selected as the **grouping column** that defines the experimental design, independently of the role assigned to it.

## Experimental design

The design states which sample belongs to which group. It can be declared in two ways.

**From a metadata column (recommended).** Each sample takes the group named in the selected column, matched by identifier. The order of the abundance columns and of the metadata rows is irrelevant to group membership, samples of the same group may be interleaved, and groups may contain **different numbers of replicates**.

**By column order, in consecutive blocks.** Samples are split into consecutive blocks of the declared replicate count, in the order the sample columns appear. This mode requires every group to have the **same** number of replicates, between one and five, and the columns to be already grouped by treatment. It cannot detect an incorrect order when the number of columns is right, so prefer the metadata column whenever a metadata file exists.

The order in which the groups appear is taken from the metadata file: groups are ordered by first appearance there. That order sets the reference level of every pairwise comparison, and therefore the sign of each log2 fold change, as well as the order of legends and heatmap columns. To control it, sort the metadata rows accordingly.

## FASTQ input

FASTQ files must use consistent sample naming. For paired-end data, forward and reverse files must form unambiguous pairs. Avoid spaces and unusual punctuation in filenames. Review the FASTQ workflow guide before running DADA2.

## Common causes of import failure

- sample identifiers differ by more than case and the separators `-`, `.`, `_` and spaces, which are the only differences tolerated;
- a selected sample has no metadata row, or an empty value in the grouping column;
- duplicated sample names or duplicated metadata rows;
- numeric columns contain text, commas used as decimal marks or missing-value labels;
- negative abundances;
- taxonomy and sample columns are assigned incorrectly;
- spreadsheet software has converted taxon names or identifiers.

Correct the source file rather than manually altering imported values whenever possible.
