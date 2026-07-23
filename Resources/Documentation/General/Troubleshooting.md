# General troubleshooting

## The launcher does not open

Run `source("Tools/Diagnostics/Verify_Installation.R")`. Install any missing launcher packages. Confirm that R is allowed to open a local browser and that security software is not blocking localhost connections.

## A child application does not start

Check `Results/Logs/` for the most recent startup log. The usual causes are a missing R package, an incomplete folder copy or a file path containing inaccessible network locations.

## An analysis is disabled

A locked analysis states what is missing and, under **How to unlock**, the concrete steps that enable it. Read that notice first: it is generated from the same rules the analysis engine enforces, so following it is sufficient.

The usual causes are missing metadata, an incompatible abundance scale, too few samples, too few groups, groups with a single sample, or a missing optional package. A locked analysis is not a defect: it means the current design cannot support that method.

The same check runs when the analysis is launched, so a stale selection cannot bypass it.

## Sample names do not match

Compare the abundance-table column names with the metadata identifier column. Differences in letter case and in the separators `-`, `.`, `_` and spaces are tolerated, so `S-1` and `s_1` match. Anything else must be identical.

Check also that every selected sample has a metadata row: a sample missing from the metadata blocks the design instead of being dropped. Either add the row or deselect that sample column.

## The grouping is wrong or the groups are unexpected

Read the design summary shown above the analysis panels, for example `3 group(s): Control (n=5), Treated (n=7)`. It states how the experiment was interpreted and should be checked before every run.

If the design is taken from the column order, samples are assigned to treatments in consecutive blocks, so an interleaved column order produces a complete analysis with every treatment mislabelled. Use a metadata grouping column instead, which matches by identifier and is insensitive to column order.

## The analysis finishes without expected taxa

Review prevalence and abundance filtering, taxonomy parsing and the selected taxonomic level. Highly strict filters or incomplete taxonomy can leave no usable features.

## Functional annotation is unavailable

Confirm internet access, inspect the cache with **Cache Manager**, and verify that the organism was classified to species or genus. Triple_A first seeks species-level evidence and may fall back to genus-level evidence when species resolution is unavailable; it does not infer above genus for this lookup.

## A download opens as HTML

Do not open the download link in a text editor. Save the file using the browser download action and verify its extension. Tabular exports should be `.tsv`, while complete result bundles are `.zip`.

## Reporting a reproducible problem

Keep the input files, the relevant startup or analysis log, the exact analysis settings and the Triple_A version. Remove confidential sample information before sharing files.
