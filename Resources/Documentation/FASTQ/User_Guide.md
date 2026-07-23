# FASTQ pipeline guide

The FASTQ application performs technical preprocessing and ASV inference with DADA2. Select single-end or paired-end layout, provide the corresponding reads and a taxonomy training FASTA, review truncation lengths and expected-error thresholds, and run the workflow. For paired reads, retain sufficient overlap after truncation. Inspect read tracking before using the exported abundance table. fastp and cutadapt are optional external preprocessing tools. The FASTQ module should be incorporated last in a project workflow because its exported table becomes the input for downstream integration or biological analysis.

## Preparing and uploading your data

### How many files can I upload?

| Field | How many |
|---|---|
| **FASTQ / forward (R1) files** | **Several** — one per sample |
| **Reverse (R2) files** | **Several** — one per sample (shown only in paired-end mode) |
| **Taxonomy training FASTA** | **Exactly one** (required) |

There is no limit on the *number* of files. The real limit is the **total upload size: 10 GB**.

### How do I upload several files? Multiple selection — not a ZIP

Press **Browse…** under *FASTQ / forward (R1) files* and select **all your R1 files at once** (Ctrl+click for individual files, Shift+click for a range). Then repeat under *Reverse (R2) files* with all the R2 files.

- The module **does not accept ZIP archives**: it does not unpack anything, the files you upload are passed straight to DADA2.
- It **does** accept each FASTQ compressed individually as `.fastq.gz` — this is the normal case and DADA2 reads gzipped files natively, so do not decompress them.
- The number of R1 files must match the number of R2 files, and they must pair by sample name (see below).

### Which file should I use when the provider sends several?

Sequencing providers often deliver three variants per sample, for example:

```
SampleX.raw_1.fastq.gz / SampleX.raw_2.fastq.gz    raw reads
SampleX_1.fastq.gz     / SampleX_2.fastq.gz        cleaned reads
SampleX.extendedFrags.fastq.gz                     already merged (FLASH)
```

- **Use the paired, non-merged reads** — either the cleaned `_1/_2` pair or the `raw_1/raw_2` pair — in **paired-end** mode.
- **Do not use a pre-merged file** (`extendedFrags` or similar). DADA2 must perform the merge itself, *after* denoising R1 and R2 separately; feeding it reads that are already merged invalidates its error model and defeats the paired-end workflow.
- Choose `raw_*` if you prefer to do all preprocessing inside Triple_A; choose the cleaned pair if you trust the provider's adapter/quality filtering.

**Adapters are not primers.** "Clean" data from a provider usually has sequencing *adapters* removed, but amplicon *primers* may still be present. If yours still carry primers, either enable **cutadapt** under *Optional external preprocessing* or remove them with truncation. Primers left in place distort DADA2's error learning and taxonomic assignment.

### File naming and pairing

Sample names are derived from the filename by removing the extension and then **one** read-direction suffix: a separator (`.`, `_` or `-`), an optional `R`, the direction digit, and an optional lane block. All of these work:

```
SampleX_1.fastq.gz   SampleX.1.fastq.gz    SampleX-1.fastq.gz
SampleX_R1.fastq.gz  SampleX_S1_L001_R1_001.fastq.gz
```

Requirements:

- Forward and reverse files of the same sample must produce the **same** name once the suffix is removed (`SampleX_1` and `SampleX_2` → `SampleX`).
- Sample names must be **unique** after that removal; otherwise the run stops with *"Paired-end sample names are ambiguous"*.
- Avoid spaces and unusual punctuation in filenames.
- Reverse files do not need to be in the same order as the forward files — they are matched by name.

Sample names that themselves end in a digit (for example `D90.4.1_1.fastq.gz` → sample `D90.4.1`) are handled correctly.

### After the run

Check the **Read tracking** tab before using the results: it shows how many reads survive filtering, denoising and chimera removal at each step. A large drop usually means the truncation lengths are too aggressive or the paired reads no longer overlap. The exported table for downstream analysis is written to `07_Export/Triple_A_input_table.csv`.
