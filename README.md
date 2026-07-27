# nf-immcantation (`main` branch)

A lightweight, computationally efficient [Nextflow](https://www.nextflow.io/) DSL2 pipeline for bulk BCR-IGH AIRR-seq analysis using the [Immcantation](https://immcantation.readthedocs.io/) framework.

Designed for paired-end Illumina BCR heavy-chain (IGH) libraries carrying a C-region (isotype) primer on one read and a V-region primer on the other. Ships two pipeline modes -- **bulletproof** (stock tools, validated baseline) and **boosted** (Rust-accelerated, ~55% faster) -- with an opt-in UMI-consensus path for libraries that carry molecular barcodes.

> **Reproducing a published study for benchmarking?** See the **[`briney`](../../tree/briney)** branch, which pins the exact `(V_gene, J_gene, CDR3_aa)` clonotype definition, ships the matching test data, and reports a gold-standard precision/recall benchmark. This branch is the general, configurable pipeline.

---

## Pipeline modes

Two profiles control which tool implementations are used. Compose them with infrastructure profiles (e.g. `-profile docker,local,boosted`):

| Profile | pRESTo steps | Read assembly | IgBLAST | Description |
|---------|-------------|---------------|---------|-------------|
| **`bulletproof`** (default) | Stock pRESTo | PEAR | Stock `AssignGenes` | Safe, validated baseline using canonical Immcantation tools |
| **`boosted`** | presto-fast (Rust) | AssemblePairsFast (Rust) | Stock `AssignGenes` | ~55% faster total compute time on Briney 2019 benchmarks |

Both modes share the same CHANGEO downstream steps (MakeDb, ParseDb, clonal analysis). Defaults are bulletproof -- no mode profile needed for the safe path.

> **Note:** An experimental k-mer pre-filtered IgBLAST (`--pre_igblast`) is available but off by default in both modes. It provides 1.2-1.5x speedup on small instances (<=4 CPUs) but is slower than stock IgBLAST on large instances (>=8 CPUs) where igblastn's native multi-threading is more efficient.

---

## Overview

```
fastp (QC + trim)
  +-- GUNZIP (decompress for pRESTo)
       +-- MaskPrimers align  (C-read -- isotype primer -> C_CALL)
       +-- MaskPrimers align  (V-read -- V-region primer)
            +-- PairSeq (synchronize the read pair; copy C_CALL onto the V-read)
                 +-- [ --umi: BuildConsensus per UMI barcode, then re-sync ]
                      +-- PEAR (bulletproof) / AssemblePairsFast (boosted)
                           +-- CollapseSeq (deduplicate)
                                +-- SplitSeq (filter by duplicate count, DUPCOUNT >= 2)
                                     +-- AssignGenes / IgBLAST (V(D)J annotation)
                                          +-- MakeDb (AIRR-format TSV)
                                               +-- ParseDb (productive filter + gene-level V/J calls)
                                                    +-- DefineClones -- SCOPer hierarchical (default) or exact
```

All IgBLAST databases and IMGT germlines are bundled in the Docker containers -- no external reference downloads needed.

---

## Requirements

| Tool | Version |
|------|---------|
| [Nextflow](https://www.nextflow.io/docs/latest/install.html) | >= 24.04.0 |
| [Docker](https://docs.docker.com/get-docker/) | any recent |

Containers are pulled automatically based on the pipeline mode.

---

## Quick start

### 1. Clone

```bash
git clone https://github.com/fmarchiano/nf-immcantation.git
cd nf-immcantation
```

### 2. Prepare a samplesheet

Edit `assets/samplesheet_template.csv` to point to your gzipped paired-end FASTQs:

```csv
sample,subject_id,fastq_1,fastq_2,species,pcr_target_locus
sampleA,subject1,/path/to/sampleA_R1.fastq.gz,/path/to/sampleA_R2.fastq.gz,human,IGH
```

Supported values: `species` = `human`; `pcr_target_locus` = `IGH`.

### 3. Prepare primer FASTAs

Two primer FASTAs are required:

- `--cprimers` -- isotype-specific C-region primers, masked on the C-read. Header IDs are propagated as the `C_CALL` (isotype) annotation.
- `--vprimers` -- V-region primers, masked on the V-read.

Example primer sets are provided under `assets/` to adapt to your own library.

### 4. Run

```bash
# Bulletproof (default -- stock tools, safe baseline)
nextflow run /path/to/nf-immcantation \
  -profile docker,local \
  --input /path/to/samplesheet.csv \
  --cprimers /path/to/cprimers.fasta \
  --vprimers /path/to/vprimers.fasta \
  --outdir /path/to/results

# Boosted (Rust-accelerated presto-fast + AssemblePairsFast)
nextflow run /path/to/nf-immcantation \
  -profile docker,local,boosted \
  --input /path/to/samplesheet.csv \
  --cprimers /path/to/cprimers.fasta \
  --vprimers /path/to/vprimers.fasta \
  --outdir /path/to/results
```

#### Test run (bundled toy data, local resources)

```bash
nextflow run /path/to/nf-immcantation \
  -profile docker,local,test \
  --outdir /path/to/results
```

`-profile test` runs a small bundled toy dataset with resource limits (`max_cpus=4`, `max_memory=6 GB`) suitable for a laptop or WSL2 environment.

---

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | -- | Samplesheet CSV (required) |
| `--cprimers` | -- | C-region (isotype) primers FASTA -- masked on the C-read (required) |
| `--vprimers` | -- | V-region primers FASTA -- masked on the V-read (required) |
| `--outdir` | -- | Output directory (required) |
| `--fastp_q` | 20 | Minimum base quality (phred) |
| `--fastp_window_size` | 5 | Sliding-window size for 3' trimming |
| `--primer_maxerror_c` | 0.3 | MaskPrimers max error rate (C-read) |
| `--primer_maxerror_v` | 0.2 | MaskPrimers max error rate (V-read) |
| `--primer_maxlen_c` | 100 | Search window on the C-read (covers a barcode + offset preamble) |
| `--primer_maxlen_v` | 35 | Search window on the V-read (covers the offset preamble) |
| `--splitseq_min_count` | 2 | Minimum duplicate count to retain a sequence |
| `--skip_clonal` | `true` | Skip clonal analysis (ParseDb is the last step) |
| `--umi` | `false` | Opt-in UMI mode: extract the UMI as `BARCODE` and build a per-UMI consensus before assembly |
| `--buildconsensus_maxerror` | 0.1 | Max error within a UMI consensus group (UMI mode) |
| `--buildconsensus_mincount` | 1 | Min reads per UMI to build a consensus (UMI mode) |
| `--buildconsensus_maxgap` | 0.5 | Max gap fraction at a consensus position (UMI mode) |
| `--cloning_method` | `hierarchical` | Clonal grouping: `hierarchical` (SCOPer hierarchicalClones, default) or `exact` (DefineClones `--model aa --dist 0`) |
| `--defineclones_model` | `aa` | DefineClones distance model when `cloning_method=exact` |
| `--defineclones_dist` | 0.0 | DefineClones distance threshold (0 = exact CDR3 aa match) |
| `--clonal_threshold` | 0.16 | SCOPer junction distance cutoff (only when `cloning_method=hierarchical`) |

---

## Output structure

```
results/
+-- fastp/{sample}/                       # QC reports + trimmed reads
+-- presto/
|   +-- 02a-maskprimers-C/{sample}/       # MaskPrimers on the C-read (isotype primers)
|   +-- 02b-maskprimers-V/{sample}/       # MaskPrimers on the V-read (VH primers)
|   +-- 03-pairseq/{sample}/
|   +-- 04-pear/{sample}/                 # PEAR-assembled reads (bulletproof)
|   +-- 04-assemblepairs/{sample}/        # AssemblePairsFast output (boosted)
|   +-- 05-collapseseq/{sample}/
|   +-- 06-splitseq/{sample}/
+-- vdj_annotation/
|   +-- 01-assigngenes/{sample}/          # IgBLAST output (.fmt7)
|   +-- 02-makedb/{sample}/               # AIRR-format TSV (db-pass.tsv)
|   +-- 03-parsedb/{sample}/              # productive-only + v_call_gene/j_call_gene added
+-- clonal_analysis/{subject}/            # clone-pass.tsv with clone_id column (if skip_clonal=false)
+-- pipeline_info/                        # Nextflow execution report, timeline, trace
```

---

## Design notes

- **Read orientation**: the C-read carries the isotype/C-region primer and the V-read carries the V-region primer. Some library layouts place the C-read on R1 and the V-read on R2 -- inverted from the common R1 = V convention -- so confirm the orientation against primer hits on your raw FASTQs.
- **Two MaskPrimers steps**:
  - C-read: isotype primers (IgM/IgG/...), preceded by a barcode + 2/4/6 nt offset preamble, `--maxlen 100`. The isotype name is written to `C_CALL` for downstream isotype assignment.
  - V-read: VH primers with a 2/4/6 nt offset, `--maxlen 35`. Removed pre-assembly.
- **UMI consensus is opt-in** (`--umi`, default off): by default reads are deduplicated by exact-match collapse + `DUPCOUNT >= 2`, which is fastest and works well when UMI bins are mostly singletons. With `--umi true` the C-read MaskPrimers step extracts the UMI into `BARCODE`, it is paired onto the V-read, each mate is consensus-built per UMI (`BuildConsensus`), the two consensus files are re-synced, and assembly proceeds as usual. Worth enabling only when UMIs have real bin depth and quantitative/error-sensitive accuracy matters. See `PLAN_umi_mode.md`.
- **Productive + allele-strip step (ParseDb)**: the clonotype definition is `(V_gene, J_gene, CDR3_aa)` on productive sequences. After MakeDb we filter `productive=T` and write gene-level columns `v_call_gene` / `j_call_gene` for DefineClones to consume.
- **Clonal grouping** (`--cloning_method`):
  - `hierarchical` (default): SCOPer hierarchicalClones at `--clonal_threshold` (0.16), single linkage (diagnosis sensitivity, Gupta et al. 2017).
  - `exact`: `DefineClones.py --model aa --dist 0 --vf v_call_gene --jf j_call_gene` -- exact `(V_gene, J_gene, CDR3_aa)` match. This is the mode used for the published-study benchmark on the [`briney`](../../tree/briney) branch.
- **Grouping for clone definition**: samples from the same `params.cloneby` value (default `subject_id`) are clonotyped together, so clones span biological + technical replicates of one subject.
- **All references bundled**: IgBLAST DB and IMGT germlines come from the Docker containers; no external download required.

---

## License

MIT
