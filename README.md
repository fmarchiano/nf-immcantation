# nf-immcantation — Briney 2019 benchmark (`briney` branch)

> **This is the benchmark branch.** It reproduces the bulk BCR-IGH processing of
> Briney et al. 2019 as faithfully as possible, so the pipeline's output can be
> compared against the published analysis. For the general-purpose, more
> computationally efficient pipeline, see the **[`main`](../../tree/main)** branch.

A [Nextflow](https://www.nextflow.io/) DSL2 pipeline for bulk BCR-IGH AIRR-seq
analysis using the [Immcantation](https://immcantation.readthedocs.io/)
framework, configured to mirror the Briney 2019 library design and clonotype
definition.

Targets the Briney 2019 SRA library layout: paired-end Illumina BCR heavy-chain
(IGH) reads where **R1 is the C-read** (12 nt UMI + 2/4/6 nt offset + isotype
C-region primer) and **R2 is the V-read** (VH1–VH6 primer) — inverted from the
common R1 = V convention. Both primer sets are present and are masked before
assembly.

---

## Overview

```
fastp (QC + trim)
  └─ GUNZIP (decompress for pRESTO)
       ├─ MaskPrimers align  (R1 / C-read — isotype primer → C_CALL)
       └─ MaskPrimers align  (R2 / V-read — VH1–VH6 primer)
            └─ PairSeq (synchronize R1/R2 pairs; copy C_CALL onto the V-read)
                 └─ PEAR (overlap assembly of the mate pair)
                      └─ CollapseSeq (deduplicate)
                           └─ SplitSeq (filter by duplicate count, DUPCOUNT ≥ 2)
                                └─ AssignGenes / IgBLAST (V(D)J annotation)
                                     └─ MakeDb (AIRR-format TSV)
                                          └─ ParseDb (productive filter + gene-level V/J calls)
                                               └─ DefineClones (exact V_gene, J_gene, CDR3_aa)
```

All IgBLAST databases and IMGT germlines are bundled in the `immcantation/suite:4.5.0` Docker image — no external reference downloads needed.

---

## Data

The pipeline was developed and validated on publicly available bulk BCR-IGH data from:

> Briney B, Inderbitzin A, Joyce C, Burton DR.  
> **Commonality despite exceptional diversity in the baseline human antibody repertoire.**  
> *Nature* 566, 393–397 (2019).  
> DOI: [10.1038/s41586-019-0879-y](https://doi.org/10.1038/s41586-019-0879-y)

The benchmark uses the three subjects with published gold-standard annotations —
**316188, 326650, 326651** — each sequenced as **18 libraries (6 biological
replicates × 3 technical replicates)**. Read depth spans a ~49× range across
subjects (≈ 0.66 M to ≈ 32 M sequences), so the benchmark probes behaviour from
shallow to near-saturating coverage. Raw FASTQs are available from the study's
SRA BioProject; the exact runs per subject and download steps are listed in the
project wiki (`Benchmark/Briney_Gold_Standard.md`).

---

## Benchmark results

This branch's exact `(V_gene, J_gene, CDR3_aa)` clonotyping was benchmarked
against Briney's own published annotations (the *gold standard*) for the three
subjects with released data — **316188, 326650, 326651**. Each gold record is a
UMI-consensus sequence annotated by AbStar and projected to the same exact
`(V_gene, J_gene, CDR3_aa)` clonotype this branch produces, so the comparison is
definition-for-definition.

**Method** (asymmetric, `MIN_SUP=2`): recall is measured against *confident*
gold (UMI support ≥ 2, the analog of the pipeline's `DUPCOUNT ≥ 2` filter), and
precision against *full* gold (≥ 1), since a call matching **any** real gold
molecule is a true positive. A match requires identical CDR3_aa **and**
overlapping V- and J-call lists (tolerating IgBLAST-vs-AbStar call ordering).
Biological replicates are pooled within a subject; subjects are weighted equally
in the macro-average.

| Subject | Gold ≥2 units | Pipeline units | Recall | Precision | F1 |
|---------|--------------:|---------------:|-------:|----------:|-----:|
| 316188 | 56,025 | 102,022 | 0.954 | 0.778 | 0.857 |
| 326650 | 633,674 | 1,259,356 | 0.949 | 0.890 | 0.919 |
| 326651 | 831,688 | 2,233,249 | 0.916 | 0.983 | 0.948 |
| **Macro-average** | — | — | **0.939** | **0.884** | **0.911** |

**Recall ≈ 0.94, precision ≈ 0.88, F1 ≈ 0.91.** The pipeline recovers ~94% of
confidently-supported gold clonotypes, and ~88% of its own calls are confirmed
against gold. Precision rises monotonically with gold depth (0.78 → 0.89 → 0.98):
a shallower gold set leaves more pipeline calls *unconfirmable* rather than
wrong, so 0.88 is a conservative floor — against the near-saturated deepest
subject (326651) precision reaches ~0.98.

Gold-standard construction, full methodology, and per-replicate detail are
documented in the project wiki (`Benchmark/Briney_Gold_Standard.md`,
`Benchmark/Comparison_Methodology.md`, `Benchmark/Annotation_vs_Gold_Results.md`).

---

## Requirements

| Tool | Version |
|------|---------|
| [Nextflow](https://www.nextflow.io/docs/latest/install.html) | ≥ 24.04.0 |
| [Docker](https://docs.docker.com/get-docker/) | any recent |
| Container: `immcantation/suite` | 4.5.0 (auto-pulled) |

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
SRR8283834,316188,/path/to/SRR8283834_R1.fastq.gz,/path/to/SRR8283834_R2.fastq.gz,human,IGH
```

Supported values: `species` = `human`; `pcr_target_locus` = `IGH`.

### 3. Prepare primer FASTAs

Two primer FASTAs are required:

- `--cprimers` — isotype-specific C-region primers, masked on the C-read (reads[0] in Briney 2019 SRA dump).
  Header IDs are propagated as the `C_CALL` (isotype) annotation.
  `assets/cprimers_briney2019.fasta` ships IgM + IgG (degenerate `S`/`R` codes cover IgG1–4).
- `--vprimers` — V-region primers, masked on the V-read (reads[1] in Briney 2019 SRA dump).
  `assets/vprimers_briney2019.fasta` ships VH1–VH6.

### 4. Run

```bash
nextflow run /path/to/nf-immcantation \
  -profile docker,local \
  --input /path/to/samplesheet.csv \
  --cprimers /path/to/cprimers.fasta \
  --vprimers /path/to/vprimers.fasta \
  --outdir /path/to/results
```

#### Quick smoke test (bundled data, local resources)

```bash
nextflow run /path/to/nf-immcantation \
  -profile docker,local,test \
  --outdir /path/to/results
```

`-profile test` sets resource limits (`max_cpus=4`, `max_memory=6 GB`) suitable for a laptop or WSL2 environment.

---

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | Samplesheet CSV (required) |
| `--cprimers` | — | C-region (isotype) primers FASTA — masked on R1 / C-read (required) |
| `--vprimers` | — | V-region primers FASTA — masked on R2 / V-read (required) |
| `--outdir` | — | Output directory (required) |
| `--fastp_q` | 20 | Minimum base quality (phred) |
| `--fastp_window_size` | 5 | Sliding-window size for 3′ trimming |
| `--primer_maxerror_c` | 0.3 | MaskPrimers max error rate (C-read) |
| `--primer_maxerror_v` | 0.2 | MaskPrimers max error rate (V-read) |
| `--primer_maxlen_c` | 100 | Search window on R1 (covers 12 nt UMI + 2/4/6 nt offset) |
| `--primer_maxlen_v` | 35 | Search window on R2 (covers 2/4/6 nt offset) |
| `--splitseq_min_count` | 2 | Minimum duplicate count to retain a sequence |
| `--cloning_method` | `exact` | Clonal grouping: `exact` (DefineClones --model aa --dist 0, Briney parity) or `hierarchical` (SCOPer) |
| `--defineclones_model` | `aa` | DefineClones distance model when `cloning_method=exact` |
| `--defineclones_dist` | 0.0 | DefineClones distance threshold (0 = exact CDR3 aa match) |
| `--clonal_threshold` | 0.16 | SCOPer junction distance cutoff (only when `cloning_method=hierarchical`) |

---

## Output structure

```
results/
├── fastp/{sample}/                       # QC reports + trimmed reads
├── presto/
│   ├── 02a-maskprimers-C/{sample}/       # MaskPrimers on R1 (C-read, isotype primers)
│   ├── 02b-maskprimers-V/{sample}/       # MaskPrimers on R2 (V-read, VH primers)
│   ├── 03-pairseq/{sample}/
│   ├── 04-pear/{sample}/                 # PEAR-assembled reads
│   ├── 05-collapseseq/{sample}/
│   └── 06-splitseq/{sample}/
├── vdj_annotation/
│   ├── 01-assigngenes/{sample}/          # IgBLAST output (.fmt7)
│   ├── 02-makedb/{sample}/               # AIRR-format TSV (db-pass.tsv)
│   └── 03-parsedb/{sample}/              # productive-only + v_call_gene/j_call_gene added
├── clonal_analysis/{subject}/            # clone-pass.tsv with clone_id column
└── pipeline_info/                        # Nextflow execution report, timeline, trace
```

---

## Design notes

- **Read orientation (Briney 2019 SRA)**: R1 = C-read, R2 = V-read — inverted from the common R1=V convention. Confirmed by primer hits on raw FASTQs.
- **Two MaskPrimers steps**:
  - R1 (C-read): isotype primers (IgM/IgG/...), 12 nt UMI + 2/4/6 nt offset preamble, `--maxlen 100`. Isotype name written to `C_CALL` header for downstream isotype assignment.
  - R2 (V-read): VH1–VH6 primers, 2/4/6 nt offset, `--maxlen 35`. Briney removed these with cutadapt post-assembly; we remove pre-assembly for the same effect.
- **No UMI consensus step**: Briney's library carries a 12 nt UMI + 2/4/6 nt offset on R1. The paper itself states most UMI bins are singletons (depth ≈ cells, ~3×10⁸ each), so `BuildConsensus` adds little benefit. Replaced by exact-match collapse + DUPCOUNT≥2.
- **Productive + allele-strip step (ParseDb)**: Briney's clonotype definition is `(V_gene, J_gene, CDR3_aa)` on productive sequences. After MakeDb we filter `productive=T` and write gene-level columns `v_call_gene` / `j_call_gene` for DefineClones to consume.
- **Clonal grouping** (`--cloning_method`):
  - `exact` (default): `DefineClones.py --model aa --dist 0 --vf v_call_gene --jf j_call_gene` — exact (V_gene, J_gene, CDR3_aa) match, matching Briney 2019.
  - `hierarchical`: SCOPer hierarchicalClones at `--clonal_threshold` (the previous default behaviour).
- **Grouping for clone definition**: samples from the same `params.cloneby` value (default `subject_id`) are clonotyped together, so clones span biological + technical replicates of one donor.
- **All references bundled**: IgBLAST DB and IMGT germlines come from `immcantation/suite:4.5.0`; no external download required.

---

## License

MIT
