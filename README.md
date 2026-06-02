# nf-immcantation

A lightweight [Nextflow](https://www.nextflow.io/) DSL2 pipeline for bulk BCR-IGH AIRR-seq analysis using the [Immcantation](https://immcantation.readthedocs.io/) framework.

Designed for paired-end Illumina BCR heavy-chain (IGH) libraries where C-region primers are visible in R2 and V-gene primers are enzymatically removed before sequencing (R1 starts directly in VH FR1).

---

## Overview

```
fastp (QC + trim)
  └─ GUNZIP (decompress for pRESTO)
       └─ MaskPrimers align (R2 only — mask C-region primer)
            └─ PairSeq (synchronize R1/R2 pairs)
                 └─ AssemblePairs sequential + blastn (overlap assembly)
                      └─ CollapseSeq (deduplicate)
                           └─ SplitSeq (filter by duplicate count)
                                └─ AssignGenes / IgBLAST (V(D)J annotation)
                                     └─ MakeDb (AIRR-format TSV)
                                          └─ DefineClones (clonal grouping)
```

All IgBLAST databases and IMGT germlines are bundled in the `immcantation/suite:4.5.0` Docker image — no external reference downloads needed.

---

## Data

The pipeline was developed and validated on publicly available bulk BCR-IGH data from:

> Briney B, Inderbitzin A, Joyce C, Burton DR.  
> **Commonality despite exceptional diversity in the baseline human antibody repertoire.**  
> *Nature* 566, 393–397 (2019).  
> DOI: [10.1038/s41586-019-0879-y](https://doi.org/10.1038/s41586-019-0879-y)

SRA accessions used for testing (25,000-read toy subsets):

| Sample | SRA accession | Subject |
|--------|--------------|---------|
| SRR11909734 | [SRR11909734](https://www.ncbi.nlm.nih.gov/sra/SRR11909734) | DONOR1 |
| SRR11909736 | [SRR11909736](https://www.ncbi.nlm.nih.gov/sra/SRR11909736) | DONOR2 |

Raw FASTQs can be downloaded with the SRA Toolkit:

```bash
prefetch SRR11909734 SRR11909736
fasterq-dump --split-files SRR11909734
fasterq-dump --split-files SRR11909736
```

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
SRR11909734,DONOR1,/path/to/SRR11909734_R1.fastq.gz,/path/to/SRR11909734_R2.fastq.gz,human,IGH
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

#### Test run (toy Briney data, local resources)

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
| `--umi` | `false` | Opt-in UMI mode: extract the UMI as BARCODE and build a per-UMI consensus before assembly |
| `--buildconsensus_maxerror` | 0.1 | Max error within a UMI consensus group (UMI mode) |
| `--buildconsensus_mincount` | 1 | Min reads per UMI to build a consensus (UMI mode) |
| `--buildconsensus_maxgap` | 0.5 | Max gap fraction at a consensus position (UMI mode) |
| `--cloning_method` | `hierarchical` | Clonal grouping: `hierarchical` (SCOPer hierarchicalClones, default on this branch) or `exact` (DefineClones --model aa --dist 0, Briney parity) |
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
- **UMI consensus is opt-in** (`--umi`, default off): Briney's library carries a 12 nt UMI + 2/4/6 nt offset on R1, but the paper states most UMI bins are singletons (depth ≈ cells), so by default `BuildConsensus` is skipped in favour of exact-match collapse + DUPCOUNT≥2. With `--umi true` the C-read MaskPrimers step extracts the UMI into `BARCODE`, it is paired onto the V-read, each mate is consensus-built per UMI (`BuildConsensus`), the two consensus files are re-synced, and assembly proceeds as usual. Worth enabling only when UMIs have real bin depth and quantitative/error-sensitive accuracy matters. See `PLAN_umi_mode.md`.
- **Productive + allele-strip step (ParseDb)**: Briney's clonotype definition is `(V_gene, J_gene, CDR3_aa)` on productive sequences. After MakeDb we filter `productive=T` and write gene-level columns `v_call_gene` / `j_call_gene` for DefineClones to consume.
- **Clonal grouping** (`--cloning_method`):
  - `hierarchical` (default): SCOPer hierarchicalClones at `--clonal_threshold` (0.16), `--scoper_linkage single` (diagnosis sensitivity, Gupta et al. 2017).
  - `exact`: `DefineClones.py --model aa --dist 0 --vf v_call_gene --jf j_call_gene` — exact (V_gene, J_gene, CDR3_aa) match, matching Briney 2019. Used on the `briney` benchmark branch.
- **Grouping for clone definition**: samples from the same `params.cloneby` value (default `subject_id`) are clonotyped together, so clones span biological + technical replicates of one donor.
- **All references bundled**: IgBLAST DB and IMGT germlines come from `immcantation/suite:4.5.0`; no external download required.

---

## License

MIT
