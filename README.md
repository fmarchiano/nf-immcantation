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

### 3. Prepare a C-region primer FASTA

Provide a FASTA file of C-region primer sequences visible in R2.  
The Briney 2019 CH1 primer (`assets/cprimers_briney2019.fasta`) is included as an example.

### 4. Run

```bash
nextflow run /path/to/nf-immcantation \
  -profile docker,local \
  --input /path/to/samplesheet.csv \
  --cprimers /path/to/cprimers.fasta \
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
| `--cprimers` | — | C-region primers FASTA (required) |
| `--outdir` | — | Output directory (required) |
| `--fastp_q` | 20 | Minimum base quality (phred) |
| `--fastp_window_size` | 5 | Sliding-window size for 3′ trimming |
| `--primer_maxerror` | 0.3 | MaskPrimers max error rate |
| `--ap_maxerror` | 0.3 | AssemblePairs max error rate |
| `--ap_alpha` | 1e-5 | AssemblePairs p-value threshold |
| `--ap_minident` | 0.5 | AssemblePairs min identity |
| `--splitseq_min_count` | 2 | Minimum duplicate count to retain a sequence |
| `--clonal_threshold` | 0.16 | Fixed Hamming distance threshold for clonal grouping |

---

## Output structure

```
results/
├── fastp/{sample}/            # QC reports + trimmed reads
├── presto/
│   ├── 02-maskprimers/{sample}/
│   ├── 03-pairseq/{sample}/
│   ├── 04-assemblepairs/{sample}/
│   ├── 05-collapseseq/{sample}/
│   └── 06-splitseq/{sample}/
├── vdj_annotation/
│   ├── 01-assigngenes/{sample}/   # IgBLAST output (.fmt7)
│   └── 02-makedb/{sample}/        # AIRR-format TSV (db-pass.tsv)
├── clonal_analysis/{subject}/     # clone-pass.tsv with clone_id column
└── pipeline_info/                 # Nextflow execution report, timeline, trace
```

---

## Design notes

- **No UMI consensus step**: Briney 2019 uses variable-length staggers (2/4/6 bp), not fixed-length UMIs. Most bins contain 1 read, so `BuildConsensus` is skipped.
- **No MaskPrimers on R1**: V-gene primers are enzymatically removed before sequencing; R1 starts directly in VH FR1.
- **Fixed clonal threshold**: `--clonal_threshold 0.16` (Hamming distance, length-normalized) — no Shazam auto-detection.
- **DefineClones groups by `subject_id`**: samples from the same subject are clonotyped together.
- **All references bundled**: IgBLAST DB and IMGT germlines come from `immcantation/suite:4.5.0`; no external download required.

---

## License

MIT
