# Plan: optional --umi mode (UMI extraction + BuildConsensus)

## Goal
Add an opt-in UMI workflow to PRESTO. Default (`--umi false`) is unchanged
(Briney exact-collapse). With `--umi true`, extract the UMI as a BARCODE,
build a per-UMI consensus on each mate, then assemble/collapse as usual.

## Read layout assumption (this repo / Briney 2019)
R1 = C-read (carries the 12 nt UMI + offset, 5' of the C primer), R2 = V-read.
UMI is extracted from the **C-read** by adding `--barcode` to MaskPrimers-C.
(Documented; other layouts would tune primers/which read carries the UMI.)

## Flow (umi = true)
```
FASTP → GUNZIP → split C-read / V-read
MaskPrimers-C(--barcode --pf C_CALL)   # extracts BARCODE (UMI) + isotype
MaskPrimers-V
PAIRSEQ_BARCODE  (-1 V, -2 C, --2f BARCODE C_CALL --coord sra)  # copy UMI+isotype onto V-read
BUILDCONSENSUS_V (-s V_pair, --bf BARCODE --cf C_CALL --act majority)
BUILDCONSENSUS_C (-s C_pair, --bf BARCODE --cf C_CALL --act majority)
PAIRSEQ_CONSENSUS (-1 Vcons, -2 Ccons, --coord presto)   # re-sync dropped UMIs
PEAR → CollapseSeq → SplitSeq
```
Default (umi = false) keeps the current single PAIRSEQ → PEAR path.

## New params (nextflow.config)
- `umi = false`
- `buildconsensus_maxerror = 0.1`
- `buildconsensus_mincount = 1`
- `buildconsensus_maxgap   = 0.5`

## Files
1. `modules/local/presto/buildconsensus/main.nf` — new `PRESTO_BUILDCONSENSUS`
   (label presto/process_medium; in (meta,reads), out consensus fastq + log + versions; stub).
2. `conf/modules.config`:
   - MaskPrimers-C args: append `${params.umi ? '--barcode' : ''}`.
   - add `withName: 'PAIRSEQ_BARCODE'`  (`--2f BARCODE C_CALL --coord sra`)
   - add `withName: 'PAIRSEQ_CONSENSUS'`(`--coord presto`)
   - add `withName: '.*BUILDCONSENSUS.*'`(`--bf BARCODE --cf C_CALL --act majority
     -n mincount --maxerror maxerror --maxgap maxgap`) + publishDir.
3. `subworkflows/local/presto.nf`:
   - include aliases: PRESTO_PAIRSEQ as PAIRSEQ_BARCODE / PAIRSEQ_CONSENSUS;
     PRESTO_BUILDCONSENSUS as BUILDCONSENSUS_V / BUILDCONSENSUS_C.
   - `if (params.umi)` branch implementing the flow above; else current path.
   - feed the chosen assembled-input channel into PEAR.

## Verify
- `-preview` umi=false: DAG identical to today (no BuildConsensus).
- `-preview` umi=true: shows PAIRSEQ_BARCODE, BUILDCONSENSUS_V/C, PAIRSEQ_CONSENSUS.
- `-stub-run` both modes.

## Non-goals
- Auto-detecting UMI presence/length (user sets `--umi`).
- Changing default behaviour or downstream VDJ/clonal steps.
- Related: [[project_buildconsensus_optional]].
