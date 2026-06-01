# Exact vs Hierarchical Clonal Assignment — Briney Subset

Comparison of two `cloning_method` settings on the same 30-sample Briney subset (10 subjects × 3 biological replicates).

## Run setup

| | Exact | Hierarchical |
|---|---|---|
| `--cloning_method` | `exact` | `hierarchical` |
| Tool | Change-O `DefineClones.py` | SCOPer `hierarchicalClones` |
| Parameters | `--model aa --dist 0 --norm len` | `method=nt`, `linkage=complete`, `threshold=0.16` |
| Clone grouping | `subject_id` | `subject_id` |
| Outdir | `/data/scratch/subset_briney/exact` | `/data/scratch/subset_briney/hierarc` |
| Wallclock | 2h 34m 39s | 9m 4s (with `-resume`, 330/340 tasks cached) |

`exact` mirrors Briney's clonotype definition (V_gene + J_gene + exact CDR3 amino acid sequence). `hierarchical` merges junctions within hamming-distance 0.16 (length-normalised) under complete linkage.

## Clone counts per subject

| Subject  | Sequences | Exact clones | Hierarchical clones | h/e ratio |
|----------|----------:|-------------:|--------------------:|----------:|
| D103     |    31,788 |       29,596 |              27,934 |     0.944 |
| S316188  |    41,123 |       33,973 |              29,740 |     0.875 |
| S326650  |    17,199 |       16,158 |              15,419 |     0.954 |
| S326651  |     5,201 |        3,326 |               2,966 |     0.892 |
| S326713  |     3,373 |        2,698 |               2,544 |     0.943 |
| S326737  |    27,944 |       26,403 |              25,101 |     0.951 |
| S326780  |    14,094 |       11,366 |              10,600 |     0.933 |
| S326797  |    17,657 |       13,461 |              12,443 |     0.924 |
| S326907  |    13,407 |       12,832 |              12,216 |     0.952 |
| S327059  |    10,965 |        8,288 |               7,473 |     0.902 |
| **TOTAL**| **182,751** | **158,101** |          **146,436** | **0.926** |

## Observations

- Hierarchical consistently produces fewer clones (≈7% fewer overall). Expected: SCOPer merges junction-similar sequences at the 0.16 hamming threshold, while `exact --dist 0` keeps them as separate clonotypes.
- Strongest compression: **S316188** (12.5% fewer), **S327059** (10% fewer) — subjects with the largest gap between sequence count and clone count, indicating closely-related sequence families that hierarchical groups together.
- Weakest compression: **S326650**, **S326737**, **S326907** (~5%) — sequences-to-clones ratio is already near 1, so little room for SCOPer to merge further.
- Sequence-to-clone ratio across all subjects: 1.156 (exact), 1.248 (hierarchical) — both modest, consistent with the singleton-dominated nature of this subsampled dataset.

## Reproduction

```bash
# Exact run
nextflow run main.nf -profile docker,test \
  --outdir /data/scratch/subset_briney/exact \
  --cloning_method exact

# Hierarchical run (reuses cached upstream)
nextflow run main.nf -profile docker,test \
  --outdir /data/scratch/subset_briney/hierarc \
  --cloning_method hierarchical \
  -resume
```

Clone-count tally (run from `/home/fabio/nf-immcantation/assets`):

```bash
for subj in D103 S316188 S326650 S326651 S326713 S326737 S326780 S326797 S326907 S327059; do
  ef="/data/scratch/subset_briney/exact/clonal_analysis/$subj/${subj}_clone-pass.tsv"
  hf="/data/scratch/subset_briney/hierarc/clonal_analysis/$subj/${subj}_clone-pass.tsv"
  e_col=$(head -1 "$ef" | tr '\t' '\n' | grep -n '^clone_id$' | cut -d: -f1)
  h_col=$(head -1 "$hf" | tr '\t' '\n' | grep -n '^clone_id$' | cut -d: -f1)
  e_rows=$(($(wc -l < "$ef") - 1))
  e_clones=$(tail -n +2 "$ef" | cut -f"$e_col" | sort -u | wc -l)
  h_clones=$(tail -n +2 "$hf" | cut -f"$h_col" | sort -u | wc -l)
  printf "%-8s  seqs=%-7d  exact=%-7d  hierarc=%-7d  ratio=%.3f\n" \
    "$subj" "$e_rows" "$e_clones" "$h_clones" \
    "$(awk -v h=$h_clones -v e=$e_clones 'BEGIN{print h/e}')"
done
```
