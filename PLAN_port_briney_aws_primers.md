# Plan: SCOPer branch = briney, but cluster with SCOPer

## Insight
`briney` is a SUPERSET: `clonal_analysis.nf` already switches on
`params.cloning_method` (`exact`=DefineClones, `hierarchical`=SCOPer), already
ships `SCOPER_HIERARCHICALCLONES`, and `bin/run_scoper.R` is byte-identical to
SCOPer's. So the SCOPer branch should just BE briney, with the clonal step
enabled and defaulting to the hierarchical (SCOPer) method.

This also fixes the AWS/Fargate + primer concerns for free: briney already has
the full `aws_batch` Fargate profile, staged-channel refs (no host bind-mount),
dual C+V primer masking, and ParseDb.

## Approach
Restore briney's versions of all pipeline-logic files onto the SCOPer working
tree, then apply a minimal SCOPer-side delta.

### Step 1 — restore briney pipeline files (git checkout briney -- <paths>)
conf/{ec2,full,modules,test}.config; nextflow.config;
modules/local/changeo/{assigngenes,makedb,parsedb}/main.nf;
modules/local/presto/{gunzip,maskprimers_align}/main.nf;
modules/local/scoper/hierarchicalclones/main.nf;
subworkflows/local/{clonal_analysis,presto,vdj_assignment}.nf;
workflows/immcantation.nf;
assets/{cprimers,vprimers}_briney2019.fasta
(Leaves untracked SCOPer extras intact: bin/find_threshold.R,
modules/local/shazam/disttonearest/. Leaves README/.gitignore/subset_test as-is.)

### Step 2 — enable clonal + default to SCOPer (deltas)
- `workflows/immcantation.nf`: uncomment `CLONAL_ANALYSIS(VDJ_ASSIGNMENT.out.airr_tab)`.
- `nextflow.config`: `cloning_method='hierarchical'` (briney stays 'exact');
  `scoper_linkage='single'` (diagnosis sensitivity, Gupta 2017);
  add `clone_threshold_mode='fixed'` + `threshold_drift_tol=0.04`.

### Step 3 — re-wire distToNearest drift-check as OPT-IN (default off)
- `subworkflows/local/clonal_analysis.nf`: in the hierarchical branch, if
  `clone_threshold_mode=='disttonearest'` also run `SHAZAM_DISTTONEAREST(ch_grouped)`.
- `conf/modules.config`: add `withName:'SHAZAM_DISTTONEAREST'` block.
- Confirm `bin/find_threshold.R` is +x; shazam module already present.

### Step 4 — sanity check
`nextflow run ... -preview` for default (fixed) and disttonearest modes; confirm
DefineClones path untouched, SCOPer path runs, distToNearest only in opt-in mode.

## Net difference vs briney
1 param default (`cloning_method='hierarchical'`), 1 enabled call, single linkage,
+ optional distToNearest mode. Branch parity otherwise.

## Non-goals
- No change to briney's AWS profile, primers, refs, or ParseDb.
- distToNearest never feeds SCOPer; drift-check only.
