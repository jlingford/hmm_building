# FeFe HydDB: Pre-`hmmbuild` Pipeline

```text
all hits
-> hydrogenase domain subsequences
-> length filtering
-> CD-HIT representative sequences
-> MSA
-> FastTree tree QC
-> candidate sequence pruning
-> re-MSA
-> MSA column QC / gap-rich column trimming
-> final trimmed MSA for hmmbuild
```

```text
yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap80.trimmed.aln.faa
```

```text
gene family: FeFe hydrogenase
domain signature: Gene3D G3DSA:3.40.50.1780
minimum domain length: 100 aa
CD-HIT identity tested: 80%, 90%, 95%
CD-HIT identity selected: 80%
CD-HIT coverage: 80%
tree outlier rule selected: Q3 + 3 * IQR
MSA gap trimming tested: gap80, gap90, gap95
MSA gap trimming selected: gap80
tree-pruned sequences removed: 584
final sequence count: 26,901
final gap80-trimmed alignment length: 370 columns
```

## Cutoff Summary

This table summarizes all major cutoff choices before `hmmbuild`.

```text
step                         cutoff(s) considered/tested                       selected cutoff
domain signature              Gene3D G3DSA:3.40.50.1780                        Gene3D G3DSA:3.40.50.1780
FeFe minimum domain length     single selected threshold: 100 aa                 >= 100 aa
CD-HIT identity                80%, 90%, 95%                                     80%
CD-HIT coverage                single selected threshold: 80% shorter/longer      >= 80%
tree terminal_branch outlier   Q3+1.5IQR, Q3+3IQR                                Q3+3IQR
tree root_to_tip outlier       Q3+1.5IQR, Q3+3IQR                                Q3+3IQR
tree parent support warning    single selected threshold: 0.5                    < 0.5
tree near-min-length warning   single selected threshold: 120 aa                 <= 120 aa
MSA gap trimming               gap80, gap90, gap95                               gap80
```

Important distinction:

```text
"considered/tested" means alternative output files or counts were generated.
"single selected threshold" means one biologically or practically motivated threshold was used, but no alternative output series was generated in this run.
```

## Why We Do This Before `hmmbuild`

`hmmbuild` builds a profile HMM from a multiple sequence alignment.

The HMM quality depends strongly on the input alignment:

- if false positives are included, the HMM becomes too broad
- if fragments are included, the HMM may learn fragmentary patterns
- if long-branch outliers are included, the alignment and HMM can become noisy
- if the MSA contains many sparse insertion columns, the HMM becomes large and less focused
- if highly redundant sequences dominate, the HMM can overrepresent common clades

So before `hmmbuild`, we want:

```text
1. sequences that really represent the target domain
2. sequences long enough to contain the domain core
3. reduced redundancy
4. an MSA of representatives
5. tree-based removal of extreme sequence outliers
6. a re-alignment after pruning
7. removal of gap-rich columns
```

The result is a more compact and biologically focused seed alignment.

## Directory Overview

Main input data:

```text
data/00_all_hits/
data/01_subseqs/
data/02_subseqs_size_filtered/
```

Yuyao pipeline outputs:

```text
yuyao/CD-HIT/
yuyao/MSA/
yuyao/Tree/
```

Important scripts:

```text
yuyao/CD-HIT/scripts/run_cdhit_one.sh
yuyao/MSA/scripts/run_famsa_one.sh
yuyao/MSA/scripts/msa_column_qc_trim.py
yuyao/Tree/scripts/make_fasttree_safe_fasta.sh
yuyao/Tree/scripts/run_fasttree_one.sh
yuyao/Tree/scripts/tree_qc_visualize.py
yuyao/Tree/scripts/apply_tree_prune.py
```

Related detailed READMEs:

```text
yuyao/Tree/README_tree_qc_pruning.md
yuyao/MSA/README_msa_column_qc_trimming.md
```

## Step 0. Starting Point: Sequence Search Hits

Input:

```text
data/00_all_hits/fefe.allhits.faa
```

Meaning:

```text
Candidate FeFe sequences obtained from the initial search.
These are full-length protein hits, not yet trimmed to the hydrogenase domain.
```

Stats:

```text
file              num_seqs  sum_len     min_len  avg_len  max_len
fefe.allhits.faa  147,687   79,631,734  125      539.2    10,808
```

Why this step exists:

```text
The initial search aims to collect candidate homologs broadly.
At this point the sequences may include full-length proteins with extra domains.
We do not want to build the final HMM directly from full-length hits.
```

Why not use full-length sequences directly:

- FeFe proteins may have extra N-terminal or C-terminal regions
- different proteins may contain additional domains
- full-length variation can make MSA noisy
- HydDB profile should focus on the conserved hydrogenase domain

## Step 1. Extract Hydrogenase Domain Subsequences

Input:

```text
data/00_all_hits/fefe.allhits.faa
```

Output:

```text
data/01_subseqs/fefe.subseqs.Gene3D_3-40-50-1780.faa
```

Domain signature used:

```text
FeFe: Gene3D G3DSA:3.40.50.1780
```

Domain signature cutoff/choice note:

```text
For this FeFe run, the selected domain annotation was Gene3D G3DSA:3.40.50.1780.
Alternative domain signatures were not generated as separate FeFe output series in this yuyao run.
```

Why this Gene3D signature was selected:

```text
1. It marks the FeFe hydrogenase core structural domain.
2. It gives consistent domain boundaries across diverse sequences.
3. It maps to structural-domain style annotation, which is useful for downstream comparison to CATH/TED-style domain resources.
4. A consistent domain signature is preferable to mixing Pfam/CDD/Gene3D boundaries in one HMM seed alignment.
```

Meaning:

```text
Only the hydrogenase core domain subsequence was extracted from each full-length hit.
```

Stats:

```text
file                                  num_seqs  sum_len     min_len  avg_len  max_len
fefe.subseqs.Gene3D_3-40-50-1780.faa  106,776   28,797,823  9        269.7    466
```

Why this step is needed:

```text
hmmbuild should learn the target domain, not unrelated flanking domains.
Extracting domain subsequences focuses the downstream MSA and HMM on the FeFe hydrogenase core.
```

Effect:

```text
147,687 full-length candidate hits
-> 106,776 FeFe domain subsequences
```

Sequences not retained here may lack the selected Gene3D annotation or may not have a confident domain boundary.

## Step 2. Length Filter Domain Subsequences

Input:

```text
data/01_subseqs/fefe.subseqs.Gene3D_3-40-50-1780.faa
```

Output:

```text
data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa
```

Threshold:

```text
minimum FeFe domain subsequence length = 100 aa
```

Command concept:

```bash
seqkit seq -m 100 \
  data/01_subseqs/fefe.subseqs.Gene3D_3-40-50-1780.faa \
  > data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa
```

Stats:

```text
file                                               num_seqs  sum_len     min_len  avg_len  max_len
fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa  104,406   28,721,111  100      275.1    466
```

Why this step is needed:

```text
Very short sequences are likely fragments.
Fragments can cause bad alignments, long branches, and noisy HMM profiles.
```

Effect:

```text
106,776 domain subsequences
-> 104,406 length-filtered subsequences
```

Threshold rationale:

```text
100 aa is deliberately permissive.
It removes very small fragments while retaining broad FeFe sequence diversity.
```

Cutoff testing note:

```text
For this FeFe run, multiple minimum-length thresholds were not generated as separate output series.
The selected threshold was 100 aa.
```

Why not use a stricter minimum length at this stage:

```text
1. The goal was to remove obvious tiny fragments, not aggressively filter diversity.
2. Some legitimate FeFe domain subsequences may be partial but still informative.
3. Later tree QC and MSA QC steps provide additional ways to detect problematic short/fragmentary sequences.
4. A stricter cutoff such as 120 aa or 150 aa could remove real divergent or partial representatives before tree-based review.
```

How short sequences were handled later:

```text
After MSA/tree QC, sequences with ungapped_len <= 120 were flagged as near_min_length.
This warning was used for manual prioritization, but it was not an automatic deletion rule by itself.
```

## Step 3. Reduce Redundancy With CD-HIT

Input:

```text
data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa
```

Output used downstream:

```text
yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa
```

Script:

```text
yuyao/CD-HIT/scripts/run_cdhit_one.sh
```

Slurm scripts tested:

```text
yuyao/CD-HIT/slurm/run_cdhit_fefe_c80.slurm
yuyao/CD-HIT/slurm/run_cdhit_fefe_c90.slurm
yuyao/CD-HIT/slurm/run_cdhit_fefe_c95.slurm
```

Main command selected for downstream work:

```bash
CDHIT_IDENTITY=0.80 CDHIT_WORD_LENGTH=5 yuyao/CD-HIT/scripts/run_cdhit_one.sh \
    data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa \
    yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80 \
    "${SLURM_CPUS_PER_TASK}" \
    64000
```

CD-HIT parameters for the selected c80 run:

```text
-c 0.80     sequence identity threshold
-n 5        word length
-G 1        global identity
-aS 0.8     shorter sequence alignment coverage >= 80%
-aL 0.8     longer sequence alignment coverage >= 80%
-g 1        accurate clustering mode
-T 16       threads
-M 64000    memory MB
-d 0        keep full sequence descriptions
```

Parameters shared across all tested CD-HIT runs:

```text
-G 1        global identity
-aS 0.8     shorter sequence alignment coverage >= 80%
-aL 0.8     longer sequence alignment coverage >= 80%
-g 1        accurate clustering mode
-T 16       threads
-M 64000    memory MB
-d 0        keep full sequence descriptions
```

The intended comparison variable was sequence identity:

```text
c80 -> -c 0.80
c90 -> -c 0.90
c95 -> -c 0.95
```

The c90 run used:

```bash
CDHIT_IDENTITY=0.90 yuyao/CD-HIT/scripts/run_cdhit_one.sh \
    data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa \
    yuyao/CD-HIT/results/03_cdhit/fefe_min100_c90_cov80 \
    "${SLURM_CPUS_PER_TASK}" \
    64000
```

The c95 run used the default identity in `run_cdhit_one.sh`:

```text
CDHIT_IDENTITY default = 0.95
```

Stats for c80:

```text
input:  104,406 sequences
output:  27,485 representative sequences
```

Detailed stats:

```text
file                                                                            num_seqs  sum_len     min_len  avg_len  max_len
data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa 104,406   28,721,111  100      275.1    466
yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa                          27,485    6,875,119   100      250.1    466
```

Why this step is needed:

```text
The raw filtered sequence set is highly redundant.
Without dereplication, abundant groups can dominate the MSA and HMM.
CD-HIT reduces redundancy while retaining representatives across diversity.
```

Why 80% identity:

```text
80% identity is a relatively strong dereplication level.
It reduces the dataset enough for practical MSA/tree building while retaining broad clade diversity.
```

### CD-HIT Cutoffs Tested

Three identity thresholds were tested:

```text
c80: identity >= 80%, coverage >= 80%
c90: identity >= 90%, coverage >= 80%
c95: identity >= 95%, coverage >= 80%
```

All three started from the same length-filtered FeFe domain input:

```text
data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa
104,406 sequences
```

Results:

```text
cutoff   representatives   sum_len      min_len   avg_len   max_len
c80      27,485            6,875,119    100       250.1     466
c90      51,686            13,656,770   100       264.2     466
c95      70,444            19,042,569   100       270.3     466
```

Interpretation:

```text
c95 keeps the most sequences.
    This preserves more near-duplicate variation, but creates a much larger MSA/tree problem.

c90 is intermediate.
    It reduces redundancy but still leaves more than 50k representatives.

c80 is the strongest dereplication.
    It gives 27,485 representatives, which is much more practical for MSA, tree building, visualization, and manual QC.
```

Why c80 was selected:

```text
1. It reduced 104,406 filtered FeFe domain sequences to 27,485 representatives.
2. This sequence count was practical for FAMSA and FastTree.
3. It reduced overrepresentation of heavily sampled clades.
4. It still retained broad FeFe diversity because 80% identity is not too low for domain-level representatives.
5. The goal here is to build a robust seed HMM, not to include every near-duplicate sequence.
```

Biological meaning of choosing c80:

```text
c80 sacrifices some within-clade redundancy to emphasize cross-clade diversity.
This is appropriate for a profile HMM seed alignment, where redundant close homologs can overweight common groups.
```

Computational meaning of choosing c80:

```text
c90 and c95 would produce larger MSAs and larger trees.
That would make tree QC and visualization harder and slower.
```

Final decision:

```text
Use c80/cov80 representatives for the current FeFe HMM-building path.
```

## Step 4. Initial MSA With FAMSA

Input:

```text
yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa
```

Output:

```text
yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.famsa.aln.faa
```

Script:

```text
yuyao/MSA/scripts/run_famsa_one.sh
```

Slurm script:

```text
yuyao/MSA/slurm/run_famsa_fefe_c80.slurm
```

Command:

```bash
yuyao/MSA/scripts/run_famsa_one.sh \
    yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa \
    yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.famsa.aln.faa \
    "${SLURM_CPUS_PER_TASK}"
```

Tool:

```text
FAMSA
```

Stats:

```text
input FASTA:
  sequences: 27,485
  min length: 100
  avg length: 250.1
  max length: 466

initial MSA:
  sequences: 27,485
  alignment length: 3,799 columns
```

Why this step is needed:

```text
hmmbuild requires a multiple sequence alignment, not unaligned FASTA.
The initial MSA is also needed for tree building and sequence outlier detection.
```

Why FAMSA:

```text
FAMSA is fast enough for tens of thousands of protein sequences and suitable for large representative datasets.
```

Important:

```text
This initial MSA was not used directly for final hmmbuild.
It was used for tree construction and tree-based sequence QC.
```

## Step 5. Convert MSA IDs To Newick-Safe IDs

Input:

```text
yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.famsa.aln.faa
```

Outputs:

```text
yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa
yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.id_map.tsv
```

Script:

```text
yuyao/Tree/scripts/make_fasttree_safe_fasta.sh
```

Command:

```bash
yuyao/Tree/scripts/make_fasttree_safe_fasta.sh \
    yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.famsa.aln.faa \
    yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa \
    yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.id_map.tsv
```

What it does:

```text
original ID: GCF_017310645___2335:205-511
safe ID:    GCF_017310645___2335_205-511__seq14386
```

The script replaces non-safe characters with `_` and appends a unique `__seqN` suffix.

Why this step is needed:

```text
Newick parsers can have trouble with characters such as ':' in sequence IDs because ':' is used for branch lengths.
Safe IDs avoid ambiguity in FastTree output and downstream parsing.
```

The ID map is essential because:

```text
tree/alignment pruning uses safe_id
original FASTA pruning uses original_id
```

## Step 6. Build Initial FastTree

Input:

```text
yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa
```

Output:

```text
yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree.nwk
```

Script:

```text
yuyao/Tree/scripts/run_fasttree_one.sh
```

Slurm script:

```text
yuyao/Tree/slurm/run_fasttree_fefe_c80.slurm
```

Command:

```bash
yuyao/Tree/scripts/run_fasttree_one.sh \
    yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa \
    yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree.nwk
```

Tree stats:

```text
leaves: 27,485
nodes:  54,968
```

Why this step is needed:

```text
The tree allows us to detect sequence-level outliers that are difficult to catch from length alone.
Long terminal branches and unusually long root-to-tip distances can indicate fragments, bad alignments, false homologs, or extremely divergent sequences.
```

Important:

```text
This tree is a QC tool.
It is not the final product of HydDB.
```

## Step 7. Tree QC: Calculate Leaf Metrics

Input files:

```text
tree:      yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree.nwk
ID map:    yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.id_map.tsv
alignment: yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa
```

Script:

```text
yuyao/Tree/scripts/tree_qc_visualize.py
```

Command:

```bash
python yuyao/Tree/scripts/tree_qc_visualize.py \
  --tree yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree.nwk \
  --id-map yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.id_map.tsv \
  --alignment yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa \
  --outdir yuyao/Tree/results/06_tree_qc \
  --prefix fefe_min100_c80_cov80
```

Outputs:

```text
yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.leaf_qc.tsv
yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.candidate_prune.tsv
yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.summary.json
yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.tree_viewer.html
```

Metrics calculated for every leaf:

```text
safe_id
original_id
terminal_branch
root_to_tip
parent_support
ungapped_len
flags
```

Meaning:

```text
terminal_branch:
    the branch length directly leading to a leaf

root_to_tip:
    total branch length from tree root to that leaf

parent_support:
    support value for the immediate parent node

ungapped_len:
    number of non-gap residues for that sequence in the MSA
```

MSA use in this step:

```text
Only used to calculate per-sequence ungapped_len.
No MSA columns were deleted at this stage.
```

Tree use in this step:

```text
Used to calculate terminal_branch, root_to_tip, and parent_support.
```

## Step 8. Tree QC Thresholds

The script uses a robust outlier cutoff:

```text
outlier cutoff = Q3 + 3 * IQR
IQR = Q3 - Q1
```

This was applied separately to:

```text
terminal_branch
root_to_tip
```

### Tree Cutoff Options Compared

Two common IQR-style outlier cutoffs are:

```text
Q3 + 1.5 * IQR
Q3 + 3.0 * IQR
```

The `1.5 * IQR` rule is more sensitive. It catches more moderate outliers.

The `3.0 * IQR` rule is stricter. It catches only stronger outliers.

For this tree, the comparison was:

```text
metric             rule          cutoff              flagged rows
terminal_branch    Q3 + 1.5IQR   0.3583242           1,690
terminal_branch    Q3 + 3.0IQR   0.509346633           584
root_to_tip        Q3 + 1.5IQR   5.325347846999999      33
root_to_tip        Q3 + 3.0IQR   7.7513541209999985      1
```

Union of terminal-branch and root-to-tip outliers:

```text
Q3 + 1.5IQR: 1,708 candidate sequences
Q3 + 3.0IQR:   584 candidate sequences
```

Why `Q3 + 3 * IQR` was selected:

```text
1. The purpose was conservative sequence pruning.
2. We wanted to remove the strongest tree outliers without deleting too much real FeFe diversity.
3. Q3 + 1.5IQR would have flagged 1,708 sequences, which is much more aggressive.
4. Q3 + 3IQR flagged 584 sequences, a smaller and more defensible first-pass removal set.
5. Low support and near-minimum length were retained as warning labels instead of automatic deletion rules.
```

Final decision:

```text
Use Q3 + 3 * IQR for terminal_branch and root_to_tip outlier detection.
```

### Terminal Branch Threshold

Observed:

```text
Q1 = 0.106620145
Q3 = 0.207301767
IQR = 0.100681622
```

Cutoff:

```text
terminal_branch_outlier_cutoff = 0.509346633
```

Rule:

```text
terminal_branch >= 0.509346633
-> terminal_branch_outlier
```

### Root-To-Tip Threshold

Observed:

```text
Q1 = 1.2820040570000002
Q3 = 2.8993415729999996
IQR = 1.6173375159999994
```

Cutoff:

```text
root_to_tip_outlier_cutoff = 7.7513541209999985
```

Rule:

```text
root_to_tip >= 7.7513541209999985
-> root_to_tip_outlier
```

### Extra Warning Flags

These flags were calculated but were not sufficient by themselves to enter the prune list:

```text
parent_support < 0.5
-> low_parent_support

ungapped_len <= 120
-> near_min_length
```

Why they are warnings only:

```text
A sequence can have low local support or be short without necessarily being wrong.
These flags help prioritize manual review, but the automatic candidate pruning list is based on branch-length outliers.
```

## Step 9. Select Candidate Sequences For Pruning

Output:

```text
yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.candidate_prune.tsv
```

Selection rule:

```text
candidate_prune.tsv includes rows where:

terminal_branch_outlier
OR
root_to_tip_outlier
```

Equivalent numeric rule:

```text
terminal_branch >= 0.509346633
OR
root_to_tip >= 7.7513541209999985
```

Result:

```text
candidate prune sequences: 584
```

Meaning:

```text
These 584 sequences were flagged as strong tree-based outliers.
They are candidates for removal before building the final HMM seed alignment.
```

Important:

```text
This step removes full sequences later.
It does not remove MSA columns.
```

## Step 10. Apply Tree-Based Sequence Pruning

Script:

```text
yuyao/Tree/scripts/apply_tree_prune.py
```

This script removes:

```text
tree leaves
FASTA records
MSA sequence rows
```

It does not remove:

```text
MSA columns
```

### Prune The Tree And Safe-ID Alignment

Command:

```bash
python yuyao/Tree/scripts/apply_tree_prune.py \
  --tree yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree.nwk \
  --fasta yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree_safe.aln.faa \
  --prune-list yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.candidate_prune.tsv \
  --out-tree yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.nwk \
  --out-fasta yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.aln.faa
```

Output:

```text
yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.nwk
yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.aln.faa
```

Result:

```text
kept:    26,901 sequences
removed:    584 sequences
```

### Prune The Original Unaligned FASTA

Input:

```text
yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa
```

Command:

```bash
python yuyao/Tree/scripts/apply_tree_prune.py \
  --fasta yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa \
  --prune-list yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.candidate_prune.tsv \
  --prune-column original_id \
  --out-fasta yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.faa
```

Output:

```text
yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.faa
```

Why `--prune-column original_id`:

```text
The original CD-HIT FASTA uses original IDs, not safe IDs.
```

Result:

```text
kept:    26,901 sequences
removed:    584 sequences
```

## Step 11. Re-Align After Sequence Pruning

Input:

```text
yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.faa
```

Output:

```text
yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.pruned.famsa.aln.faa
```

Slurm script:

```text
yuyao/MSA/slurm/run_famsa_fefe_c80_pruned.slurm
```

Command in Slurm:

```bash
yuyao/MSA/scripts/run_famsa_one.sh \
    yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.faa \
    yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.pruned.famsa.aln.faa \
    "${SLURM_CPUS_PER_TASK}"
```

Stats:

```text
input pruned FASTA:
  sequences: 26,901
  sum length: 6,757,423
  min length: 100
  avg length: 251.2
  max length: 466

new pruned MSA:
  sequences: 26,901
  alignment length: 3,764 columns
```

Why re-align after pruning:

```text
The original MSA was built with the outlier sequences included.
After removing 584 outliers, the alignment should be rebuilt so the remaining sequences define the alignment without those outliers.
```

Important:

```text
The re-MSA changed alignment length from 3,799 columns to 3,764 columns.
This is a natural result of re-aligning, not manual column trimming.
```

## Step 12. MSA Column QC

Input:

```text
yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.pruned.famsa.aln.faa
```

Script:

```text
yuyao/MSA/scripts/msa_column_qc_trim.py
```

Main output QC table:

```text
yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.column_qc.tsv
```

For each alignment column, the script calculates:

```text
column_1based
gap_count
gap_fraction
non_gap_count
non_gap_fraction
dominant_residue
dominant_residue_fraction_of_non_gap
```

Formula:

```text
gap_fraction = gap_count / number_of_sequences
```

Gap characters counted:

```text
-
.
```

Why this step is needed:

```text
Large MSAs often contain sparse insertion columns.
These are columns where only a few sequences have residues and most sequences have gaps.
Such columns can make the HMM unnecessarily large and noisy.
```

Observed gap distribution:

```text
input alignment length: 3,764 columns
minimum gap fraction:  0.007583361213337794
Q1 gap fraction:       0.9995910932679083
median gap fraction:   0.9999256533214379
Q3 gap fraction:       0.9999628266607189
maximum gap fraction:  0.9999628266607189
```

This means most alignment columns are nearly all gaps.

## Step 13. Trim Gap-Rich MSA Columns

Three thresholds were generated for comparison:

```text
gap80: remove columns with gap_fraction > 0.80
gap90: remove columns with gap_fraction > 0.90
gap95: remove columns with gap_fraction > 0.95
```

All three used the same input alignment:

```text
yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.pruned.famsa.aln.faa
26,901 sequences x 3,764 columns
```

The comparison variable was the maximum allowed gap fraction per column.

Meaning:

```text
gap80:
    remove columns where more than 80% of sequences are gaps

gap90:
    remove columns where more than 90% of sequences are gaps

gap95:
    remove columns where more than 95% of sequences are gaps
```

These thresholds test how aggressive column trimming should be:

```text
lower threshold -> more aggressive trimming
higher threshold -> less aggressive trimming
```

### Recommended Main Threshold: gap80

Command:

```bash
python yuyao/MSA/scripts/msa_column_qc_trim.py \
  --alignment yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.pruned.famsa.aln.faa \
  --outdir yuyao/MSA/results/05_column_qc \
  --prefix fefe_min100_c80_cov80.pruned.famsa \
  --gap-threshold 0.8
```

Rule:

```text
remove if gap_fraction > 0.8
keep if gap_fraction <= 0.8
```

Output:

```text
yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap80.trimmed.aln.faa
```

Stats:

```text
input columns:   3,764
kept columns:      370
removed columns: 3,394
removed fraction: 0.9017003188097769
```

Why gap80 is recommended:

```text
It removes sparse insertion columns and keeps the shared core alignment.
The FeFe domain core should be represented by columns present in a substantial fraction of sequences.
```

### Comparison: gap90

Output:

```text
yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap90.trimmed.aln.faa
```

Stats:

```text
kept columns:      372
removed columns: 3,392
```

### Comparison: gap95

Output:

```text
yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap95.trimmed.aln.faa
```

Stats:

```text
kept columns:      378
removed columns: 3,386
```

Interpretation:

```text
gap80, gap90, and gap95 are very similar here.
This shows that most removed columns are not borderline; they are almost completely gaps.
```

### Why gap80 Was Selected

Comparison:

```text
threshold   kept columns   removed columns   final alignment length
gap80       370            3,394             370
gap90       372            3,392             372
gap95       378            3,386             378
```

The differences are small:

```text
gap90 keeps only 2 more columns than gap80.
gap95 keeps only 8 more columns than gap80.
```

This tells us:

```text
Most removed columns have extremely high gap fractions.
The choice between 80%, 90%, and 95% does not greatly change the core alignment.
```

Why gap80 is the recommended first `hmmbuild` input:

```text
1. It is the most conservative about excluding sparse insertion columns.
2. It keeps columns represented in at least 20% of sequences.
3. It produces the most compact alignment: 370 columns.
4. Since gap90 and gap95 add only a few columns, gap80 is unlikely to remove a large amount of shared core signal.
5. A compact core-domain HMM is a good first profile for downstream validation.
```

Why keep gap90 and gap95:

```text
They are useful comparison alignments.
If the gap80 HMM is too strict or loses sensitivity, gap90 or gap95 can be used to build alternative HMMs.
```

Final decision:

```text
Use gap80 trimmed alignment for the first FeFe hmmbuild.
Keep gap90 and gap95 as alternative inputs for comparison if needed.
```

## Final Pre-`hmmbuild` Output

Recommended alignment for first HMM build:

```text
yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap80.trimmed.aln.faa
```

Final alignment stats:

```text
sequences: 26,901
columns:     370
```

This alignment has gone through:

```text
domain extraction
length filtering
CD-HIT redundancy reduction
initial MSA
tree QC
sequence pruning
re-MSA
gap-rich column trimming
```

## What Has Not Yet Been Done

At this point, we have not yet:

```text
run hmmbuild
calibrated the HMM
run hmmsearch validation
selected score cutoffs
selected coverage cutoffs
benchmarked false positives
compared against HydDB_v1 or curated references
```

These are post-alignment HMM validation steps.

## Next Command: `hmmbuild`

Suggested next command:

```bash
mkdir -p yuyao/HMM/results

hmmbuild \
  yuyao/HMM/results/fefe_min100_c80_cov80.pruned.gap80.hmm \
  yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap80.trimmed.aln.faa
```

Expected input:

```text
trimmed aligned FASTA
```

Expected output:

```text
HMMER profile HMM
```

After this, the HMM should be tested with:

```text
hmmsearch
hmmscan
score distribution analysis
domain coverage analysis
false-positive review
```

## Full FeFe Count Summary

```text
full-length all hits:
  147,687 sequences

domain subsequences:
  106,776 sequences

length-filtered FeFe domain subsequences:
  104,406 sequences

CD-HIT c80 representatives:
  27,485 sequences

tree-prune candidates removed:
  584 sequences

pruned FASTA:
  26,901 sequences

pruned re-MSA:
  26,901 sequences x 3,764 columns

gap80 final trimmed MSA:
  26,901 sequences x 370 columns
```

## Reproducibility Checklist

Before running `hmmbuild`, confirm these files exist:

```bash
ls -lh \
  data/02_subseqs_size_filtered/fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa \
  yuyao/CD-HIT/results/03_cdhit/fefe_min100_c80_cov80.faa \
  yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.famsa.aln.faa \
  yuyao/Tree/results/05_trees/fefe_min100_c80_cov80.fasttree.nwk \
  yuyao/Tree/results/06_tree_qc/fefe_min100_c80_cov80.candidate_prune.tsv \
  yuyao/Tree/results/07_pruned/fefe_min100_c80_cov80.candidate_pruned.faa \
  yuyao/MSA/results/04_alignments/fefe_min100_c80_cov80.pruned.famsa.aln.faa \
  yuyao/MSA/results/05_column_qc/fefe_min100_c80_cov80.pruned.famsa.gap80.trimmed.aln.faa
```

If all files exist, the pipeline is ready for `hmmbuild`.
