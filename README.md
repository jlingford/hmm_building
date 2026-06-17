# Building new HMM profiles for HydDB_v2

Overarching goal: 

1. build new HMM profiles for the FeFe and NiFe to serve as a broad hmmsearch tool
    - use this to search for new remote hyd homologs that may be missed by the HydDB_v1 profiles

2. build new HMM profiles for each orthogroup within each hyd family, and use these as hmmscan group classification tools
    - e.g., build 4 profiles for the 4 NiFe groups --> run 1 query against all 4 --> whichever is the top scoring profile is the classification assignment to the query

## 0. Background -- Sequence search

Sequences obtained by mmseqs search of GlobDB:

Notes:

- query fastas = HydDB_v1 sequences
- target fastas = each genome .faa file from the GlobDB individually
- used `--num-iteration 2` for PSI-BLAST style search
- `-c 0.8`
- `-s 7.5`

Stats:
```
file              format  type     num_seqs      sum_len  min_len  avg_len  max_len
fefe.allhits.faa  FASTA   Protein   147,687   79,631,734      125    539.2   10,808
nife.allhits.faa  FASTA   Protein   353,837  155,084,992       63    438.3    1,354
```

## 1. Strategy for going from sequence to MSA to HMM profile:

#### Strategy 1: 

Start with subsequences of the hydrogenase domain, extracted from full length proteins.
Makes downstream alignment simple. 

1. Extract core hyd domain from full-length seqs by HMM scan --> get domain coordinates --> bedtools extract subseq
2. Basic filtering
    - remove very short sequences
    - remove non CxxC...CxxC seqs for NiFe
3. Reduce sequence redundancy to min-seq-id 80% and cov 80% with mmseqs cluster
4. Align with FAMSA
    - Check tree with FastTree

IF NEEDED:

5. Coarse grain filter to remove/prune sequences that significantly disrupt MSA
6. Align with FAMSA
    - Check tree with FastTree

7. HMMbuild
8. Test HMM
    - Check that HMM recovers the query set with a hmmsearch
    - search GlobDB, inspect low e-value hits for false positives
    - compare with old HMMs

## 2. Extraction of hydrogenase domain subseqs from full length seqs

### Annotation of sequences with InterProScan

Slurm script for running InterProScan on HPC cluster

```bash
#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J interproscan
#SBATCH --mem=60000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --account=rp24
#SBATCH --partition=genomicsb
#SBATCH --qos=genomicsbq
#SBATCH --time=0-24:00:00
#SBATCH --error=logs/%j.split.err
#SBATCH --output=logs/%j.split.out

# #SBATCH --mail-user=james.lingford@monash.edu
# #SBATCH --mail-type=BEGIN,END,FAIL,TIME_OUT

# # get data and install (do once)
# curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/alt/interproscan-data-5.73-104.0.tar.gz
# curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/alt/interproscan-data-5.73-104.0.tar.gz.md5
# md5sum -c interproscan-data-5.73-104.0.tar.gz.md5
# tar -pxzf interproscan-data-5.73-104.0.tar.gz
# mkir input temp output
# wget -O input/e-coli.fa 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29'
# singularity pull docker://interpro/interproscan:latest
#
# # for first time setup...
# inputpath="e-coli.fa"
# module load miniforge3
# cd ./interproscan-5.73-104.0
# python3 setup.py -f interproscan.properties
# cd ..

# input
# WARN: must supply the INPUTFASTA path with the "./" prefix in "./path/to/fasta.faa", otherwise won't work
INPUTFASTA=$1
JOBNAME=${2:-$(date +%F)}
# INPUTFASTA='./input/ivysaur/nife/nife_bighmmrivysaur.minus_bulbasaur_seqs.no_ambig_seqs.faa'
# INPUTFASTA='./input/ivysaur/fefe/fefe_ivysaur.minus_bulbasaur_seqs.no_ambig_seqs.faa'
# JOBNAME="$(date +%F)"

# WARN: must supply the INPUTFASTA path with the "./" prefix in "./path/to/fasta.faa", otherwise won't work
# replace "./" with "/input/" in filepath
inputpath=${INPUTFASTA/.\/input\//}

# set outdir name
stem=${inputpath##*/}
stem=${stem%.*}
outpath=${stem}
# WARN: double check this
dir="${JOBNAME}/${outpath}"
mkdir -p ./output/"${dir}"

# run container
singularity exec \
    -B $PWD/interproscan-5.73-104.0/data:/opt/interproscan/data \
    -B $PWD/input:/input \
    -B $PWD/temp:/temp \
    -B $PWD/output:/output \
    interproscan_5.73-104.0.sif \
    /opt/interproscan/interproscan.sh \
    --input /input/${inputpath} \
    --disable-precalc \
    --output-dir /output/${dir} \
    --tempdir /temp \
    --cpu 8
```

Split the giant fasta file into many smaller files containing 2000 seqs each, and send each of these as input.
This speeds up annotation speed, allows some level of parallelism.

E.g., batch InterProScan job launcher on split fasta files:

```bash
#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J james
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --account=rp24
#SBATCH --partition=genomics
#SBATCH --qos=genomics
#SBATCH --time=4:00:00
#SBATCH --mail-user=james.lingford@monash.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_OUT
#SBATCH --error=log-%j.err
#SBATCH --output=log-%j.out

# ---
# WARN: change these each time
# INPUTDIR='./input/NiFe_hyddb1_cell2023-globdb_covm2-iter2-mmseqs_search-UNIQUE_hits_SPLIT_FASTAS'
# INPUTDIR='./input/Guaymas_faa_combined/Guaymas_MAGs_combined_SPLIT_FASTAS'
# jobname='Guaymas_MAGs_faa'
INPUTDIR="./input/ivysaur/nife/nife_bighmmrivysaur.minus_bulbasaur_seqs.no_ambig_seqs_SPLIT_FASTAS"
jobname="nife_bighmmrivysaur_minus_bulbasaur_noabigseqs_splitinput"

count=-1
for file in ${INPUTDIR}/*.fa*; do
    ((count++))
    sbatch -J proscan_${count} slurm_docker_proscan.sh $file $jobname
done
```

Combine all individual .tsv output files into one .tsv file.

Domain boundaries can be extracted from this file.

### Choosing domain signature

The 4th and 5th columns in the interproscan.tsv output are the annotation analysis type (i.e., Pfam, CDD, Gene3D) and annotation code (e.g., PF02906).
The 7th and 8th columns are the START and STOP boundaries of the domain, respectively

Best to stick with one choice of annotation type + code per gene when extracting subseqs, for consistency.

I've chosen the following annotation codes to focus on:

- NiFe: Gene3D - G3DSA:1.10.645.10
- FeFe: Gene3D - G3DSA:3.40.50.1780

From experience, these have good coverage of each hydrogenase core structural domain, and have an equivalent code in the CATH50 structural database, which can be searched with Foldseek.

These domain codes are also found in the TED/CATH structural database.

### Create bed file from InterProScan tsv file

To extract domain subseqs, we can use `bedtools`.
First, we need to create a .bed file in the bed file format (3 columns: gene_name, start_position, stop_position)

```bash
# just an example
tsv=interproscan_output_file.tsv

# for FeFe
grep "G3DSA:3.40.50.1780" $tsv | cut -d$'\t' -f1,7,8 > ${tsv/.tsv/.bed}

# for NiFe
grep "G3DSA:3.40.50.1780" $tsv | cut -d$'\t' -f1,7,8 > ${tsv/.tsv/.bed}
```

### Extract domain subseqs

Now we can extract the subseqs with bedtools, like so:

```bash
bedtools getfasta \
    -fi input.faa \
    -bed bedfile.bed \
    -fo output.subseq.faa
```

Sequence stats after subseq extraction:

```
file                                  format  type     num_seqs     sum_len  min_len  avg_len  max_len
fefe.subseqs.Gene3D_3-40-50-1780.faa  FASTA   Protein   106,776  28,797,823        9    269.7      466
nife.subseqs.Gene3D_1-10-645-10.faa   FASTA   Protein   262,363  56,251,123       24    214.4      719
```

### Considerations

Not all the hits in the allhits.faa file will have a hydrogenase annotation from InterProScan, so a large chunk of sequences are during the subseq extraction.
This is likely due to the falling below the significance threshold for hydrogenase domain annotation.
So we are probably losing some of the sequence diversity here, but this might be a worthwhile tradeoff since the downstream technical challenges of MSA generation will probably be easier.

A different approach to try might be to use mmseqs clusters with `--min-seq-id 0` to gather highly divergent sequences, generate MSAs, and find domain boundaries that way.

## 3. Some basic filtering

Some of the subseqs extracted result in seqs that are far too small to be complete hydrogenase domains.
E.g., a sequence as small as 9 or 24 amino acids is just a fragment of a hydrogenase, and very small seqs will make MSA generation more complicated

Filtered seqs to have a minimum sequence length of:

- FeFe: 100 min length
- NiFe: 200 min length

These are arbitrarily chosen, and quite liberal, minumum lengths.

Filter length with:

```bash
#e.g.,
faa=nife.subseq.faa

# for NiFe
seqkit seq -m 200 $faa > ${faa/.faa/.minseqlen200.faa}

# for FeFe
seqkit seq -m 100 $faa > ${faa/.faa/.minseqlen100.faa}
```

Stats post-filtering for length

```
file                                               format  type     num_seqs     sum_len  min_len  avg_len  max_len
fefe.subseqs.Gene3D_3-40-50-1780.minseqlen100.faa  FASTA   Protein   104,406  28,721,111      100    275.1      466
nife.subseqs.Gene3D_1-10-645-10.minseqlen200.faa   FASTA   Protein   107,850  49,053,374      200    454.8      719
```

### Filtering NiFe by CxxC...CxxC motif

- even with filtering, there are sequences that do not contain the double CxxC motif in NiFe subseq set
- want putatively active NiFe hyds for this HMM, since we want to search for active hyds

A liberal filter could be: CxxC...(at least 50 amino acids)...CxxC

```bash
seqkit grep -s -r -p '"C..C.{50,}C..C"' nife.subseqs.Gene3D_1-10-645-10.minseqlen200.faa | grep -c ">"
99503
```

This gives us:

```text
file                                                          format  type     num_seqs     sum_len  min_len  avg_len  max_len
nife.subseqs.Gene3D_1-10-645-10.minseqlen200.CxxC_filter.faa  FASTA   Protein    99,503  45,642,942      291    458.7      719
nife.subseqs.Gene3D_1-10-645-10.minseqlen200.faa              FASTA   Protein   107,850  49,053,374      200    454.8      719
```

So a removal of 8,347 sequences

Looking into these sequences:

```bash
grep ">" nife.subseqs.Gene3D_1-10-645-10.minseqlen200.CxxC_filter.faa | sed 's/>//' > nifepostfilterlist.txt

seqkit grep -v -n -f nifepostfilterlist.txt nife.subseqs.Gene3D_1-10-645-10.minseqlen200.faa > nife.seq_removed_in_CxxC_filter.faa                 

seqkit stats *.faa

# file                                                          format  type     num_seqs     sum_len  min_len  avg_len  max_len
# nife.seq_removed_in_CxxC_filter.faa                           FASTA   Protein     8,347   3,410,432      200    408.6      690
# nife.subseqs.Gene3D_1-10-645-10.minseqlen200.CxxC_filter.faa  FASTA   Protein    99,503  45,642,942      291    458.7      719
# nife.subseqs.Gene3D_1-10-645-10.minseqlen200.faa              FASTA   Protein   107,850  49,053,374      200    454.8      719
```

## 4. Clustering to reduce sequence redundancy

Clustering this subset of CxxC filtered seqs with `mmseqs cluster` (Using sbatch script on M3).

Parameters:

- `--cov-mode 0 --cluster-mode 0`: this sets it to 'greedy' clustering
- `--cov 0.8`: must have >=80% seq coverage over b/w query and target
- `--min-seq-id 0.8`: reduce sequence redundancy

This produced:

- 20,385 cluster representative seqs

```
file                                                                                   format  type     num_seqs    sum_len  min_len  avg_len  max_len
nife.subseqs.Gene3D_1-10-645-10.minseqlen200.CxxC_filter.mmseqs_mode00_cov80_id80.faa  FASTA   Protein    20,385  9,030,229      291      443      715
```

## 5. FAMSA alignment

We compute MSAs using FAMSA, which is the fastest alignment tool available. Its accuracy also exceeds that of MAFFT, ClustalO, Kalign (MUSCLE5 beats it on inputs <12,000 seqs): <https://www.nature.com/articles/s41587-026-03095-3>.

Running famsa:

```bash
#!/usr/bin/bash

# make sure famsa binary is available on PATH
#
# famsa usage:
# famsa [options] infile outfile
#
# see `famsa -h`
#
# seems like each file needs to be generated with separate runs... won't output all in one go
#

# positional args
infaa=$1    # WARN: expects `.faa` extension
CPU=${2:-0} # default to 0 cores if not set

outaln=${infaa%.*}.famsa.aln
guidetree=${infaa%.*}.guidetree.dnd
log="$(dirname "$infaa")/famsa.log"

########### Just get the MSA ############

famsa -v -t $CPU $infaa $outaln 2>$log

####### Saving the guide tree #############

famsa -gt sl -gt_export -v -t $CPU $infaa $guidetree 2>$log
```

Completes alignment in ~14 min on my laptop with 11 cpu

### Tree building

With the MSA and guidetree, we can make a phylogenetic tree.
The phylogenetic tree can help identify any 'spurious' sequences.
Using 'VeryFastTree' for this (which is a re-write of FastTree):

```bash
VeryFastTree -threads 11 \
    -intree1 $guidetree \
    $inaln \
    > tree.nwk
```

Completes in ~14 min on my laptop.

## 6. MSA QC / Coarse filtering

Coarse filtering of 'spurious' sequences might be needed to help reduce number of gaps.

>NOTE:
This is different from 'trimming', which removes columns (e.g., trimAl or bmge).
We don't want to remove columns, as this changes the underlying sequence of all sequences in the MSA.
Trimming is needed for building accurate phylogenetic trees, but not necessarily fit for building HMM profiles for searching and classification of real sequences.
Coarse filtering involves removing whole sequences, or 'pruning'

Testing two different option for coarse filtering:

1. Based on IQR of tree branches ("IQR_prune")
2. Based on MSA occupancy ("Kay_prune")

- IQR_prune is an idea developed by Yuyao. Basic idea is:
    - Calculate tree metrics:
        - terminal branch length
        - root-to-tip distance
        - parent node support
        - ungapped sequence length in the MSA
    - Then use a conservative Q3 + 3 * IQR (interquartile range) threshold on terminal branch length and root-to-tip distance to find sequences for pruning

- Kay_prune is based on the description of a coarse filtering method in Kay et al. 2025 Nature <https://www.nature.com/articles/s41586-025-09808-z>
    - their script is unavailable, but according to their description:

    >The coarse filter is applied first, it is based on transforming the protein sequence alignment
    into a column occupancy matrix, which has the format
        - If there is a column which only has one row with an amino acid position, its value is one.
        - If there is a column where there is occupancy in every row, its value regardless of amino acid, is equal to the number of rows. The coarse filter then operates on this transformed matrix, to remove rows where
        - Have less than 80% the average row sum
        - Have less than 80% the average column score
        - Have more than 10% of columns having less than 10% occupancy
        - Has more than 10% of columns that break those of more than 90% occupancy
        - After poorly aligned sequences were removed, the sequences were realigned with MAFFT as described above.

Yuyao has made python scripts for both.

Running IQR_prune:

```bash
# 1. build tree viz and output stats
python tree_qc_visualize.py \
    --tree tree.nwk \
    --outdir iqr_filter_output

# 2. Extract sequence list from candidate_prune.tsv
awk 'NR>1{print $1}' $prunetsv > prunelist.txt

# 3. run partner script, which prunes the seqs from the tree and MSA
python apply_tree_prune.py \
    --tree $tree \
    --fasta $msa \
    --out-tree prune.nwk prune.aln


```

Stats output:

```json
{
  "n_leaves": 20385,
  "n_nodes": 40768,
  "terminal_branch_q1": 0.12614,
  "terminal_branch_q3": 0.20514,
  "terminal_branch_outlier_cutoff": 0.44214,
  "root_to_tip_q1": 2.8839699999999997,
  "root_to_tip_q3": 4.11093,
  "root_to_tip_outlier_cutoff": 7.79181,
  "candidate_prune_count": 371
}
```


Running Kay_prune:

```bash
python msa_coarse_occupancy_filter.py \
    --alignment $msa \
    --outdir msa_coarse_filter
```

Then rerun FAMSA and VeryFastTree on pruned MSAs

### Comparing MSAs

Comparison of:

1. MSA without pruning
2. MSA after IQR_prune and FAMSA rerun
3. MSA after Kay_prune and FAMSA rerun

Comparison (note: files have been renamed for brevity)

```
file                                             num_seqs      sum_len  min_len  avg_len  max_len
nife.cluster_reps.famsa.aln                        20,385  114,747,165    5,629    5,629    5,629
nife.cluster_reps.famsa.IQR_prune.famsa.aln        20,014  111,077,700    5,550    5,550    5,550
nife.cluster_reps.famsa.Kay_prune.famsa.aln        19,755   94,744,980    4,796    4,796    4,796
```

Observations:

- visually inspecting the MSAs using `termal`, I can see that the NiFe L1 and L2 motifs are well aligned.
- the middle portion of the alignment is mostly gaps, as expected since this is the least well conserved region between different NiFe groups
- either pruning method does little to remove this large middle gap region...
- but the aligned regions look sensible to me, and don't have that "tetris" look you see in bad MSAs
- so this alignment is probably as good as it can get when including all these remote orthologs in the same alignment?

Stats about these MSAs:

```bash
for aln in ./*.aln; do
    echo "esl-alistat for ${aln}"
    esl-alistat $aln
done
```

```
esl-alistat for ./nife.cluster_reps.famsa.aln
Alignment number:    1
Format:              aligned FASTA
Number of sequences: 20385
Alignment length:    5629
Total # residues:    9030229
Smallest:            291
Largest:             715
Average length:      443.0
Average identity:    26%
//
esl-alistat for ./nife.cluster_reps.famsa.IQR_prune.famsa.aln
Alignment number:    1
Format:              aligned FASTA
Number of sequences: 20014
Alignment length:    5550
Total # residues:    8874351
Smallest:            291
Largest:             715
Average length:      443.4
Average identity:    26%
//
esl-alistat for ./nife.cluster_reps.famsa.Kay_prune.famsa.aln
Alignment number:    1
Format:              aligned FASTA
Number of sequences: 19755
Alignment length:    4796
Total # residues:    8710751
Smallest:            324
Largest:             587
Average length:      440.9
Average identity:    26%
//
```

## 7. Build HMMs

Building HMM profiles of all three alignments; will choose best at a later stage after some comparisons

```bash
for aln in ./*.aln; do
    hmmbuild ${aln%.*} $aln
done
```

Stats about these hmm profiles with `hmmstat`

```
hmmstat for ./nife.cluster_reps.famsa.hmm
# idx  name                 accession        nseq eff_nseq      M relent   info p relE compKL
# ---- -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
1      nife.cluster_reps.famsa -               20385     6.94    464   0.59   0.57   0.43   0.01


hmmstat for ./nife.cluster_reps.famsa.IQR_prune.famsa.hmm
# idx  name                 accession        nseq eff_nseq      M relent   info p relE compKL
# ---- -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
1      nife.cluster_reps.famsa.IQR_prune.famsa -               20014     6.78    464   0.59   0.57   0.43   0.01


# idx  name                 accession        nseq eff_nseq      M relent   info p relE compKL
# ---- -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
1      nife.cluster_reps.famsa.Kay_prune.famsa -               19755     6.34    449   0.59   0.57   0.44   0.01
```


