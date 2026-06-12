# Building new HMM profiles for HydDB_v2

## 1. Sequence search

Sequences obtained by mmseqs search of GlobDB:

Notes:

- query fastas = HydDB_v1 sequences
- target fastas = each genome .faa file from the GlobDB individually
- used `--num-iteration 2` for PSI-BLAST style search
- `-c` 0.8

Stats:
```
file              format  type     num_seqs      sum_len  min_len  avg_len  max_len
fefe.allhits.faa  FASTA   Protein   147,687   79,631,734      125    539.2   10,808
nife.allhits.faa  FASTA   Protein   353,837  155,084,992       63    438.3    1,354
```

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

From experience, these have good coverage of each hydrogenase core structural domain, and have an equivalent code in the CATH structural database, which can be searched with Foldseek.
These domain codes are also found in the TED structural database.

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

## 4. Clustering

Can either use mmseqs cluster or cd-hit.
