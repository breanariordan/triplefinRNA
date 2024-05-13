# Getting basic assembly metrics
You should now have a folder labelled RF_trinity_output (or whatever you named your output file) as well as a general output file RF_trinity_output.Trinity.fasta
For these steps, we are going back to the smaller 2x 51 bp paired-end reads that we perform quality control and trimming on in [Quality_control_and_trimming](https://github.com/breanariordan/triplefinRNA/blob/main/Quality_control_and_trimming.md)
For this, we need to make another samples_file with all of the smaller samples (see [samples_file_small.txt](https://github.com/breanariordan/triplefinRNA/blob/main/samples_file_small.txt) ).

```
cd nobackup/denovosourcefiles/scheduler/
/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/util/TrinityStats.pl RF_trinity_output/Trinity.tmp.fasta > trinityStats.log
```

# Obtaining the iscountmatrix

```
nano iscountmatrix.sl    # create slurm script
```

The following should be input as a slurm script

```
#!/bin/bash -e

#SBATCH --job-name=trinity-iscountmatrix
#SBATCH --account=uoo03946
#SBATCH --time=7:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --hint=nomultithread

# user specific environment

module purge
module load Trinity/2.14.0-gimkl-2022a
module load Miniconda3
module load Perl/5.38.2-GCC-12.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load Bowtie/1.3.1-GCC-11.3.0
module load RSEM/1.3.3-gimkl-2022a

# Trinity command

srun perl /opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/util/align_and_estimate_abundance.pl \
  --transcripts RF_trinity_output.Trinity.fasta \
  --seqType fq \
  --samples_file samples_file_small.txt \
  --est_method RSEM --aln_method bowtie2 \
  --trinity_mode \
  --prep_reference \
  --thread_count 16 \
  --coordsort_bam > bowtie-rsem_align_and_estimate_abundance.log
```

```
sbatch iscountmatrix.sl    # submit slurm job
```

# Organise data
The output will give you a folder for each of your samples with RSEM.gene.results and RSEM.isoform.results. We want to use the **gene** results

```
mkdir RSEM_results
cd RSEM_results
ln -s  ../10L1heart/RSEM.genes.results 10L1heart.txt
ln -s  ../10L2heart/RSEM.genes.results 10L2heart.txt
ln -s  ../10L3heart/RSEM.genes.results 10L3heart.txt
ln -s  ../10L4heart/RSEM.genes.results 10L4heart.txt
ln -s  ../10L5heart/RSEM.genes.results 10L5heart.txt
ln -s  ../10L1brain/RSEM.genes.results 10L1brain.txt
ln -s  ../10L2brain/RSEM.genes.results 10L2brain.txt
ln -s  ../10L3brain/RSEM.genes.results 10L3brain.txt
ln -s  ../10L4brain/RSEM.genes.results 10L4brain.txt
ln -s  ../10L5brain/RSEM.genes.results 10L5brain.txt
ln -s  ../10N1heart/RSEM.genes.results 10N1heart.txt
ln -s  ../10N2heart/RSEM.genes.results 10N2heart.txt
ln -s  ../10N3heart/RSEM.genes.results 10N3heart.txt
ln -s  ../10N4heart/RSEM.genes.results 10N4heart.txt
ln -s  ../10N5heart/RSEM.genes.results 10N5heart.txt
ln -s  ../10N1brain/RSEM.genes.results 10N1brain.txt
ln -s  ../10N2brain/RSEM.genes.results 10N2brain.txt
ln -s  ../10N3brain/RSEM.genes.results 10N3brain.txt
ln -s  ../10N4brain/RSEM.genes.results 10N4brain.txt
ln -s  ../10N5brain/RSEM.genes.results 10N5brain.txt
ln -s  ../22L1heart/RSEM.genes.results 22L1heart.txt
ln -s  ../22L2heart/RSEM.genes.results 22L2heart.txt
ln -s  ../22L3heart/RSEM.genes.results 22L3heart.txt
ln -s  ../22L4heart/RSEM.genes.results 22L4heart.txt
ln -s  ../22L5heart/RSEM.genes.results 22L5heart.txt
ln -s  ../22L1brain/RSEM.genes.results 22L1brain.txt
ln -s  ../22L2brain/RSEM.genes.results 22L2brain.txt
ln -s  ../22L3brain/RSEM.genes.results 22L3brain.txt
ln -s  ../22L4brain/RSEM.genes.results 22L4brain.txt
ln -s  ../22L5brain/RSEM.genes.results 22L5brain.txt
ln -s  ../22N1brain/RSEM.genes.results 22N1brain.txt
ln -s  ../22N2brain/RSEM.genes.results 22N2brain.txt
ln -s  ../22N3brain/RSEM.genes.results 22N3brain.txt
ln -s  ../22N4brain/RSEM.genes.results 22N4brain.txt
ln -s  ../22N5brain/RSEM.genes.results 22N5brain.txt
ln -s  ../22N1heart/RSEM.genes.results 22N1heart.txt
ln -s  ../22N2heart/RSEM.genes.results 22N2heart.txt
ln -s  ../22N3heart/RSEM.genes.results 22N3heart.txt
ln -s  ../22N4heart/RSEM.genes.results 22N4heart.txt
ln -s  ../22N5heart/RSEM.genes.results 22N5heart.txt
```

The next step, merging the same into a single matrix, is done in R Studio

```R
# set working directory to RSEM_results
library("tximport")
files <- dir()
txi_rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi-rsem$counts)
colnames(txi.rsem$counts) <- gsub(".txt","",files)
head(txi.rsem$counts)
write.table(txi.rsem$counts, "RSEM_gene_counts.txt", row.names=T, col.names=T, sep="\t")
```
