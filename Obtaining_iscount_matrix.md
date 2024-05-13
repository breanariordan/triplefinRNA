# Getting basic assembly metrics
You should now have a folder labelled RF_trinity_output (or whatever you named your output file) as well as a general output file RF_trinity_output.Trinity.fasta

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
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=30G
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
  --samples_file samples_file.txt \
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
ln -s  ../10L2brain/RSEM.genes.results 10L2brain.txt
ln -s  ../10N2heart/RSEM.genes.results 10N2heart.txt
ln -s  ../10N3brain/RSEM.genes.results 10N3brain.txt
ln -s  ../22L1heart/RSEM.genes.results 22L1heart.txt
ln -s  ../22L4brain/RSEM.genes.results 22L4brain.txt
ln -s  ../22N1brain/RSEM.genes.results 22N1brain.txt
ln -s  ../22N1heart/RSEM.genes.results 22N1heart.txt
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
