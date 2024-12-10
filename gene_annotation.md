# Transdecoder annotation
Download the transdecoder from online if you need to (https://github.com/TransDecoder/TransDecoder/wiki)
I would recommend first loading up all of the modules you require, such as follows:

```
module load Trinotate/3.2.2-GCC-9.2.0
module load Miniconda3
module load Infernal
module load HMMER
module load BLAST
module load Trinity
module load SQLite
module load Infernal
module load Perl
```

If you do not wish to do this, I have included most of the individual modules required for each step.

```
cd /nobackup/denovosourcefiles/scheduler
module load Miniconda3
source activate transdecoder  	# self install
TransDecoder.Long0rfs -t RF_trinity_output.Trinity.fasta
TransDecoder.Predict -t RF_trinity_output.Trinity.fasta
conda deactivate
```

# Creating the SQLite Database
Note that for this to work you **must** use version 3.2.2 of Trinotate

```
module load Trinotate/3.2.2-GCC-9.2.0
Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
```

# Adding BLAST hits
## BlastX

```
nano blastx.sl    # create slurm script
```
```
#!/bin/bash -e

#SBATCH --job-name=trinotate-blast
#SBATCH --account=uoo03946
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --hint=nomultithread

module load Miniconda3
source activate transdecoder
module load BLAST

srun blastx -query RF_trinity_output.Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
```
```
sbatch blastx.sl    # submit slurm script
```

## BlastP

```
nano blastp.sl    # create slurm script
```
```
#!/bin/bash -e

#SBATCH --job-name=trinotate-blast
#SBATCH --account=uoo03946
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --hint=nomultithread

module load Miniconda3
source activate transdecoder
module load BLAST

srun blastp -query RF_trinity_output.Trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
```
```
sbatch blastp.sl    # submit slurm script
```

# Loading the SQLite database

```
module load Trinotate/3.2.2-GCC-9.2.0

Trinotate Trinotate.sqlite init \ --gene_trans_map RF_trinity_output.Trinity.fasta.gene_trans_map --transcript_fasta RF_trinity_output.Trinity.fasta --transdecoder_pep RF_trinity_output.Trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6

Trinotate Trinotate.sqlite report > Trinotate.xls
```

# Generating GO annotations

```
module load Miniconda3
source activate transdecoder
/opt/nesi/CS400_centos7_bdw/Trinotate/3.2.2-GCC-9.2.0/util/extract_GO_assignments_from_Trinotate.xls.pl --gene --Trinotate_xls Trinotate.xls -G --include_ancestral_terms > go_annotations.txt
conda deactivate
```
# Extract gene lengths for samples
```
/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022b/trinityrnaseq-Trinity-v2.14.0/util/misc/fasta_seq_length.pl RF_trinity_output.Trinity.fasta > Trinity.fasta.seq_lens
```
# Aligning samples using salmon 
This line of code specifically works to generate a matrix of gene abundance estimates for all of the samples
```
/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/util/align_and_estimate_abundance.pl --transcripts RF_trinity_output.Trinity.fasta --seqType fq --samples_file samples_file.txt --est_method salmon --trinity_mode --prep_reference > salmon_align_and_estimate_abundance.log 2>&1 &

/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl \
> --est_method salmon \
> --name_sample_by_basedir \
> --gene_trans_map none \
> 10L1brain/quant.sf \
> 10L1heart/quant.sf \
> 10L2brain/quant.sf \
> 10L2heart/quant.sf \
> 10L3brain/quant.sf \
> 10L3heart/quant.sf \
> 10L4brain/quant.sf \
> 10L4heart/quant.sf \
> 10L5brain/quant.sf \
> 10L5heart/quant.sf \
> 10N1brain/quant.sf \
> 10N1heart/quant.sf \
> 10N2brain/quant.sf \
> 10N2heart/quant.sf \
> 10N3brain/quant.sf \
> 10N3heart/quant.sf \
> 10N4brain/quant.sf \
> 10N4heart/quant.sf \
> 10N5brain/quant.sf \
> 10N5heart/quant.sf \
> 22L1brain/quant.sf \
> 22L1heart/quant.sf \
> 22L2brain/quant.sf \
> 22L2heart/quant.sf \
> 22L3brain/quant.sf \
> 22L3heart/quant.sf \
> 22L4brain/quant.sf \
> 22L4heart/quant.sf \
> 22L5brain/quant.sf \
> 22L5heart/quant.sf \
> 22N1brain/quant.sf \
> 22N1heart/quant.sf \
> 22N2brain/quant.sf \
> 22N2heart/quant.sf \
> 22N3brain/quant.sf \
> 22N3heart/quant.sf \
> 22N4brain/quant.sf \
> 22N4heart/quant.sf \
> 22N5brain/quant.sf \
> 22N5heart/quant.sf
```
# Incorporating gene lengths into the matrix
```
/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/util/misc/TPM_weighted_gene_length.py --gene_trans_map RF_trinity_output.Trinity.fasta.gene_trans_map --trans_lengths Trinity.fasta.seq_lens --TPM_matrix salmon.isoform.TMM.EXPR.matrix > Trinity.gene_lengths.txt
```
# Factor labelling using DE genes
***Note***, the DE files you require for the following steps are generated in [DE_Analysis](https://github.com/breanariordan/triplefinRNA/blob/main/DE_Analysis.md). This step is performed separately for the two species with files put into separate folders created for _F. nigripenne_ and _F. lapillum_. **This example is for the lapillum brain samples**. 
```
wc -l lapbrain_DE_results.txt
cut -f 1 lapbrain_DE_results.txt | tail -n 7616 | awk '{print "diff\t",$0}' | sed -s 's/\"//g' > factor_labelinglapbrain.txt
```
Download the relevant modules and packages
```
module load Miniconda3

conda install -c bioconda bioconductor-goseq

module load R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16", lib="/home/riobr412/R/gimkl-2022a/4.2")

BiocManager::install("qvalue", lib="/home/riobr412/R/gimkl-2022a/4.2", force = TRUE)
BiocManager::install("goseq", lib="/home/riobr412/R/gimkl-2022a/4.2", force = TRUE)
```
Perform all factor labelling. I moved all generated files into the folder of the designated species.
```
/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling DE_files/lapbrain/factor_labelinglapbrain.txt --GO_assignments go_annotations.txt --lengths Trinity.gene_lengths.txt --background DE_files/lapbrain/GOlapbrain_background.txt

# Upregulated genes

cat ~/nobackup/denovosourcefiles/scheduler/DE_files/lapbrain/GOlapbrain_upregulated.txt |  awk '{print "diff\t",$0}' > factor_labelingupregulated.txt

/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling DE_files/lapbrain/factor_labelingupregulated.txt --GO_assignments go_annotations.txt --lengths Trinity.gene_lengths.txt --background DE_files/lapbrain/GOlapbrain_background.txt

# Downregulated genes

cat ~/nobackup/denovosourcefiles/scheduler/DE_files/lapbrain/GOlapbrain_downregulated.txt |  awk '{print "diff\t",$0}' > factor_labelingdownregulated.txt

/opt/nesi/CS400_centos7_bdw/Trinity/2.14.0-gimkl-2022a/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling DE_files/lapbrain/factor_labelingdownregulated.txt --GO_assignments go_annotations.txt --lengths Trinity.gene_lengths.txt --background DE_files/lapbrain/GOlapbrain_background.txt
```
