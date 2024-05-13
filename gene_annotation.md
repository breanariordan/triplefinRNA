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

## PFAM

```
nano pfam.sl
```
```
**TO ADD**
```
```
sbatch pfam.sl
```

# Loading the SQLite database

```
module load Trinotate/3.2.2-GCC-9.2.0
Trinotate Trinotate.sqlite init \ --gene_trans_map RF_trinity_output.Trinity.fasta.gene_trans_map --transcript_fasta RF_trinity_output.Trinity.fasta --transdecoder_pep RF_trinity_output.Trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam PFAM.out

Trinotate Trinotate.sqlite report > Trinotate.xls
```

# Generating GO annotations

```
module load Miniconda3
source activate transdecoder
/opt/nesi/CS400_bdw/Trinity/2.14.0/trinityrnaseq/util/extract_GO_assignments_from_Trinotate.xls.pl --gene --Trinotate_xls Trinotate.xls -G --include_ancestral_terms > go_annotations.txt
conda deactivate
```
