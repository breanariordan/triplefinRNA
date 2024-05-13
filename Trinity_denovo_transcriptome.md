# Preparing the files

This file contains the code for making the _de_novo_ transcriptome using the trimmed and quality checked 2x 151 bp pair-end read triplefin data.

```
cd /nesi/nobackup/denovosourcefiles
mkdir scheduler      # make a directory for submitting all of your requests
cd scheduler
```

make a new file, named samples_file.txt (see [samples_file.txt](https://github.com/breanariordan/triplefinRNA/blob/main/samples_file.txt))

To use Trinity to make the _de_novo_ transcriptome, we submitted slurm jobs to the NeSI interface using scripts. The Trinity run was separated into two phases to make it easier and quicker.

# Trinity Phase 1

```
cd scheduler              # if not already there
nano trinityslurmp1.sl    # this creates a slurm script you can submit
```
Using the nano code will open up a script that you can then input the following commands into...

```
#!/bin/bash -e

#SBATCH --job-name=trinity-phase1
#SBATCH --account=uoo03946
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --hint=nomultithread

# load a trinity module
module purge
module load Trinity/2.14.0-gimkl-2022a

# run trinity, stop before phase 2
srun Trinity --no_distributed_trinity_exec \
  --CPU ${SLURM_CPUS_PER_TASK} --max_memory 200G \
  --seqType fq --samples_file samples_file.txt --SS_lib_type RF --output RF_trinity_output
```

Once the script is created you can send it away to be run
```
sbatch trinityslurmp1.sl
```

# Trinity Phase 2

```
nano trinityslurmp2.sl      # generate script
```

```
#!/bin/bash -e

#SBATCH --job-name=trinity-phase2grid
#SBATCH --account=uoo03946
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --hint=nomultithread
#SBATCH --partition=bigmem

# load a trinity module
module purge
module load Trinity/2.14.0-gimkl-2022a
module load HpcGridRunner/20210803

# run Trinity phase 2
srun Trinity --CPU ${SLURM_CPUS_PER_TASK} --max_memory 200G \
  --grid_exec "hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/SLURM.conf -c" \
  --seqType fq --samples_file samples_file.txt --SS_lib_type RF --output RF_trinity_output
```

```
sbatch trinityslurmp2.sl    # send script away to be run
```

After this, the next step is to obtain the iscountmatrix. The code for this can be found at [Obtaining_iscount_matrix](https://github.com/breanariordan/triplefinRNA/blob/main/Obtaining_iscount_matrix.md)
