For the longer 2x 151 bp pair-end reads used to generate the de novo transcriptome, the code used to trim the data is different to that used for shorter paired-end reads. To see the code for the shorter paired-end reads, see the file [Quality_control_and_trimming](https://github.com/breanariordan/triplefinRNA/blob/main/Quality_control_and_trimming.md)

# Moving raw files

```
pwd                             # past working directory, find where you are
cd OG7591-691868150             # identify/move to to directory where raw RNASeq files are
cd ../nobackup                  # move back to nobackup folder
mkdir denovosourcefiles         # make a new directory called 'denovosourcefiles'
cd denovosourcefiles
mkdir denovorun                 # make a new directory for the raw reads
cd ..                           # move back to nobackup directory
cp ../OG7591-691868150/AA*/*gz denovosourcefiles/denovorun/ # take all the files containing the letters/numbers in the asterix and copy them into the new denovorun directory
```

# Quality control

```
cd denovosourcefiles          # if you aren't already in it
mkdir QC                      # make a folder called QC
cd denovorun
module load FastQC/0.12.1     # load fastqc for quality control
fastqc -o QC/ denovorun/*     # perform quality control on the raw data
cd ../QC                      # move into the QC directory
for filename in *.zip         # unzip fastqc text file outputs
    do
    unzip $filename
    done
cat */summary.txt > ~/nobackup.sourcefiles/QC/fastqc_summaries.txt # generate summary of QC stat tests to see where failures have occurred
module load MultiQC/1.13-gimkl-2022a-Python-3.10.5 # load MultiQC
cd ..                         # move back to denovosourcefiles
mkdir MultiQC                 # make a directory called 'MultiQC'
cd MultiQC
cp ../QC/* ./                 # move all QC files for individual samples into MultiQC directory
multiqc .                     # perform multiqc analysis 
```

# Trimming
this will use the cutadapt programme to trim low quality reads out of the files and put them into a new folder called 'Trimmed'
note that as this is paired-end reads R1 and R2 represent the two directional runs of the reads

```
pwd
cd ..                         # move back to denovosourcefiles directory
mkdir Trimmed                 # make directory 'Trimmed'
module load cutadapt/4.1-gimkl-2022a-Python-3.10.5    # load latest cutadapt version
cd denovorun                  # go to the raw data folder
for filename in *R1_001.fastq.gz
    do
    echo $filename
    base=$(basename ${filename} R1_100.fastq.gz)
    cutadapt -j 4 -m 20 -q 20 -a AGATCGGAAGAG -A AGATCGGAAGAG -o ../Trimmed/${base}R1_100.fastq.gz -p ../Trimmed/${base}R2_001.fastq.gz ${base}R1_001.fastq.gz ${base}R2_001.fastq.gz
    done
cd ../MultiQC                 # return to the MultiQC folder
cp ../Trimmed/*fastq.gz .     # take all the trimmed log files
multiqc .                     # and run the multiqc again
```
