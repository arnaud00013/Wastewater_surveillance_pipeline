#!/bin/bash
#SBATCH --time=01-00:00
#SBATCH --account=rrg-shapiro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20480M
#SBATCH --cpus-per-task=10
chmod u+x /home/p1211536/scratch/Wastewater/Illumina/*.py /home/p1211536/scratch/Wastewater/Illumina/*.r /home/p1211536/scratch/Wastewater/Illumina/*.sh
module load nixpkgs/16.09 gcc/7.3.0 nixpkgs/16.09 intel/2018.3 samtools/1.10 varscan/2.4.1
module load samtools/1.10
ls *preprocessed_sorted.bam > lst_samples_bamfiles.txt
#depth report
samtools depth *preprocessed_sorted.bam -o common_depth_report.csv
