#!/bin/bash
#SBATCH --time=01-20:00
#SBATCH --account=rrg-shapiro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=204800M
#SBATCH --cpus-per-task=2
module load  nixpkgs/16.09 gcc/7.3.0 intel/2018.3 gcc/7.3.0 blast+/2.9.0 prinseq/0.20.4 fastp/0.20.0 bwa/0.7.17 picard/2.20.6 samtools/1.10 varscan/2.4.1 python/3.6
chmod u+x WORKSPACE/*.py WORKSPACE/*.sh
./run_iPMVC_in_parallel.py WORKSPACE lst_samples.txt 2
