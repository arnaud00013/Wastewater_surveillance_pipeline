# Wastewater_surveillance_pipeline
@Author: Arnaud N'Guessan

## Overview
This repository contains a suite of tools for analyzing the within-sample genetic diversity of SARS-CoV-2 in Wastewater (WW) samples. It is designed to work on a Compute canada server but the scripts can be run on any linux system cluster with a slurm Work manager.

## Dependencies
Linux packages: nixpkgs/16.09 gcc/7.3.0 intel/2018.3 gcc/7.3.0 blast+/2.9.0 prinseq/0.20.4 fastp/0.20.0 bwa/0.7.17 picard/2.20.6 samtools/1.10 varscan/2.4.1 python/3.6
Python3.6 modules: "sys", "time", "multiprocessing" and "os"
R (version 3.5.2+) packages: "ggplot2", "seqinr", "grid", "RColorBrewer", "cowplot", "randomcoloR", "gplots", "lmPerm", "ggpubr", "gridExtra", "RColorBrewer", "tidyr", "dendextend", "VennDiagram", "Cairo", "UpSetR", "parallel", "foreach", "doParallel", "infotheo", "igraph", "glmnet", "FD", "vegan", "indicspecies", "ConsReg", "MASS", "leaps", "caret", "mgcv" and "session"

## Pipeline modules
The pipeline is separated in two main modules that sould be run in the following order:
1. Variant calling module
For unning this module of the pipeline, you need to submit a job to your cluster slurm queue using the following command 'sbatch' command and the *.sh files. 
a) Inputs (files you need to copy into the Variant_calling workspace): 
-2 Paired-end .fastq files for each sample (extension should be _R1.fastq and _R2.fastq). However, you can use single-end fastq but you have to modify the run_iPMVC_in_parallel.py script at lines 15-25 appropriately. The fastq files should be in the Variant_calling workspace and you need to specify the absolute path of the Variant_calling workspace in the 2 *.sh files where 'WORKSPACE' is indicated.
-A text file with the list of samples named "lst_samples.txt"

b) Output: VarScan .tab files, which are tab-delimited files.

2. Post-variant-calling analysis module
This module run various analysis of SARS-CoV-2 within host diversity (Sample coverage, lineage detection, estimation of the lineages within-sample frequency, etc)
a) Inputs (files you need to copy into the Post_variant_calling_analysis workspace):

