#!/bin/python3.6

#Script for individual samples pre-processing before multisample variant calling
#import libraries
import sys
import time
import multiprocessing as mp
import os

#time feedback
start_time = time.time()

#A function that execute the recurring step of this pipeline (Reads preprocessing + Mapping + Variant Calling)
def run_core_iPMVC(wp_path,current_sample, nb_t_profiling):
	#unzip fastq.gz
	os.system("gunzip -k {0}{1}_R1_001.fastq.gz".format(wp_path, current_sample))
	os.system("gunzip -k {0}{1}_R2_001.fastq.gz".format(wp_path, current_sample))
	# quality control, trimming, automatic adapters removal and low complexity reads removal
	os.system("fastp -i {0}_R1_001.fastq -I {0}_R2_001.fastq -o {0}_trimmed_1.fastq -O {0}_trimmed_2.fastq -l 70 -x --cut_tail --cut_tail_mean_quality 20 --detect_adapter_for_pe --thread {1}".format(current_sample, nb_t_profiling))
	os.system("prinseq-lite.pl -fastq {0}_trimmed_1.fastq -fastq2 {0}_trimmed_2.fastq -lc_method entropy -lc_threshold 60 -out_good {0}_trimmed_lcfiltered -out_bad {0}_trimmed_bad -out_format 1".format(current_sample))
	os.system("rm {0}{1}_R1_001.fastq".format(wp_path, current_sample))
	os.system("rm {0}{1}_R2_001.fastq".format(wp_path, current_sample))
	
	#if sample contains whole metagenomes shotgun sequences
	#print(("blastn -task megablast -db {0}Betacoronaviruses.fasta -query {1}_trimmed_lcfiltered_1.fasta -out {1}_trimmed_lcfiltered_1.xml -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -qcov_hsp_perc 60 -perc_identity 60 -outfmt 5 -num_threads {2}".format(wp_path, current_sample, nb_t_profiling)))
	#os.system("blastn -task megablast -db {0}Betacoronaviruses.fasta -query {1}_trimmed_lcfiltered_1.fasta -out {1}_trimmed_lcfiltered_1.xml -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -qcov_hsp_perc 60 -perc_identity 60 -outfmt 5 -num_threads {2}".format(wp_path, current_sample, nb_t_profiling))
	#print(("blastn -task megablast -db {0}Betacoronaviruses.fasta -query {1}_trimmed_lcfiltered_2.fasta -out {1}_trimmed_lcfiltered_2.xml -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -qcov_hsp_perc 60 -perc_identity 60 -outfmt 5 -num_threads {2}".format(wp_path, current_sample, nb_t_profiling)))
	#os.system("blastn -task megablast -db {0}Betacoronaviruses.fasta -query {1}_trimmed_lcfiltered_2.fasta -out {1}_trimmed_lcfiltered_2.xml -evalue 1e-10 -max_target_seqs 10 -max_hsps 1 -qcov_hsp_perc 60 -perc_identity 60 -outfmt 5 -num_threads {2}".format(wp_path, current_sample, nb_t_profiling))
	#print(("./extract_Betacoronavirus_reads_from_sample.py {0} {1}".format(wp_path, current_sample)))
	##os.system("./extract_Betacoronavirus_reads_from_sample.py {0} {1}".format(wp_path, current_sample))

	#index reference genome if it has not been done yet
	#print("bwa index {0}MN908947_3.fasta".format(wp_path))
	#os.system("bwa index {0}MN908947_3.fasta".format(wp_path))
	
	#Align reads to reference genome
	#print(("bwa mem {0}MN908947_3.fasta {0}{1}_trimmed_lcfiltered_onlyBetacoronavirus_1.fasta {0}{1}_trimmed_lcfiltered_onlyBetacoronavirus_2.fasta -t {2} > {0}{1}_preprocessed.sam".format(wp_path, current_sample, nb_t_profiling)))
	os.system("bwa mem {0}MN908947_3.fasta {0}{1}_trimmed_lcfiltered_1.fasta {0}{1}_trimmed_lcfiltered_2.fasta -t {2} > {0}{1}_preprocessed.sam".format(wp_path, current_sample, nb_t_profiling))
		#if you need to target Betacoronaviruses reads, map the filtered reads with the following command
	##os.system("bwa mem {0}MN908947_3.fasta {0}{1}_trimmed_lcfiltered_onlyBetacoronavirus_1.fasta {0}{1}_trimmed_lcfiltered_onlyBetacoronavirus_2.fasta -t {2} > {0}{1}_preprocessed.sam".format(wp_path, current_sample, nb_t_profiling))
	
	#Compress, sort and filter alignment
	#print("samtools view -bS {0}{1}_preprocessed.sam > {0}{1}_preprocessed.bam".format(wp_path, current_sample))
	os.system("samtools view -bS {0}{1}_preprocessed.sam > {0}{1}_preprocessed.bam".format(wp_path, current_sample))
	#print(("samtools sort {0}{1}_preprocessed.bam -o {0}{1}_preprocessed_sorted.bam".format(wp_path, current_sample)))
	os.system("samtools sort {0}{1}_preprocessed.bam -o {0}{1}_preprocessed_sorted.bam".format(wp_path, current_sample))
	#Remove PCR duplicates only if the method was shotgun sequencing
	#print(("java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={0}{1}_preprocessed_sorted.bam O={0}{1}_preprocessed_sorted_noduplicates.bam M={0}last_marked_dupl.txt REMOVE_DUPLICATES=true".format(wp_path, current_sample)))
	#os.system("java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={0}{1}_preprocessed_sorted.bam O={0}{1}_preprocessed_sorted_noduplicates.bam M={0}last_marked_dupl.txt REMOVE_DUPLICATES=true".format(wp_path, current_sample))
	#print(("samtools mpileup -B -f {0}MN908947_3.fasta -d 0 {0}{1}_preprocessed_sorted_noduplicates.bam > {0}{1}.pileup".format(wp_path, current_sample)))
	#os.system("samtools mpileup -B -f {0}MN908947_3.fasta {0}{1}_preprocessed_sorted_noduplicates.bam > {0}{1}.pileup".format(wp_path, current_sample))
	os.system("samtools mpileup -B -f {0}MN908947_3.fasta {0}{1}_preprocessed_sorted.bam > {0}{1}.pileup".format(wp_path, current_sample))
	#print(("java -jar $EBROOTVARSCAN/VarScan.v2.4.1.jar pileup2snp {0}{1}.pileup --min-coverage 400 --min-reads2 5 --min-var-freq 0.03 > {0}variants_{1}.tab".format(wp_path, current_sample)))
	os.system("java -jar $EBROOTVARSCAN/VarScan.v2.4.1.jar pileup2snp {0}{1}.pileup --min-coverage 100 --min-reads2 5 --min-var-freq 0.05 > {0}variants_{1}.tab".format(wp_path, current_sample))
	#Add sample name in Already_analyzed_samples.txt
	print("echo {1} >> {0}Already_analyzed_samples.txt".format(wp_path, current_sample))
	os.system("echo {1} >> {0}Already_analyzed_samples.txt".format(wp_path, current_sample))

#create a function that test if a string is numeric
def is_numeric(the_str):
	try:
		float_conv = float(the_str)
		return True
	except (ValueError, TypeError):
		return False

if __name__ == "__main__":
	#Initializing global variables
	#Path of the Workspace where there are all samples files
	workspace_path =sys.argv[1] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((workspace_path[len(workspace_path)-1]) != "/"):
		workspace_path = "{0}/".format(workspace_path)
	print("Chosen Workspace is : '{0}'".format(workspace_path))

	# create the variable lst_samples that will contain the list of the samples in the workspace
	samples_list_file = sys.argv[2]
	lst_samples=[]
	with open("{0}{1}".format(workspace_path,samples_list_file)) as f:
		lst_samples = f.readlines()

	the_idx = 0
	for sample_name in lst_samples:
		lst_samples[the_idx] = sample_name.replace("\n", "")
		the_idx = the_idx + 1
	#Number of samples being analysed simultaneously
	nb_sim_process = 8

	#Number of threads to use for multithread tasks (fastp and bwa)
	nb_threads = int(sys.argv[3])

	#create file that records samples that have already been analyzed
	os.system("/bin/rm -f {0}Already_analyzed_samples.txt".format(workspace_path))
	os.system("/bin/touch {0}Already_analyzed_samples.txt".format(workspace_path))

	#execute nb_sim_process samples at a time
	liste_index_sample = range(0,len(lst_samples),nb_sim_process)
	for i in liste_index_sample:
		# create a list that will record the samples that have already been analysed
		already_analysed_samples_list=[]
		with open("{0}Already_analyzed_samples.txt".format(workspace_path)) as f:
			already_analyzed_samples_list = f.readlines()
		index_line = 0
		for line in already_analyzed_samples_list:
			already_analyzed_samples_list[index_line] = line.replace("\n", "")
			index_line = index_line + 1

		processes = []
		#define all processes that will be executed parallelly (run only those that have not been executed yet)
		if (((i)<len(lst_samples)) and (not (lst_samples[i] in already_analyzed_samples_list))):
			p1 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i], nb_threads,))
			processes.append(p1)
			p1.start()
		if (((i+1)<len(lst_samples)) and (not (lst_samples[i+1] in already_analyzed_samples_list))):
			p2 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+1], nb_threads,))
			processes.append(p2)
			p2.start()
		if (((i+2)<len(lst_samples)) and (not (lst_samples[i+2] in already_analyzed_samples_list))):
			p3 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+2], nb_threads,))
			processes.append(p3)
			p3.start()
		if (((i+3)<len(lst_samples)) and (not (lst_samples[i+3] in already_analyzed_samples_list))):
			p4 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+3], nb_threads,))
			processes.append(p4)
			p4.start()
		if (((i+4)<len(lst_samples)) and (not (lst_samples[i+4] in already_analyzed_samples_list))):
			p5 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+4], nb_threads,))
			processes.append(p5)
			p5.start()
		if (((i+5)<len(lst_samples)) and (not (lst_samples[i+5] in already_analyzed_samples_list))):
			p6 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+5], nb_threads,))
			processes.append(p6)
			p6.start()
		if (((i+6)<len(lst_samples)) and (not (lst_samples[i+6] in already_analyzed_samples_list))):
			p7 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+6], nb_threads,))
			processes.append(p7)
			p7.start()
		if (((i+7)<len(lst_samples)) and (not (lst_samples[i+7] in already_analyzed_samples_list))):
			p8 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+7], nb_threads,))
			processes.append(p8)
			p8.start()
		#if (((i+8)<len(lst_samples)) and (not (lst_samples[i+8] in already_analyzed_samples_list))):
		 	#p9 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+8], nb_threads,))
		 	#processes.append(p9)
		 	#p9.start()
		#if (((i+9)<len(lst_samples)) and (not (lst_samples[i+9] in already_analyzed_samples_list))):
		 	#p10 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+9], nb_threads,))
		 	#processes.append(p10)
		 	#p10.start()
		#if (((i+10)<len(lst_samples)) and (not (lst_samples[i+10] in already_analyzed_samples_list))):
		 	#p11 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+10], nb_threads,))
		 	#processes.append(p11)
		 	#p11.start()
		#if (((i+11)<len(lst_samples)) and (not (lst_samples[i+11] in already_analyzed_samples_list))):
		 	#p12 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+11], nb_threads,))
		 	#processes.append(p12)
			#p12.start()
		#if (((i+12)<len(lst_samples)) and (not (lst_samples[i+12] in already_analyzed_samples_list))):
		 	#p13 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+12], nb_threads,))
		 	#processes.append(p13)
			#p13.start()
		#if (((i+13)<len(lst_samples)) and (not (lst_samples[i+13] in already_analyzed_samples_list))):
		 	#p14 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+13], nb_threads,))
		 	#processes.append(p14)
		 	#p14.start()
		#if (((i+14)<len(lst_samples)) and (not (lst_samples[i+14] in already_analyzed_samples_list))):
		 	#p15 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+14], nb_threads,))
		 	#processes.append(p15)
		 	#p15.start()
		#if (((i+15)<len(lst_samples)) and (not (lst_samples[i+15] in already_analyzed_samples_list))):
		 	#p16 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+15], nb_threads,))
		 	#processes.append(p16)
		 	#p16.start()
		# if (((i+16)<len(lst_samples)) and (not (lst_samples[i+16] in already_analyzed_samples_list))):
		# 	p17 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+16], nb_threads,))
		# 	processes.append(p17)
		# 	p17.start()
		# if (((i+17)<len(lst_samples)) and (not (lst_samples[i+17] in already_analyzed_samples_list))):
		# 	p18 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+17], nb_threads,))
		# 	processes.append(p18)
		# 	p18.start()
		# if (((i+18)<len(lst_samples)) and (not (lst_samples[i+18] in already_analyzed_samples_list))):
		# 	p19 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+18], nb_threads,))
		# 	processes.append(p19)
		# 	p19.start()
		# if (((i+19)<len(lst_samples)) and (not (lst_samples[i+19] in already_analyzed_samples_list))):
		# 	p20 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+19], nb_threads,))
		# 	processes.append(p20)
		# 	p20.start()
		#wait for all processes to finish
		if (len(processes) != 0):
			for p in processes:
				p.join()
			for p in processes:
				while (p.is_alive()):
					time.sleep(5)

	print("--- Pipeline Execution started at %s ---" % (start_time))
	print("--- Pipeline Execution finished at %s ---" % (time.time()))
	print("--- Number of samples analysed : %s ---" % (len(lst_samples)))
	print("--- Runtime for the pipeline with this machine is {0} seconds using {1} processors only for the preprocessing step!!! ---".format(time.time() - start_time, nb_threads))

