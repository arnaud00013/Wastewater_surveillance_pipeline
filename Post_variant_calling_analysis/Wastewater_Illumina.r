#@Author=Arnaud NG
#This script analyses wastewater samples (Coverage, sensistivity of detection, concordance with clinical data (TO UPDATE), etc)

#import libraries
library("ggplot2")
library("seqinr")
library("grid")
library("RColorBrewer")
library("cowplot")
library("randomcoloR")
library("gplots")
library("lmPerm")
library("ggpubr")
library("gridExtra")
library("RColorBrewer")
library("tidyr")
library("dendextend")
library("VennDiagram")
library("Cairo")
library("UpSetR")
library("parallel")
library("foreach")
library("doParallel")
library("infotheo")
library("glmnet")
library("FD")
library("vegan")
library("ConsReg")
library("MASS")
library("leaps")
library("caret")
library("mgcv")
library("session")
#import script arguments
#workspace with all the .bam and .tab files
output_workspace <- as.character(commandArgs(TRUE)[1]) #TO UPDATE
#name of the depth report file 
depth_report_filename <- "common_depth_report.csv"
#name of the reference genome fasta file
fasta_refseq_filename <- "MN908947_3.fasta"
#Number of cores to use during parallel processing
nb_cores <- 16 #TO UPDATE

palette_mutations_of_interest <- RColorBrewer::brewer.pal(length(c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S")),"Set1")
names(palette_mutations_of_interest) <- c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S")

v_lineages_of_interest <- c("B.1.1.7","B.1.351","P.1","B.1.427_and_B.1.429","B.1.160","B.1.177","B.1.617.X","B.1.525","B.1.526","C.37","P.3","P.2","A.2.5","B.1.1.318","B.1.1.519","B.1.466.2","B.1.621","B.1.214.2","AV.1","AT.1","C.36.3","R.1","R.2")
v_lineages_of_interest_with_who_desgnation <- c("Alpha (B.1.1.7)","Beta (B.1.351)","Gamma (P.1)","Epsilon (B.1.427_and_B.1.429)","B.1.160","B.1.177","Kappa+Delta (B.1.617.X)","Eta (B.1.525)","Iota (B.1.526)","Lambda (C.37)","Theta (P.3)","Zeta (P.2)","A.2.5","B.1.1.318","B.1.1.519","B.1.466.2","B.1.621","B.1.214.2","AV.1","AT.1","C.36.3","R.1","R.2")
names(v_lineages_of_interest_with_who_desgnation) <- v_lineages_of_interest

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(208)
color <- sample(col_vector, length(v_lineages_of_interest))
pie(1:length(color), col=color,labels = color)
palette_PANGO_lineages_of_interest <- sample(color,length(v_lineages_of_interest))
names(palette_PANGO_lineages_of_interest) <- v_lineages_of_interest
palette_PANGO_lineages_of_interest_with_who_designation <- palette_PANGO_lineages_of_interest
names(palette_PANGO_lineages_of_interest_with_who_designation) <- v_lineages_of_interest_with_who_desgnation[names(palette_PANGO_lineages_of_interest)]

#Set language as English for date formatting
Sys.setlocale("LC_ALL","English")

#Get list of samples
lst_samples <- unname(sapply(X = (sapply(X = read.csv2(file = paste0(output_workspace,"lst_samples.txt"),sep = "\t",header = FALSE,stringsAsFactors = FALSE)[1], FUN=function(x) gsub(pattern = output_workspace,replacement = "",x,fixed = TRUE))),FUN=function(x) gsub(pattern = "_preprocessed_sorted.bam",replacement = "",x,fixed = TRUE)))
lst_samples_original <- lst_samples
#calculate total number of samples
nb_samples <- length(lst_samples)
nb_samples_original <- length(lst_samples_original)
threshold_nb_individuals_sharing_mut <- ceiling(nb_samples/10)
#import reference fasta 
genome_refseq <- seqinr::getSequence(object = toupper(read.fasta(paste0(output_workspace,fasta_refseq_filename),seqtype = "DNA",as.string = TRUE,forceDNAtolower = FALSE)),as.string = TRUE)[[1]]

#import and concatenate samples variants calling file
df_variants <- read.csv2(file = paste0(output_workspace,"variants_",lst_samples[1],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
if (nrow(df_variants)!=0){
  df_variants$Sample <- lst_samples[1]
}
for (i in 2:nb_samples){
  df_to_add_to_variants_df <- read.csv2(file = paste0(output_workspace,"variants_",lst_samples[i],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  if (nrow(df_to_add_to_variants_df)==0){
    next()
  }
  df_to_add_to_variants_df$Sample <- lst_samples[i]
  df_variants <- rbind(df_variants,df_to_add_to_variants_df)
  
}
#Do not confuse nucleotide T with the alias for the boolean value TRUE
df_variants$Ref <- gsub(pattern = "TRUE", replacement = "T",x = df_variants$Ref,fixed=TRUE) 
df_variants$VarAllele <- gsub(pattern = "TRUE", replacement = "T",x = df_variants$VarAllele,fixed=TRUE) 
df_variants$VarAllele <- gsub(pattern = "+", replacement = "",x = df_variants$VarAllele,fixed=TRUE) 
df_variants$VarAllele <- gsub(pattern = "-", replacement = "",x = df_variants$VarAllele,fixed=TRUE)
df_variants <- subset(df_variants,nchar(VarAllele)==1)

df_variants <- subset(df_variants,subset = !duplicated(paste0(Position,VarAllele,Sample,sep="")))
#make sure that the variant frequency filter is applied (freq variant >=0.05 on at least one strand)
df_variants$VarFreq <- (as.numeric(gsub(pattern = "%",replacement = "",x = df_variants$VarFreq,fixed = TRUE))/100)
df_variants <- subset(df_variants, (Reads1+Reads2>=50)&(VarFreq>=0.25))
#apply strand bias filter on variants (freq variant >=0.02 on each strand)
df_variants <- df_variants[(((df_variants$Reads2Plus/(df_variants$Reads1Plus+df_variants$Reads2Plus))>=0.02)&((df_variants$Reads2Minus/(df_variants$Reads1Minus+df_variants$Reads2Minus))>=0.02)),]


#Eliminate site with no variants
df_variants <- subset(df_variants,!is.na(VarAllele))
#remove control samples 
df_variants <- subset(df_variants,unname(vapply(X = toupper(df_variants$Sample),FUN = function(x) !grepl(pattern = "CTRL",x = x,fixed = T),FUN.VALUE = c(F))))
df_variants <- subset(df_variants,unname(vapply(X = toupper(df_variants$Sample),FUN = function(x) !grepl(pattern = "CTL",x = x,fixed = T),FUN.VALUE = c(F))))
df_variants <- subset(df_variants,unname(vapply(X = toupper(df_variants$Sample),FUN = function(x) !grepl(pattern = "CONTROL",x = x,fixed = T),FUN.VALUE = c(F))))
#Get list of genomic region and positions
v_orfs <- c("5'UTR", "orf1a", "orf1b", "S","ORF3a","ORF3b","ORF3c","E","M","ORF6","ORF7a", "ORF7b","ORF8", "N", "ORF9c","ORF10","3'UTR")
v_start_orfs <- c(1, 266, 13468, 21563, 25393, 25814, 25524,26245, 26523, 27202, 27394, 27756,27894, 28274, 28734,29558, 29675)
names(v_start_orfs) <- v_orfs
v_end_orfs <- c(265, 13483, 21555, 25384, 26220, 25882, 25697, 26472, 27191, 27387, 27759, 27887,28259, 29533, 28955,29674, 29903)
names(v_end_orfs) <- v_orfs
find_ORF_of_mutation <- function(the_site_position){
  indx <- which((v_start_orfs<=the_site_position)&(v_end_orfs>=the_site_position))[1]
  if (length(indx)==0){
    return(NA)
  }else{
    return(v_orfs[indx])
  }
}
v_orfs_length <- v_end_orfs - v_start_orfs + 1 
palette_orfs_epitopes <- c("orf1a"="red","orf1b"="blue","S"="green3","E"="orange","M"="grey","N"="purple")

v_genes_with_unique_product <- c(paste0("NSP",1:10),paste0("NSP",12:16), "S","ORF3a","ORF3b","ORF3c","E","M","ORF6","ORF7a", "ORF7b","ORF8", "N", "ORF9c", "ORF10")
v_start_genes <- c(266,806,2720,8555,10055,10973,11843,12092,12686,13025,13442,16237,18040,19621,20659,21563, 25393, 25814, 25524, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 28734, 29558)
names(v_start_genes) <- v_genes_with_unique_product
v_end_genes <- c(805,2719,8554,10054,10972,11842,12091,12685,13024,13441,16237,18039,19620,20658,21552,25384, 26220,25882, 25697, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 28955, 29674)
names(v_end_genes) <- v_genes_with_unique_product
find_gene_of_mutation <- function(the_site_position){
  indx <- which((v_start_genes<=the_site_position)&(v_end_genes>=the_site_position))[1]
  if (length(indx)==0){
    return(NA)
  }else{
    return(v_genes_with_unique_product[indx])
  }
}

v_genes_length <- v_end_genes - v_start_genes + 1

#import metdata to analyze stratifications of interest
df_sample_stratifications_of_interest <- read.csv2(file = paste0(output_workspace,"Table_Sample_stratifications_of_interest.csv"),sep = ";",header = TRUE,stringsAsFactors = FALSE)
df_sample_stratifications_of_interest$X <- NULL
colnames(df_sample_stratifications_of_interest)[1] <- "Sample"
df_sample_stratifications_of_interest$date <- ifelse(test = df_sample_stratifications_of_interest$date=="",yes=NA,no=df_sample_stratifications_of_interest$date)
for (i in (1:nrow(df_sample_stratifications_of_interest))){
  if ((!is.na(df_sample_stratifications_of_interest$date[i]))){
    v_pos_character_delim <- gregexpr(pattern = "-",text = df_sample_stratifications_of_interest$date[i],fixed = T)[[1]]
    first_term <- substr(df_sample_stratifications_of_interest$date[i],1,v_pos_character_delim[1]-1)
    if (nchar(first_term)==4){
      next()
    }
    second_term <- substr(df_sample_stratifications_of_interest$date[i],v_pos_character_delim[1]+1,v_pos_character_delim[2]-1)
    third_term <- substr(df_sample_stratifications_of_interest$date[i],v_pos_character_delim[2]+1,nchar(df_sample_stratifications_of_interest$date[i]))
    df_sample_stratifications_of_interest$date[i] <- paste0(third_term,"-",second_term,"-",first_term)
  }
}
get_location_from_sample_name_with_location <- function(the_sample_name){
  return(toupper(subset(df_sample_stratifications_of_interest,Sample==the_sample_name)$City[1]))
}
#get sample site
get_site_of_sample <- function(the_sample_name){
  return(subset(df_sample_stratifications_of_interest,Sample==the_sample_name)$site_name[1])
}

#get the type of site with the sample name
get_site_type_of_sample <- function(the_sample_name){
  return(subset(df_sample_stratifications_of_interest,Sample==the_sample_name)$Site_type[1])
}
#get the sample location in the network
get_Location_in_network_of_sample <- function(the_sample_name){
  return(subset(df_sample_stratifications_of_interest,Sample==the_sample_name)$Location_in_network[1])
}
#get the type of colelctor with the sample name
get_collector_type_of_sample <- function(the_sample_name){
  return(subset(df_sample_stratifications_of_interest,Sample==the_sample_name)$Collector_type[1])
}
df_variants$location <- unname(vapply(X = df_variants$Sample,FUN = get_location_from_sample_name_with_location,FUN.VALUE = c("")))
df_variants$location <- ifelse(test=df_variants$location=="",yes=NA,no=df_variants$location)
df_variants <- subset(df_variants,(location!="ON")&(location!="OTTAWA"&(location!="ONT")))
df_variants$ORF <- vapply(X = df_variants$Position,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_variants$gene <- vapply(X = df_variants$Position,FUN = function(x) return(find_gene_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
#duplicate variants of orf3a if they occur also in orf3b or orf3c
df_subset_orf3b_orf3c <- subset(df_variants,subset=(Position>=v_start_orfs["ORF3a"])&(Position<=v_end_orfs["ORF3a"]))
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"])),"ORF"] <- "ORF3b"
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"])),"ORF"] <- "ORF3c"
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"])),"gene"] <- "ORF3b"
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"])),"gene"] <- "ORF3c"
df_subset_orf3b_orf3c <- subset(df_subset_orf3b_orf3c,ORF!="ORF3a")
df_variants <- rbind(df_variants,df_subset_orf3b_orf3c)
#duplicate variants of ORF7a if they occur also in ORF7b
df_subset_orf7b <- subset(df_variants,subset=(Position>=v_start_orfs["ORF7a"])&(Position<=v_end_orfs["ORF7a"]))
df_subset_orf7b[which((df_subset_orf7b$Position>=v_start_orfs["ORF7b"])&(df_subset_orf7b$Position<=v_end_orfs["ORF7b"])),"ORF"] <- "ORF7b"
df_subset_orf7b <- subset(df_subset_orf7b,ORF!="ORF7a")
df_variants <- rbind(df_variants,df_subset_orf7b)
#duplicate variants of N if they occur also in orf9c
df_subset_ORF9c <- subset(df_variants,subset=(Position>=v_start_orfs["N"])&(Position<=v_end_orfs["N"]))
df_subset_ORF9c[which((df_subset_ORF9c$Position>=v_start_orfs["ORF9c"])&(df_subset_ORF9c$Position<=v_end_orfs["ORF9c"])),"ORF"] <- "ORF9c"
df_subset_ORF9c <- subset(df_subset_ORF9c,ORF!="N")
df_variants <- rbind(df_variants,df_subset_ORF9c)
#Check if CODING regions (not UTR!!!) length are a multiple of 3
#vapply(X = v_orfs,FUN = function(x) return((((v_end_orfs[x]-v_start_orfs[x])+1)%%3)==0),FUN.VALUE = c(TRUE))

#color palettes
set.seed(1234)
palette_orfs <- distinctColorPalette(length(v_orfs))
names(palette_orfs) <- v_orfs
# pie(rep(1, length(v_orfs)), col=palette_orfs)
palette_genes <- distinctColorPalette(length(v_genes_with_unique_product))
names(palette_genes) <- v_genes_with_unique_product
# pie(rep(1, length(v_genes_with_unique_product)), col=palette_genes)

#function to get mutation in genomic format
find_candidate_genomic_mutation_causing_aa_change_in_orf <- function(current_mutation_name){
  the_orf <- strsplit(x = current_mutation_name,split = ":")[[1]][1]
  the_mut <-strsplit(x = current_mutation_name,split = ":")[[1]][2]
  old_aa <- substr(the_mut,1,1)
  new_aa <- substr(the_mut,nchar(the_mut),nchar(the_mut))
  pos_in_orf_prot_seq <- as.integer(substr(the_mut,2,nchar(the_mut)-1))
  pos_start_codon_in_orf <- ((pos_in_orf_prot_seq-1)*3)+1
  pos_middle_codon_in_orf <- pos_start_codon_in_orf + 1
  pos_stop_codon_in_orf <- pos_start_codon_in_orf + 2
  pos_start_codon_in_genome <- v_start_orfs[the_orf] + pos_start_codon_in_orf - 1
  pos_middle_codon_in_genome <- v_start_orfs[the_orf] + pos_middle_codon_in_orf - 1
  pos_stop_codon_in_genome <- v_start_orfs[the_orf] + pos_stop_codon_in_orf - 1
  ref_codon <- substr(genome_refseq,pos_start_codon_in_genome,pos_stop_codon_in_genome)
  if (old_aa!=translate_seq(the_codon = ref_codon)){
    stop("Logical error: extracted reference codon does not translate into old aa!")
  }
  v_candidate_genomic_mutations <- NULL
  for (pos_in_codon in 1:3){
    for (current_new_nucl in setdiff(c("A","T","C","G"),substr(ref_codon,pos_in_codon,pos_in_codon))){
      if (pos_in_codon ==1){
        new_codon <- paste0(current_new_nucl,substr(ref_codon,2,3))
      }else if (pos_in_codon ==2){
        new_codon <- paste0(substr(ref_codon,1,1),current_new_nucl,substr(ref_codon,3,3))
      }else if (pos_in_codon ==3){
        new_codon <- paste0(substr(ref_codon,1,1),substr(ref_codon,2,2),current_new_nucl)
      }else{
        stop("Logical error: position in codon cannot be different than 1, 2 or 3!")
      }
      if (translate_seq(new_codon)==new_aa){
        v_candidate_genomic_mutations <- c(v_candidate_genomic_mutations,paste0(substr(ref_codon,pos_in_codon,pos_in_codon),pos_start_codon_in_genome+pos_in_codon-1,current_new_nucl))
      }
    }
  }
  return(v_candidate_genomic_mutations)
}

#compute Shannon entropy
get_entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

#function that find the original and mutated codons of a variant
get_ref_and_mutated_codon <- function(the_position,ref_nucl,new_nucl){
  the_orf <- find_ORF_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    the_ref_codon <- NA
    the_mut_codon <- NA
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      the_ref_codon <- paste0(ref_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2),sep="")
      the_mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2),sep="")
    }else if (pos_in_codon==2){
      the_ref_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),ref_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1),sep="")
      the_mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1),sep="")
    }else if (pos_in_codon==3){
      the_ref_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),ref_nucl,sep="")
      the_mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl,sep="")
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  return(list(ref_codon=the_ref_codon,mutated_codon=the_mut_codon))
}
#build function that determines whether a mutation is synonymous or not 
is_mutation_synonymous <- function(the_reference_codon,the_mutated_codon){
  if (the_reference_codon %in% c("TAA","TAG","TGA")){
    return(NA)
  }else{
    return(seqinr::translate(seq = unlist(strsplit(the_reference_codon,"")))==seqinr::translate(seq = unlist(strsplit(the_mutated_codon,""))))
  }
}
#build function that determines whether a mutation is synonymous or not 
translate_seq <- function(the_codon){
  if (is.na(the_codon)){
    return(NA)
  }else if (the_codon %in% c("TAA","TAG","TGA")){
    return("Stop")
  }else{
    return(seqinr::translate(seq = unlist(strsplit(the_codon,""))))
  }
}

#function for plotting linear model
ggplotRegression <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit)$coefficients[,4][2]), format = "e", digits = 3)))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P: ",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm")
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = unname(summary(fit)$coefficients[,4][2])))
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}
ggplotRegression_export_eps <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit)$coefficients[,4][2]), format = "e", digits = 3)))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P: ",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm", device = cairo_ps)
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = unname(summary(fit)$coefficients[,4][2])))
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}

# get_amplicon_of_position <- function(the_position){
#   if (nrow(subset(df_amplicon_positions,(the_position>=start)&(the_position<=stop)))>0){
#     return(subset(df_amplicon_positions,(the_position>=start)&(the_position<=stop))$amplicon_name)
#   }else{
#     return("")
#   }
# }
#Original codon and mutated codon
df_variants$ref_codon <- NA
df_variants$mut_codon <- NA
df_variants$pos_in_ORF <- NA
df_variants$pos_in_gene <- NA
df_variants$pos_in_protein <- NA
for (i in 1:nrow(df_variants)){
  df_variants$ref_codon[i] <-(get_ref_and_mutated_codon(the_position = df_variants$Position[i],ref_nucl = df_variants$Ref[i],new_nucl = df_variants$VarAllele[i]))$ref_codon
  df_variants$mut_codon[i] <-(get_ref_and_mutated_codon(the_position = df_variants$Position[i],ref_nucl = df_variants$Ref[i],new_nucl = df_variants$VarAllele[i]))$mutated_codon
  df_variants$old_aa[i] <- translate_seq(the_codon = df_variants$ref_codon[i] )
  df_variants$new_aa[i] <- translate_seq(the_codon = df_variants$mut_codon[i] )
  df_variants$pos_in_ORF[i] <- df_variants$Position[i] - v_start_orfs[df_variants$ORF[i]] + 1
  df_variants$pos_in_gene[i] <- df_variants$Position[i] - v_start_genes[df_variants$gene[i]] + 1
  df_variants$pos_in_protein[i] <- ceiling(df_variants$pos_in_gene[i]/3)
  
}
df_variants$mutation_name <- paste0(paste0(df_variants$Ref,df_variants$Position,df_variants$VarAllele,""),";",paste0(df_variants$old_aa,df_variants$pos_in_protein,df_variants$new_aa),";",df_variants$ORF,";",df_variants$gene)
#Define Nonsense and non-coding mutations
df_variants$is_nonsense <- (df_variants$new_aa=="Stop")
df_variants$is_UTR <- (is.na(df_variants$new_aa))

#import metadata 
df_Metadata_samples <- read.csv2(file = paste0(output_workspace,"Metadata_samples.csv"),sep = ";",header = TRUE,stringsAsFactors = FALSE)
colnames(df_Metadata_samples)[1] <- "Sample"
df_Metadata_samples$Ct <- as.numeric(df_Metadata_samples$Ct)
df_Metadata_samples$label_cardinal_point <- NA
df_Metadata_samples$label_search_sample <- NA
df_Metadata_samples$label_search_location <- NA
df_Metadata_samples$label_search_site_id_1 <- NA
df_Metadata_samples$label_search_site_id_2 <- NA
df_Metadata_samples$label_search_date_format_1 <- NA
df_Metadata_samples$label_search_date_format_2 <- NA
for (i in 1:nrow(df_Metadata_samples)){
  current_v_pos_sep <- gregexpr(pattern = "_",text = df_Metadata_samples$Sample[i],fixed = T)[[1]]
  df_Metadata_samples$label_search_sample[i] <- paste0(df_Metadata_samples$Sample[i],"_2")
  df_Metadata_samples$label_search_location[i] <- substr(df_Metadata_samples$Sample[i],1,current_v_pos_sep[1]-1)
  df_Metadata_samples$label_search_site_id_1[i] <- paste0("_",as.character(as.integer(substr(df_Metadata_samples$Sample[i],current_v_pos_sep[1]+1,current_v_pos_sep[2]-1))),"_")
  df_Metadata_samples$label_search_site_id_2[i] <- paste0(substr(df_Metadata_samples$Sample[i],current_v_pos_sep[1]+1,current_v_pos_sep[2]-1),"_")
  df_Metadata_samples$label_search_date_format_1[i] <- as.character(format(as.Date(df_Metadata_samples$Sampling_date[i],"%Y-%m-%d"),"%d-%m-%Y"))
  df_Metadata_samples$label_search_date_format_2[i] <- as.character(format(as.Date(df_Metadata_samples$Sampling_date[i],"%Y-%m-%d"),"%d-%m-%y"))
  if (grepl(pattern = "NORTH",x = toupper(df_Metadata_samples$Sampling_site[i]),fixed = T)){
    df_Metadata_samples$label_cardinal_point[i] <- "NORTH"
  }else if (grepl(pattern = "SOUTH",x = toupper(df_Metadata_samples$Sampling_site[i]),fixed = T)){
    df_Metadata_samples$label_cardinal_point[i] <- "SOUTH"
  }else if (grepl(pattern = "EAST",x = toupper(df_Metadata_samples$Sampling_site[i]),fixed = T)){
    df_Metadata_samples$label_cardinal_point[i] <- "EAST"
  }else if (grepl(pattern = "WEST",x = toupper(df_Metadata_samples$Sampling_site[i]),fixed = T)){
    df_Metadata_samples$label_cardinal_point[i] <- "WEST"
  }
}

#calculate the number of SNVs per samples 
v_nb_SNVs_per_sample <- as.vector(table(df_variants$Sample)[lst_samples_original])
names(v_nb_SNVs_per_sample) <- lst_samples_original
v_nb_SNVs_per_sample <- ifelse(test = is.na(v_nb_SNVs_per_sample),yes=0,no=v_nb_SNVs_per_sample)
#coverage plot 
#average depth per sample
df_raw_depth_report <- read.csv2(file = paste0(output_workspace,depth_report_filename),sep = "\t",header = FALSE,stringsAsFactors = FALSE)
v_avg_depth_samples <- colMeans(df_raw_depth_report[,3:ncol(df_raw_depth_report)],na.rm = TRUE)
df_cov_avg_depth_per_sample = data.frame(sample=lst_samples_original,nb_SNVs=v_nb_SNVs_per_sample,avg_depth=v_avg_depth_samples,Ct=NA,stringsAsFactors = FALSE)
for (i in 1:nrow(df_cov_avg_depth_per_sample)){
  if (sum(unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_sample[j]),x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T)&(grepl(pattern = df_Metadata_samples$label_search_sample[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))&(!is.na(df_Metadata_samples$Sampling_date[j])), FUN.VALUE = c(T))))==1){
    df_cov_avg_depth_per_sample$Ct[i] <- subset(df_Metadata_samples,unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_sample[j]),x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T)&(grepl(pattern = df_Metadata_samples$label_search_sample[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))&(!is.na(df_Metadata_samples$Sampling_date[j])), FUN.VALUE = c(T))))$Ct[1]
  }else{
    if (any(unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_location[j]),x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T)&((grepl(pattern = df_Metadata_samples$label_search_date_format_1[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))|(grepl(pattern = df_Metadata_samples$label_search_date_format_2[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T)))&(!is.na(df_Metadata_samples$Sampling_date[j]))&((grepl(pattern = df_Metadata_samples$label_search_site_id_1[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))|((grepl(pattern = df_Metadata_samples$label_cardinal_point[j],x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T))&(!is.na(df_Metadata_samples$label_cardinal_point[j])))|(grepl(pattern = df_Metadata_samples$label_search_site_id_2[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))), FUN.VALUE = c(T))))){
      if ((grepl(pattern = "GRB",x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T))|(grepl(pattern = "GRAB",x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T))){
        df_cov_avg_depth_per_sample$Ct[i] <- subset(df_Metadata_samples, unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_location[j]),x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T)&((grepl(pattern = df_Metadata_samples$label_search_date_format_1[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))|(grepl(pattern = df_Metadata_samples$label_search_date_format_2[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T)))&(!is.na(df_Metadata_samples$Sampling_date[j]))&((grepl(pattern = df_Metadata_samples$label_search_site_id_1[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))|((grepl(pattern = df_Metadata_samples$label_cardinal_point[j],x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T))&(!is.na(df_Metadata_samples$label_cardinal_point[j])))|(grepl(pattern = df_Metadata_samples$label_search_site_id_2[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T)))&((grepl(pattern = "GRB",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))|(grepl(pattern = "GRAB",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))), FUN.VALUE = c(T))))$Ct[1]
      }else{
        df_cov_avg_depth_per_sample$Ct[i] <- subset(df_Metadata_samples, unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_location[j]),x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T)&((grepl(pattern = df_Metadata_samples$label_search_date_format_1[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))|(grepl(pattern = df_Metadata_samples$label_search_date_format_2[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T)))&(!is.na(df_Metadata_samples$Sampling_date[j]))&((grepl(pattern = df_Metadata_samples$label_search_site_id_1[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T))|((grepl(pattern = df_Metadata_samples$label_cardinal_point[j],x = toupper(df_cov_avg_depth_per_sample$sample[i]),fixed=T))&(!is.na(df_Metadata_samples$label_cardinal_point[j])))|(grepl(pattern = df_Metadata_samples$label_search_site_id_2[j],x = df_cov_avg_depth_per_sample$sample[i],fixed=T)))&((grepl(pattern = "CPT",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))|(grepl(pattern = "COMPOSITE",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))), FUN.VALUE = c(T))))$Ct[1]
      }
    }else{
      df_cov_avg_depth_per_sample$Ct[i] <- NA
    }
  }
}
df_cov_avg_depth_per_sample$Ct <- as.numeric(df_cov_avg_depth_per_sample$Ct)
#df_cov_avg_depth_per_sample$log10_nb_SNVs <- log10(df_cov_avg_depth_per_sample$nb_SNVs)
df_cov_avg_depth_per_sample$log10_avg_depth <- log10(df_cov_avg_depth_per_sample$avg_depth)
ggplotRegression(fit = lmp(nb_SNVs~avg_depth,data = df_cov_avg_depth_per_sample,center = FALSE,Iter=9999),ggsave_path = output_workspace,the_filename = "nb_SNVs_vs_Coverage_Qc.png",xlabl = "Sample average coverage",ylabl = "Number of SNVs per sample")
ggplotRegression(fit = lmp(nb_SNVs~Ct,data = df_cov_avg_depth_per_sample,center = FALSE,Iter=9999),ggsave_path = output_workspace,the_filename = "nb_SNVs_vs_Ct_Qc.png",xlabl = "Ct",ylabl = "Number of SNVs per sample")
ggplotRegression(fit = lmp(log10_avg_depth~Ct,data = subset(df_cov_avg_depth_per_sample,is.finite(log10_avg_depth)),center = FALSE,Iter=9999),ggsave_path = output_workspace,the_filename = "log10_Coverage_vs_Ct_Qc.png",xlabl = "Ct",ylabl = "log10(Sample average coverage)")

df_depth <- data.frame(sample=lst_samples_original[1],position=df_raw_depth_report[,2],depth=df_raw_depth_report[,3],stringsAsFactors = FALSE)
for (i in 2:nb_samples_original){
  df_depth <- rbind(df_depth,data.frame(sample=lst_samples_original[i],position=df_raw_depth_report[,2],depth=df_raw_depth_report[,2+i],stringsAsFactors = FALSE))
}
df_depth$ORF <- unname(vapply(X = df_depth$position,FUN = find_ORF_of_mutation,FUN.VALUE = c("")))
df_depth <- unique(df_depth)

# #amplicon positions
# df_amplicon_positions <- read.csv2(file = paste0(output_workspace,"amplicon_positions.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
# 
# #mean coverage and breadth of coverage per Sample for each amplicon
# v_lst_unique_amplicons <- unique(df_amplicon_positions$amplicon_name)
# df_plot_mean_cov_per_sample_by_amplicon <- data.frame(Sample=rep(sort(unique(df_depth$sample)),length(v_lst_unique_amplicons)),amplicon=rep(v_lst_unique_amplicons,nb_samples_original),avg_cov=NA,breadth_of_cov=NA,stringsAsFactors = F)
# 
# lst_splits <- split(1:nrow(df_plot_mean_cov_per_sample_by_amplicon), ceiling(seq_along(1:nrow(df_plot_mean_cov_per_sample_by_amplicon))/(nrow(df_plot_mean_cov_per_sample_by_amplicon)/nb_cores)))
# the_f_parallel <- function(i_cl){
#   the_vec<- lst_splits[[i_cl]]
#   df_plot_mean_cov_per_sample_by_amplicon_current_subset <- df_plot_mean_cov_per_sample_by_amplicon[the_vec,]
#   count_iter <- 0
#   for (the_i in 1:nrow(df_plot_mean_cov_per_sample_by_amplicon_current_subset)){
#     start_current_amplicon <- subset(df_amplicon_positions,amplicon_name==df_plot_mean_cov_per_sample_by_amplicon_current_subset$amplicon[the_i])$start
#     stop_current_amplicon <- subset(df_amplicon_positions,amplicon_name==df_plot_mean_cov_per_sample_by_amplicon_current_subset$amplicon[the_i])$stop
#     length_current_amplicon <- stop_current_amplicon - start_current_amplicon + 1
#     df_plot_mean_cov_per_sample_by_amplicon_current_subset$avg_cov[the_i] <- mean(subset(df_depth,(sample==df_plot_mean_cov_per_sample_by_amplicon_current_subset$Sample[the_i])&(position>=start_current_amplicon)&(position<=stop_current_amplicon))$depth,na.rm=T)
#     df_plot_mean_cov_per_sample_by_amplicon_current_subset$breadth_of_cov[the_i] <- length(unique(subset(df_depth,(sample==df_plot_mean_cov_per_sample_by_amplicon_current_subset$Sample[the_i])&(position>=start_current_amplicon)&(position<=stop_current_amplicon)&(depth>=50)&(!is.na(depth)))$position))/length_current_amplicon
#     count_iter <- count_iter + 1
#     print(paste0("Core ",i_cl,": Step ",count_iter," done out of ",nrow(df_plot_mean_cov_per_sample_by_amplicon_current_subset),"!"))
#   }
#   return(df_plot_mean_cov_per_sample_by_amplicon_current_subset)
# }
# cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_df_plot_mean_cov_per_sample_by_amplicon.txt"))
# registerDoParallel(cl)
# df_plot_mean_cov_per_sample_by_amplicon <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel(i_cl)
# stopCluster(cl)
# df_plot_mean_cov_per_sample_by_amplicon <- unique(df_plot_mean_cov_per_sample_by_amplicon)
# saveRDS(df_plot_mean_cov_per_sample_by_amplicon,paste0(output_workspace,"df_plot_mean_cov_per_sample_by_amplicon.rds"))
# ggplot(data = df_plot_mean_cov_per_sample_by_amplicon,mapping=aes(x=factor(amplicon,levels=v_lst_unique_amplicons),y=log10(avg_cov))) + geom_violin(fill="red3") + geom_boxplot(width=0.05) + geom_jitter() + xlab("amplicon") + ylab("log10(Mean coverage per sample)") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=7,angle = 60,hjust=1))
# ggsave(filename = "amplicon_mean_coverage_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")
# ggplot(data = df_plot_mean_cov_per_sample_by_amplicon,mapping=aes(x=factor(amplicon,levels=v_lst_unique_amplicons),y=((breadth_of_cov)))) + geom_violin(fill="red3") + geom_jitter() + xlab("amplicon") + ylab("Breadth of coverage per sample") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=7,angle = 60,hjust=1))
# ggsave(filename = "amplicon_breadth_of_coverage_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")
# #heatmap mean cov per amplicon in each sample
# mtx_hmp_mean_cov_per_sample_for_each_amplicon <- reshape2::acast(data = df_plot_mean_cov_per_sample_by_amplicon, formula = Sample~amplicon, value.var="avg_cov")
# #heatmap breadth of cov per amplicon in each sample
# mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon <- reshape2::acast(df_plot_mean_cov_per_sample_by_amplicon, Sample~amplicon, value.var="breadth_of_cov")
# 
# #proportion of samples in which mean depth of coverage is inferior to 100x + proportion of samples in which breadth of coverage is inferior to 0.5
# df_amplicons_proportion_samples_with_low_coverage <- df_amplicon_positions
# df_amplicons_proportion_samples_with_low_coverage$ORF <- unname(vapply(X = df_amplicons_proportion_samples_with_low_coverage$start,FUN = find_ORF_of_mutation,FUN.VALUE = c("")))
# df_amplicons_proportion_samples_with_low_coverage$proportion_of_samples_with_mean_depth_lower_than_100x <- unname(vapply(X = df_amplicons_proportion_samples_with_low_coverage$amplicon_name,FUN = function(x) sum(mtx_hmp_mean_cov_per_sample_for_each_amplicon[,x]<100,na.rm=T)/nb_samples_original,FUN.VALUE = c(0.0)))
# df_amplicons_proportion_samples_with_low_coverage$proportion_of_samples_with_breadth_of_cov_lower_than_50pct <- unname(vapply(X = df_amplicons_proportion_samples_with_low_coverage$amplicon_name,FUN = function(x) sum(mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon[,x]<0.5,na.rm=T)/nb_samples_original,FUN.VALUE = c(0.0)))
# df_amplicons_proportion_samples_with_low_coverage$ORF <- unname(vapply(X = df_amplicons_proportion_samples_with_low_coverage$start,FUN = find_ORF_of_mutation,FUN.VALUE = c("")))
# write.table(x=df_amplicons_proportion_samples_with_low_coverage,file = paste0(output_workspace,"Table_amplicons_metadata_and_coverage_analysis_results.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

#per genomic region
#ggplot(data = subset(df_amplicons_proportion_samples_with_low_coverage,!is.na(ORF)),mapping=aes(x=factor(ORF,levels=v_orfs),y=proportion_of_samples_with_mean_depth_lower_than_100x)) + geom_violin(fill="red3") + geom_boxplot(width=0.05) + geom_jitter() + xlab("Genomic region") + ylab("Proportion of samples in which amplicon mean depth < 100") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=10,angle = 60,hjust=1))
#ggsave(filename = "Distribution_of_proportion_samples_in_which_amplicon_has_low_coverage_per_Genomic_region.png", path=output_workspace, width = 20, height = 15, units = "cm")

#PCA of samples by SNVs VAF
mtx_VAF_per_sample <- (as.matrix((reshape2::acast(df_variants, Sample~mutation_name, value.var="VarFreq"))))
mtx_VAF_per_sample[is.na(mtx_VAF_per_sample)] <- 0
mtx_pres_abs_snvs_in_samples <- ifelse(mtx_VAF_per_sample>0,yes=1,no=0)
write.table(x=mtx_VAF_per_sample,file = paste0(output_workspace,"mtx_VAF_per_sample.csv"),sep = ",",na = "NA",row.names = TRUE,col.names = TRUE)
location_samples <- as.vector(unname(vapply(X = rownames(mtx_VAF_per_sample),FUN = function(x) subset(df_variants,Sample==x)$location[1],FUN.VALUE = c(""))))
location_samples <- ifelse(test = location_samples=="NA",yes = NA,no=location_samples)
site_samples <- unname(vapply(X = rownames(mtx_VAF_per_sample),FUN = get_site_of_sample,FUN.VALUE = c("")))
site_samples <- ifelse(test = site_samples=="NA",yes = NA,no=site_samples)

#PCA VAF across Samples
# the_pca_WW_samples <- prcomp(mtx_VAF_per_sample, center = T,scale. = T)
# ggbiplot::ggbiplot(the_pca_WW_samples,ellipse=T,labels=rownames(mtx_VAF_per_sample),var.axes = F)
#tSNE
res_tsne <- Rtsne::Rtsne(mtx_VAF_per_sample, dims = 2, perplexity=30, verbose=TRUE, max_iter = 9999)
ggplot(mapping=aes(x=res_tsne$Y[,1],y=res_tsne$Y[,2],col=location_samples)) + geom_point(size=2) + xlab("tsne1") + ylab("tsne2") + labs(col="Location") + theme_bw()
ggsave(filename = "Dimensional_reduction_tSNE_WW_samples_by_location.png", path=output_workspace, width = 20, height = 15, units = "cm")
ggplot(mapping=aes(x=res_tsne$Y[!is.na(site_samples),1],y=res_tsne$Y[!is.na(site_samples),2],col=site_samples[!is.na(site_samples)])) + geom_point(size=2) + xlab("tsne1") + ylab("tsne2") + labs(col="Site") + theme_bw()
ggsave(filename = "Dimensional_reduction_tSNE_WW_samples_by_sites.png", path=output_workspace, width = 20, height = 15, units = "cm")
#PHATE
# library("phateR")
#tree.phate <- phate(mtx_VAF_per_sample, gamma=0, t=120, init=phate(mtx_VAF_per_sample)$data)
# plot embedding
#palette(rainbow(10))
#plot(tree.phate, col = tree.data$branches)
#UMAP
# library("uwot")
# umap_ww <- umap(mtx_VAF_per_sample,n_neighbors = 15,min_dist = 1, spread = 5)
# 
# umap_ww <- data.frame(
#   UMAP1 = umap_ww[, 1],
#   UMAP2 = umap_ww[, 2],
#   classification = location_samples
# )
# 
# ggplot(umap_ww, aes(x = UMAP1, y = UMAP2, col = classification)) +geom_point()

#PERMANOVA
mtx_sample_dist <- vegan::vegdist(mtx_VAF_per_sample, method='bray')
set.seed(1234)#reproducible results
permanova_sample_vaf_profile_vs_location <- adonis2(mtx_sample_dist~location_samples, permutations = 9999, method="bray", strata="PLOT")
permanova_sample_vaf_profile_vs_location
mtx_sample_dist_for_sites <- vegan::vegdist(mtx_VAF_per_sample[!is.na(site_samples),], method='bray')
site_samples_for_sites <- site_samples[!is.na(site_samples)]
permanova_sample_vaf_profile_vs_site <- adonis2(mtx_sample_dist_for_sites~site_samples_for_sites, permutations = 9999, method="bray", strata="PLOT")
permanova_sample_vaf_profile_vs_site

#Heatmap_mean_cov_per_sample_for_each_amplicon
# svg(paste0(output_workspace,"Heatmap_mean_cov_per_sample_for_each_amplicon.svg"),width=20/2.54,height=16/2.54)
# heatmap.2(x = log10(mtx_hmp_mean_cov_per_sample_for_each_amplicon+(1)), key.xlab = "log10(mean coverage+1)", labRow = "",cex=3, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(7.5, 5),cexCol = 0.1,cexRow = 1, ColSideColors =palette_orfs[vapply(X = df_amplicon_positions$start,FUN =find_ORF_of_mutation ,FUN.VALUE = c(""))],xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = F,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Heatmap_mean_cov_per_sample_for_each_amplicon.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = log10(mtx_hmp_mean_cov_per_sample_for_each_amplicon+(1)), key.xlab = "log10(mean coverage+1)", labRow = "",cex=3, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(7.5, 5),cexCol = 0.1,cexRow = 1, ColSideColors =palette_orfs[vapply(X = df_amplicon_positions$start,FUN =find_ORF_of_mutation ,FUN.VALUE = c(""))],xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = F,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Clustered_Heatmap_mean_cov_per_sample_for_each_amplicon.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = log10(mtx_hmp_mean_cov_per_sample_for_each_amplicon+(1)), key.xlab = "log10(mean coverage+1)", labRow = "",cex=3, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(7.5, 5),cexCol = 0.1,cexRow = 1, xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()

#Heatmap_breadth_of_cov_per_sample_for_each_amplicon
# svg(paste0(output_workspace,"Heatmap_breadth_of_cov_per_sample_for_each_amplicon.svg"),width=20/2.54,height=16/2.54)
# heatmap.2(x = mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon, key.xlab = "Breadth of coverage",labRow = "",cex=2, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(6, 5),cexCol = 0.1,cexRow = 1, ColSideColors =palette_orfs[vapply(X = df_amplicon_positions$start,FUN =find_ORF_of_mutation ,FUN.VALUE = c(""))],xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = F,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Heatmap_breadth_of_cov_per_sample_for_each_amplicon.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon, key.xlab = "Breadth of coverage",labRow = "",cex=2, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(6, 5),cexCol = 0.1,cexRow = 1, ColSideColors =palette_orfs[vapply(X = df_amplicon_positions$start,FUN =find_ORF_of_mutation ,FUN.VALUE = c(""))],xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = F,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Clustered_Heatmap_breadth_of_cov_per_sample_for_each_amplicon.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon, key.xlab = "Breadth of coverage",labRow = "",cex=2, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(6, 5),cexCol = 0.1,cexRow = 1, xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()

# #QIAseq replicates 
# lst_QIAseq_samples <-  read.csv2(file = paste0( "D:/Mirror/Covid19/Wastewater/QIAseq/Nanopore/","lst_samples.txt"),sep = ",",header = F,stringsAsFactors = FALSE)[,1]
# is_QIAseq_sample_replicate <- function(the_sample){
#   out <- F
#   for (current_QIAseq_sample in lst_QIAseq_samples){
#     if (grepl(pattern = toupper(current_QIAseq_sample),x = toupper(the_sample),fixed = T)){
#       out <- T
#       break()
#     }
#   }
#   return(out)
# }
# mtx_hmp_mean_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates <- mtx_hmp_mean_cov_per_sample_for_each_amplicon[unname(vapply(X = rownames(mtx_hmp_mean_cov_per_sample_for_each_amplicon),FUN = is_QIAseq_sample_replicate,FUN.VALUE = F)),]
# mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates <- mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon[unname(vapply(X = rownames(mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon),FUN = is_QIAseq_sample_replicate,FUN.VALUE = F)),]
# png(paste0(output_workspace,"Heatmap_mean_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = log10(mtx_hmp_mean_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates+(1)), key.xlab = "log10(mean coverage+1)", labRow = "",cex=3, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(7.5, 5),cexCol = 0.1,cexRow = 1, ColSideColors =palette_orfs[vapply(X = df_amplicon_positions$start,FUN =find_ORF_of_mutation ,FUN.VALUE = c(""))],xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = F,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Clustered_Heatmap_mean_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = log10(mtx_hmp_mean_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates+(1)), key.xlab = "log10(mean coverage+1)", labRow = "",cex=3, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(7.5, 5),cexCol = 0.1,cexRow = 1, xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Heatmap_breadth_of_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates, key.xlab = "Breadth of coverage",labRow = "",cex=2, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(6, 5),cexCol = 0.1,cexRow = 1, ColSideColors =palette_orfs[vapply(X = df_amplicon_positions$start,FUN =find_ORF_of_mutation ,FUN.VALUE = c(""))],xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = F,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# png(paste0(output_workspace,"Clustered_Heatmap_breadth_of_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates.png"),width=20/2.54,height=16/2.54,res = 1200,units = "in")
# heatmap.2(x = mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates, key.xlab = "Breadth of coverage",labRow = "",cex=2, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(6, 5),cexCol = 0.1,cexRow = 1, xlab = "Amplicon", ylab = "Sample",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36))
# dev.off()
# ggplot(data = subset(df_plot_nb_mutations_per_sample,(Sample%in%rownames(mtx_hmp_breadth_of_cov_per_sample_for_each_amplicon_QIAseq_Illumina_replicates))&(!grepl(pattern="Ctrl",x=Sample,fixed=T)))) + geom_col(aes(x=reorder(Sample,-x),y=x),fill="black") + xlab("Sample") + ylab("Number of mutations")  + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=4.7, angle = 60,hjust=1)) 
# ggsave(filename = "Nb_mutations_per_sample_QIAseq_Illumina_replicates.png", path=output_workspace, width = 20, height = 20, units = "cm")


# #mean coverage per Sample for each ORF
df_plot_mean_cov_per_sample_by_ORF <- data.frame(Sample=rep(sort(unique(df_depth$sample)),length(v_orfs[!v_orfs%in%c("ORF3b","ORF3c","ORF9c")])),ORF=rep(v_orfs[!v_orfs%in%c("ORF3b","ORF3c","ORF9c")],nb_samples_original),stringsAsFactors = F)
df_plot_mean_cov_per_sample_by_ORF$avg_cov <- unname(vapply(X = 1:nrow(df_plot_mean_cov_per_sample_by_ORF),FUN = function(i) mean(subset(df_depth,(sample==df_plot_mean_cov_per_sample_by_ORF$Sample[i])&(ORF==df_plot_mean_cov_per_sample_by_ORF$ORF[i]))$depth,na.rm=T),FUN.VALUE = c(0.0)))
#[subset(df_depth,(sample==df_plot_mean_cov_per_sample_by_ORF$Sample[i])&(ORF==df_plot_mean_cov_per_sample_by_ORF$ORF[i]))$depth>0]
ggplot(data = df_plot_mean_cov_per_sample_by_ORF,mapping=aes(x=factor(ORF,levels=v_orfs[!v_orfs%in%c("ORF3b","ORF3c","ORF9c")]),y=log10(avg_cov))) + geom_violin(fill="red3") + geom_jitter() + geom_boxplot(width=0.1)  + xlab("Genomic region") + ylab("log10(Mean coverage per sample)") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(angle = 60,hjust=1))
ggsave(filename = "ORF_mean_coverage_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")

df_plot_mean_cov_per_sample_by_ORF$breadth_of_cov <- unname(vapply(X = 1:nrow(df_plot_mean_cov_per_sample_by_ORF),FUN = function(i) length(unique(subset(df_depth,(sample==df_plot_mean_cov_per_sample_by_ORF$Sample[i])&(ORF==df_plot_mean_cov_per_sample_by_ORF$ORF[i])&(depth>=50)&(!is.na(depth)))$position))/(v_orfs_length[df_plot_mean_cov_per_sample_by_ORF$ORF[i]]) ,FUN.VALUE = c(0.0)))
ggplot(data = df_plot_mean_cov_per_sample_by_ORF,mapping=aes(x=factor(ORF,levels=v_orfs[!v_orfs%in%c("ORF3b","ORF3c","ORF9c")]),y=((breadth_of_cov)))) + geom_jitter() + geom_violin(fill="red3") + xlab("Genomic region") + ylab("Breadth of coverage per sample") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(angle = 60,hjust=1))
ggsave(filename = "ORF_breadth_of_coverage_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")

# ggplot(data = df_depth,aes(x = position,y=log10(depth+1),fill=as.factor(sample)), width=.01) + geom_bar(position = "identity",stat="identity",alpha=0.6)+theme_bw() + ylab("Site sequencing depth (log10 scale)") + xlab("Genomic position (bp)") + theme(legend.position = "none") + xlim(c(0,30000)) + ylim(0,log10(max(df_depth$depth,na.rm=TRUE))+0.1) + geom_hline(yintercept = log10(400),lty=2,color="black")
# ggsave(filename = "depth_report_across_samples_barplot_Qc.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)
# ggplot() + geom_line(mapping = aes(x = position,y=log10(depth),col=sample),data = df_depth,alpha=0.4)+theme_bw() + ylab("Site sequencing depth (log10 scale)") + xlab("Genomic position (bp)") + theme(legend.position = "none") + xlim(c(0,30000)) + ylim(0,log10(max(df_depth$depth,na.rm=TRUE))+0.1) + geom_hline(mapping = aes(y=log10(df_depth$depth)),yintercept = log10(400),lty=2,color="black")
# ggsave(filename = "depth_report_across_samples_lines_Qc.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)
# ggplot(mapping = aes(x=log10(depth),y=log10(depth),col=sample,fill=sample),data = df_depth) + geom_area() + theme_bw() + ylab("Count") + xlab("Site sequencing depth (log10 scale)") + theme(legend.position = "none") #+ geom_vline(mapping = aes(xintercept = log10(400)),lty=2,color="black")
# ggsave(filename = "depth_distr_across_samples_Qc.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)
# #ggplot(data = df_depth,mapping = aes(x=position,y=depth,fill=as.factor(sample))) + geom_area(alpha=0.6 , size=1, colour="black") + scale_fill_viridis(discrete = TRUE)  +  theme_ipsum() + ylab("Count") + xlab("Site sequencing depth (log10 scale)") + theme(legend.position = "right") #+ geom_vline(mapping = aes(xintercept = log10(400)),lty=2,color="black")
# 
# ggplot() + geom_col(aes(x=lst_samples, y=log10(colMeans(df_raw_depth_report[,3:ncol(df_raw_depth_report)]))),fill="black") + theme_bw() + theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=8,angle=60,hjust=1)) + ylab("log10(average depth)") + xlab("Sample") + scale_y_continuous(breaks = seq(0,4,0.5),limits=c(0,4)) + geom_hline(yintercept = log10(400),lty=2,col="red")
# ggsave(filename = "avg_log10_depth_across_samples_Qc.png", path=output_workspace, width = 40, height = 20, units = "cm")

df_raw_depth_report$ORF <- vapply(X = df_raw_depth_report[,2],FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_raw_depth_report_orf3bc <- subset(df_raw_depth_report,((df_raw_depth_report$V2>=v_start_orfs["ORF3a"]) & (df_raw_depth_report$V2<=v_end_orfs["ORF3a"])))
df_raw_depth_report_orf3bc$ORF <- ifelse(test = ((df_raw_depth_report_orf3bc$V2>=v_start_orfs["ORF3b"]) & (df_raw_depth_report_orf3bc$V2<=v_end_orfs["ORF3b"])),yes="ORF3b",no=df_raw_depth_report_orf3bc$ORF)
df_raw_depth_report_orf3bc$ORF <- ifelse(test = ((df_raw_depth_report_orf3bc$V2>=v_start_orfs["ORF3c"]) & (df_raw_depth_report_orf3bc$V2<=v_end_orfs["ORF3c"])),yes="ORF3c",no=df_raw_depth_report_orf3bc$ORF)
df_raw_depth_report_orf3bc <- subset(df_raw_depth_report_orf3bc,ORF!="ORF3a")
df_raw_depth_report <- rbind(df_raw_depth_report,df_raw_depth_report_orf3bc)
df_raw_depth_report_ORF9c <- subset(df_raw_depth_report,((df_raw_depth_report$V2>=v_start_orfs["N"]) & (df_raw_depth_report$V2<=v_end_orfs["N"])))
df_raw_depth_report_ORF9c$ORF <- ifelse(test = ((df_raw_depth_report_ORF9c$V2>=v_start_orfs["ORF9c"]) & (df_raw_depth_report_ORF9c$V2<=v_end_orfs["ORF9c"])),yes="ORF9c",no=df_raw_depth_report_ORF9c$ORF)
df_raw_depth_report_ORF9c <- subset(df_raw_depth_report_ORF9c,ORF!="N")
df_raw_depth_report <- rbind(df_raw_depth_report,df_raw_depth_report_ORF9c)
v_avg_cov_ORF <- vapply(X = v_orfs,FUN = function(x) return(mean((rowMeans(subset(df_raw_depth_report,ORF==x)[,3:(3+nb_samples-1)])),na.rm=TRUE)),FUN.VALUE = c(0.1))
v_orfs_length <- v_end_orfs - v_start_orfs + 1 
df_raw_depth_report$gene <- vapply(X = df_raw_depth_report[,2],FUN = function(x) return(find_gene_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_raw_depth_report_orf3bc <- subset(df_raw_depth_report,((df_raw_depth_report$V2>=v_start_genes["ORF3a"]) & (df_raw_depth_report$V2<=v_end_genes["ORF3a"])))
df_raw_depth_report_orf3bc$gene <- ifelse(test = ((df_raw_depth_report_orf3bc$V2>=v_start_genes["ORF3b"]) & (df_raw_depth_report_orf3bc$V2<=v_end_genes["ORF3b"])),yes="ORF3b",no=df_raw_depth_report_orf3bc$gene)
df_raw_depth_report_orf3bc$gene <- ifelse(test = ((df_raw_depth_report_orf3bc$V2>=v_start_genes["ORF3c"]) & (df_raw_depth_report_orf3bc$V2<=v_end_genes["ORF3c"])),yes="ORF3c",no=df_raw_depth_report_orf3bc$gene)
df_raw_depth_report_orf3bc <- subset(df_raw_depth_report_orf3bc,gene!="ORF3a")
df_raw_depth_report <- rbind(df_raw_depth_report,df_raw_depth_report_orf3bc)
df_raw_depth_report_gene9c <- subset(df_raw_depth_report,((df_raw_depth_report$V2>=v_start_genes["N"]) & (df_raw_depth_report$V2<=v_end_genes["N"])))
df_raw_depth_report_gene9c$gene <- ifelse(test = ((df_raw_depth_report_gene9c$V2>=v_start_genes["gene9c"]) & (df_raw_depth_report_gene9c$V2<=v_end_genes["gene9c"])),yes="gene9c",no=df_raw_depth_report_gene9c$gene)
df_raw_depth_report_gene9c <- subset(df_raw_depth_report_gene9c,gene!="N")
df_raw_depth_report <- rbind(df_raw_depth_report,df_raw_depth_report_gene9c)
v_avg_cov_gene <- vapply(X = v_genes_with_unique_product,FUN = function(x) return(mean((rowMeans(subset(df_raw_depth_report,gene==x)[,3:(3+nb_samples-1)])),na.rm=TRUE)),FUN.VALUE = c(0.1))
v_genes_length <- v_end_genes - v_start_genes + 1 

#Number of mutations per sample
v_nb_SNVs_per_sample_with_all_original_samples <- as.vector(table(df_variants$Sample)[lst_samples_original])
names(v_nb_SNVs_per_sample_with_all_original_samples) <- lst_samples_original
v_nb_SNVs_per_sample_with_all_original_samples <- ifelse(test = is.na(v_nb_SNVs_per_sample_with_all_original_samples),yes=0,no=v_nb_SNVs_per_sample_with_all_original_samples)
df_plot_nb_mutations_per_sample <- data.frame(Sample=names(v_nb_SNVs_per_sample_with_all_original_samples),x=v_nb_SNVs_per_sample_with_all_original_samples,stringsAsFactors = F)
Nb_mutations_per_sample_gg <- ggplot(data = subset(df_plot_nb_mutations_per_sample,!grepl(pattern="Ctrl",x=Sample,fixed=T))) + geom_col(aes(x=reorder(Sample,-x),y=x),fill="black") + xlab("Sample") + ylab("Number of mutations")  + theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.title = element_text(size=16),legend.text = element_text(size=12),axis.text.x = element_blank()) 
Nb_mutations_per_sample_gg
ggsave(filename = "Nb_mutations_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")
summary(df_plot_nb_mutations_per_sample$x);sd(df_plot_nb_mutations_per_sample$x);

#add collection date to df_variants 
get_date_of_sample <- function(the_sample_name){
  if (sum(unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_sample[j]),x = toupper(the_sample_name),fixed=T)&(grepl(pattern = df_Metadata_samples$label_search_sample[j],x = the_sample_name,fixed=T))&(!is.na(df_Metadata_samples$Sampling_date[j])), FUN.VALUE = c(T))))==1){
    return(subset(df_Metadata_samples,unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_sample[j]),x = toupper(the_sample_name),fixed=T)&(grepl(pattern = df_Metadata_samples$label_search_sample[j],x = the_sample_name,fixed=T))&(!is.na(df_Metadata_samples$Sampling_date[j])), FUN.VALUE = c(T))))$Sampling_date[1])
  }else{
    if (any(unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_location[j]),x = toupper(the_sample_name),fixed=T)&((grepl(pattern = df_Metadata_samples$label_search_date_format_1[j],x = the_sample_name,fixed=T))|(grepl(pattern = df_Metadata_samples$label_search_date_format_2[j],x = the_sample_name,fixed=T)))&(!is.na(df_Metadata_samples$Sampling_date[j]))&((grepl(pattern = df_Metadata_samples$label_search_site_id_1[j],x = the_sample_name,fixed=T))|((grepl(pattern = df_Metadata_samples$label_cardinal_point[j],x = toupper(the_sample_name),fixed=T))&(!is.na(df_Metadata_samples$label_cardinal_point[j])))|(grepl(pattern = df_Metadata_samples$label_search_site_id_2[j],x = the_sample_name,fixed=T))), FUN.VALUE = c(T))))){
      if ((grepl(pattern = "GRB",x = toupper(the_sample_name),fixed=T))|(grepl(pattern = "GRAB",x = toupper(the_sample_name),fixed=T))){
        return(subset(df_Metadata_samples, unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_location[j]),x = toupper(the_sample_name),fixed=T)&((grepl(pattern = df_Metadata_samples$label_search_date_format_1[j],x = the_sample_name,fixed=T))|(grepl(pattern = df_Metadata_samples$label_search_date_format_2[j],x = the_sample_name,fixed=T)))&(!is.na(df_Metadata_samples$Sampling_date[j]))&((grepl(pattern = df_Metadata_samples$label_search_site_id_1[j],x = the_sample_name,fixed=T))|((grepl(pattern = df_Metadata_samples$label_cardinal_point[j],x = toupper(the_sample_name),fixed=T))&(!is.na(df_Metadata_samples$label_cardinal_point[j])))|(grepl(pattern = df_Metadata_samples$label_search_site_id_2[j],x = the_sample_name,fixed=T)))&((grepl(pattern = "GRB",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))|(grepl(pattern = "GRAB",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))), FUN.VALUE = c(T))))$Sampling_date[1])
      }else{
        return(subset(df_Metadata_samples, unname(vapply(X = 1:nrow(df_Metadata_samples),FUN = function(j) grepl(pattern = toupper(df_Metadata_samples$label_search_location[j]),x = toupper(the_sample_name),fixed=T)&((grepl(pattern = df_Metadata_samples$label_search_date_format_1[j],x = the_sample_name,fixed=T))|(grepl(pattern = df_Metadata_samples$label_search_date_format_2[j],x = the_sample_name,fixed=T)))&(!is.na(df_Metadata_samples$Sampling_date[j]))&((grepl(pattern = df_Metadata_samples$label_search_site_id_1[j],x = the_sample_name,fixed=T))|((grepl(pattern = df_Metadata_samples$label_cardinal_point[j],x = toupper(the_sample_name),fixed=T))&(!is.na(df_Metadata_samples$label_cardinal_point[j])))|(grepl(pattern = df_Metadata_samples$label_search_site_id_2[j],x = the_sample_name,fixed=T)))&((grepl(pattern = "CPT",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))|(grepl(pattern = "COMPOSITE",x = toupper(df_Metadata_samples$Sample[j]),fixed=T))), FUN.VALUE = c(T))))$Sampling_date[1])
      }
    }else{
      return("NA")
    }
  }
}

df_variants$date <- unname(vapply(X = df_variants$Sample,FUN = get_date_of_sample,FUN.VALUE = c("")))
df_variants$date <- ifelse(test=df_variants$date=="NA",yes=NA,no=df_variants$date)

get_date_from_sample_name_with_date_en <- function(the_sample_name){
  pos_year <- gregexpr(pattern = "_202",text = the_sample_name,fixed = T)[[1]][1]
  pos_start_date <- pos_year + 1
  pos_end_date <- pos_year + 10
  return(substr(x = the_sample_name,start=pos_start_date,stop=pos_end_date))
}
get_date_from_sample_name_with_date_en2 <- function(the_sample_name){
  pos_year <- gregexpr(pattern = "-202",text = the_sample_name,fixed = T)[[1]][1]
  pos_start_date <- pos_year + 1
  pos_end_date <- pos_year + 10
  return(substr(x = the_sample_name,start=pos_start_date,stop=pos_end_date))
}
df_variants$date <- ifelse(test = (is.na(df_variants$date))&(unname(vapply(X = df_variants$Sample,FUN = function(x) grepl(pattern = "_202",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_variants$Sample,FUN=get_date_from_sample_name_with_date_en,FUN.VALUE = c(""))),no = df_variants$date)
df_variants$date <- gsub(pattern = "2021_",replacement = "2021-",x = df_variants$date,fixed = T)
df_variants$date <- ifelse(test = (is.na(df_variants$date))&(unname(vapply(X = df_variants$Sample,FUN = function(x) grepl(pattern = "-202",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_variants$Sample,FUN=get_date_from_sample_name_with_date_en2,FUN.VALUE = c(""))),no = df_variants$date)
# 
get_date_from_ON_20_samples_with_date <- function(the_sample_name){
  pos_year <- gregexpr(pattern = "-20_",text = the_sample_name,fixed = T)[[1]][1]
  pos_start_date <- pos_year - 5
  pos_end_date <- pos_year + 2
  str_date <- paste0(substr(x = the_sample_name,start=pos_start_date,stop=pos_end_date-2),"20",substr(x = the_sample_name,start=pos_end_date-1,stop=pos_end_date))
  return(format(as.Date(str_date,format="%d-%m-%Y"),"%Y-%m-%d"))
}
df_variants$date <- ifelse(test = (is.na(df_variants$date))&(unname(vapply(X = df_variants$Sample,FUN = function(x) grepl(pattern = "-20_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_variants$Sample,FUN=get_date_from_ON_20_samples_with_date,FUN.VALUE = c(""))),no = df_variants$date)

get_date_from_ON_21_samples_with_date <- function(the_sample_name){
  pos_year <- gregexpr(pattern = "-21_",text = the_sample_name,fixed = T)[[1]][1]
  pos_start_date <- pos_year - 5
  pos_end_date <- pos_year + 2
  str_date <- paste0(substr(x = the_sample_name,start=pos_start_date,stop=pos_end_date-2),"20",substr(x = the_sample_name,start=pos_end_date-1,stop=pos_end_date))
  return(format(as.Date(str_date,format="%d-%m-%Y"),"%Y-%m-%d"))
}
df_variants$date <- ifelse(test = (is.na(df_variants$date))&(unname(vapply(X = df_variants$Sample,FUN = function(x) grepl(pattern = "-21_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_variants$Sample,FUN=get_date_from_ON_21_samples_with_date,FUN.VALUE = c(""))),no = df_variants$date)

get_date_from_sample_name_with_date_fr <- function(the_sample_name){
  pos_year <- gregexpr(pattern = "-2020_",text = the_sample_name,fixed = T)[[1]][1]
  pos_start_date <- pos_year - 5
  pos_end_date <- pos_year + 4
  return(format(as.Date(substr(x = the_sample_name,start=pos_start_date,stop=pos_end_date),format="%d-%m-%Y"),"%Y-%m-%d"))
}
df_variants$date <- ifelse(test = (unname(vapply(X = df_variants$Sample,FUN = function(x) grepl(pattern = "-2020_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_variants$Sample,FUN=get_date_from_sample_name_with_date_fr,FUN.VALUE = c(""))),no = df_variants$date)

get_date_from_sample_name_with_date_fr2 <- function(the_sample_name){
  pos_year <- gregexpr(pattern = "-2021_",text = the_sample_name,fixed = T)[[1]][1]
  pos_start_date <- pos_year - 5
  pos_end_date <- pos_year + 4
  
  return(format(as.Date(substr(x = the_sample_name,start=pos_start_date,stop=pos_end_date),format="%d-%m-%Y"),"%Y-%m-%d"))
}
df_variants$date <- ifelse(test = (unname(vapply(X = df_variants$Sample,FUN = function(x) grepl(pattern = "-2021_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_variants$Sample,FUN=get_date_from_sample_name_with_date_fr2,FUN.VALUE = c(""))),no = df_variants$date)

v_samples_date <- rep(NA, length(lst_samples_original))
names(v_samples_date) <- lst_samples_original
v_samples_date <- unname(vapply(X = names(v_samples_date),FUN = get_date_of_sample,FUN.VALUE = c("")))
v_samples_date <- ifelse(test=v_samples_date=="NA",yes=NA,no=v_samples_date)
names(v_samples_date) <- lst_samples_original
v_samples_date <- ifelse(test = (is.na(v_samples_date))&(unname(vapply(X = names(v_samples_date),FUN = function(x) grepl(pattern = "_202",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = names(v_samples_date),FUN=get_date_from_sample_name_with_date_en,FUN.VALUE = c(""))),no = v_samples_date)
names(v_samples_date) <- lst_samples_original
v_samples_date <- ifelse(test = (is.na(v_samples_date))&(unname(vapply(X = names(v_samples_date),FUN = function(x) grepl(pattern = "-202",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = names(v_samples_date),FUN=get_date_from_sample_name_with_date_en2,FUN.VALUE = c(""))),no = v_samples_date)
names(v_samples_date) <- lst_samples_original
v_samples_date <- gsub(pattern = "2021_",replacement = "2021-",x = v_samples_date,fixed = T)
names(v_samples_date) <- lst_samples_original
v_samples_date <- ifelse(test = (is.na(v_samples_date))&(unname(vapply(X = lst_samples_original,FUN = function(x) grepl(pattern = "-20_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = lst_samples_original,FUN=get_date_from_ON_20_samples_with_date,FUN.VALUE = c(""))),no = v_samples_date)
v_samples_date <- ifelse(test = (is.na(v_samples_date))&(unname(vapply(X = lst_samples_original,FUN = function(x) grepl(pattern = "-21_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = lst_samples_original,FUN=get_date_from_ON_21_samples_with_date,FUN.VALUE = c(""))),no = v_samples_date)
names(v_samples_date) <- lst_samples_original
v_samples_date <- ifelse(test = (unname(vapply(X = lst_samples_original,FUN = function(x) grepl(pattern = "-2020_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = lst_samples_original,FUN=get_date_from_sample_name_with_date_fr,FUN.VALUE = c(""))),no = v_samples_date)
v_samples_date <- ifelse(test = (unname(vapply(X = lst_samples_original,FUN = function(x) grepl(pattern = "-2021_",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = lst_samples_original,FUN=get_date_from_sample_name_with_date_fr2,FUN.VALUE = c(""))),no = v_samples_date)
names(v_samples_date) <- lst_samples_original
write.table(x = names(v_samples_date)[is.na(v_samples_date)],file = paste0(output_workspace,"List_samples_with_missing_sampling_date_ALL_ILLUMINA_SAMPLES.csv"),sep = ",",row.names = F,col.names = F)
paste0("Number of samples with a MISSING sampling date = ",nb_samples_original - sum(!is.na(v_samples_date))," out of ",nb_samples_original,"!")

v_samples_location <- unname(vapply(X = lst_samples_original,FUN = get_location_from_sample_name_with_location,FUN.VALUE = c("")))
v_samples_location <- ifelse(test=v_samples_location=="",yes=NA,no=v_samples_location)
names(v_samples_location) <- lst_samples_original

#function to get the position of a mutation from the mutation name
get_position_from_mut_name <- function(the_mut_name){
  pos_init <- 2
  pos_end <- gregexpr(pattern = ";",text = the_mut_name,fixed = T)[[1]][1] - 2
  return(as.integer(substr(the_mut_name,pos_init,pos_end)))
}

df_plot_time_series_mutation_prevalence_per_location <- aggregate(x = subset(df_variants, table(df_variants$mutation_name)[mutation_name]>1)$Sample, by=list(mutation_name=subset(df_variants, table(df_variants$mutation_name)[mutation_name]>1)$mutation_name,Sampling_date=subset(df_variants, table(df_variants$mutation_name)[mutation_name]>1)$date,location=subset(df_variants, table(df_variants$mutation_name)[mutation_name]>1)$location),FUN=function(x) length(unique(x)))
df_plot_time_series_mutation_prevalence_per_location$nb_samples_in_which_SNV_is_detectable <- unname(vapply(X = 1:nrow(df_plot_time_series_mutation_prevalence_per_location),FUN = function(i) length(unique(intersect(subset(df_depth,(position==get_position_from_mut_name(the_mut_name = df_plot_time_series_mutation_prevalence_per_location$mutation_name[i]))&(depth>=50))$sample,names(table(names((v_samples_date[v_samples_date==df_plot_time_series_mutation_prevalence_per_location$Sampling_date[i]]))))))),FUN.VALUE = c(0)))
df_plot_time_series_mutation_prevalence_per_location$prevalence_proportion <- df_plot_time_series_mutation_prevalence_per_location$x/df_plot_time_series_mutation_prevalence_per_location$nb_samples_in_which_SNV_is_detectable
ggplot(data = df_plot_time_series_mutation_prevalence_per_location,mapping=aes(x=Sampling_date,y=x)) + geom_col(mapping=aes(fill=factor(mutation_name,levels=names(sort(table(df_variants$mutation_name)[table(df_variants$mutation_name)>1],decreasing = F))))) + xlab("Sampling date") + ylab(paste0("Prevalence (Number of samples in which \nthe mutation is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=8),legend.text = element_text(size=6),legend.position="bottom",axis.text.x = element_text(size=8, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="Mutations detected in more than 1 sample")
ggsave(filename = "LEGEND_Time_series_mutation_prevalence_per_location.png", path=output_workspace, width = 40, height = 40, units = "cm",dpi=1200)
ggplot(data = df_plot_time_series_mutation_prevalence_per_location,mapping=aes(x=Sampling_date,y=x)) + geom_col(mapping=aes(fill=factor(mutation_name,levels=names(sort(table(df_variants$mutation_name)[table(df_variants$mutation_name)>1],decreasing = F))))) + xlab("Sampling date") + ylab(paste0("Prevalence (Number of samples in which \nthe mutation is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=8),legend.text = element_text(size=6),legend.position="none",axis.text.x = element_text(size=8, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="Mutations detected in more than 1 sample")
ggsave(filename = "Time_series_mutation_prevalence_per_location.png", path=output_workspace, width = 40, height = 40, units = "cm",dpi=1200)
# ggplot(data = df_plot_time_series_mutation_prevalence_per_location,mapping=aes(x=Sampling_date,y=prevalence_proportion)) + geom_col(mapping=aes(fill=factor(mutation_name,levels=names(sort(table(df_variants$mutation_name)[table(df_variants$mutation_name)>1],decreasing = F)))),position = position_dodge(width = 0.9)) + xlab("Sampling date") + ylab(paste0("Normalized prevalence \n(Number of samples in which the mutation is detected/Number of samples in which the mutation is detectable)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=8),legend.text = element_text(size=6),axis.text.x = element_text(size=8, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="Mutations detected in more than 1 sample")
# ggsave(filename = "Time_series_mutation_normalized_prevalence_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm",dpi=1200)

df_plot_time_series_S_mutation_of_interest_prevalence_per_location <- aggregate(x = subset(df_variants,mutation_name%in%c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S"))$Sample, by=list(mutation_name=subset(df_variants,mutation_name%in%c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S"))$mutation_name,Sampling_date=subset(df_variants,mutation_name%in%c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S"))$date,location=subset(df_variants,mutation_name%in%c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S"))$location),FUN=function(x) length(unique(x)))
df_plot_time_series_S_mutation_of_interest_prevalence_per_location$nb_samples_in_which_SNV_is_detectable <- unname(vapply(X = 1:nrow(df_plot_time_series_S_mutation_of_interest_prevalence_per_location),FUN = function(i) length(unique(intersect(subset(df_depth,(position==get_position_from_mut_name(the_mut_name = df_plot_time_series_S_mutation_of_interest_prevalence_per_location$mutation_name[i]))&(depth>=50))$sample,names(table(names((v_samples_date[v_samples_date==df_plot_time_series_S_mutation_of_interest_prevalence_per_location$Sampling_date[i]]))))))),FUN.VALUE = c(0)))
df_plot_time_series_S_mutation_of_interest_prevalence_per_location$prevalence_proportion <- df_plot_time_series_S_mutation_of_interest_prevalence_per_location$x/df_plot_time_series_S_mutation_of_interest_prevalence_per_location$nb_samples_in_which_SNV_is_detectable
ggplot(data = subset(df_plot_time_series_S_mutation_of_interest_prevalence_per_location,mutation_name%in%c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S")),mapping=aes(x=Sampling_date,y=prevalence_proportion)) + geom_col(mapping=aes(fill=factor(mutation_name,levels=names(sort(table(df_variants$mutation_name)[table(df_variants$mutation_name)>0],decreasing = F)))),position = position_dodge(width = 0.9)) + xlab("Sampling date") + ylab(paste0("Normalized prevalence (Number of detected\noccurences / Number of detectable occurences ;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=8),legend.text = element_text(size=6),legend.position="right",axis.text.x = element_text(size=8, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="S protein mutations of interest") + scale_fill_manual(values=palette_mutations_of_interest)
ggsave(filename = "Time_series_S_mutations_of_interest_normalized_prevalence_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm",dpi=1200)
ggplot(data = subset(df_plot_time_series_S_mutation_of_interest_prevalence_per_location,mutation_name%in%c("A23063T;N501Y;S;S","T22917G;L452R;S;S","G23012A;E484K;S;S","A23403G;D614G;S;S","C23604A;P681H;S;S","G22992A;S477N;S;S","A22812C;K417T;S;S")),mapping=aes(x=Sampling_date,y=x)) + geom_col(mapping=aes(fill=factor(mutation_name,levels=names(sort(table(df_variants$mutation_name)[table(df_variants$mutation_name)>0],decreasing = F))))) + xlab("Sampling date") + ylab(paste0("Prevalence (Number of samples in which \nthe mutation is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=8),legend.text = element_text(size=6),legend.position="right",axis.text.x = element_text(size=8, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="S protein mutations of interest") + scale_fill_manual(values=palette_mutations_of_interest)
ggsave(filename = "Time_series_S_mutations_of_interest_prevalence_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm",dpi=1200)

#Presence of lineage marker-mutations (>=90% prevalent only in a certain lineage consensus sequences) and inference of lineages presence
#function that converts genomic mutation into protein mutation
convert_genomic_mut_to_prot_mut <- function(the_genomic_mut){
  ref_nucl <- substr(the_genomic_mut,1,1)
  the_position <- as.integer(substr(the_genomic_mut,2,nchar(the_genomic_mut)-1))
  new_nucl <- substr(the_genomic_mut,nchar(the_genomic_mut),nchar(the_genomic_mut))
  the_orf <- find_ORF_of_mutation(the_position)
  the_gene <- find_gene_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    return(NA)
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
      mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
    }else if (pos_in_codon==2){
      ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
    }else if (pos_in_codon==3){
      ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  if (nchar(ref_codon)!=3){
    stop("codon length should be 3!")
  }
  pos_in_prot <- unname(ceiling((the_position - v_start_genes[the_gene] + 1)/3))
  ref_aa <- translate_seq(the_codon = ref_codon)
  new_aa <- translate_seq(the_codon = mut_codon)
  return(as.character(paste0(the_gene,":",ref_aa,pos_in_prot,new_aa)))
}
#function that converts genomic mutation into ORF protein seq mutation
convert_genomic_mut_to_ORF_prot_mut <- function(the_genomic_mut){
  ref_nucl <- substr(the_genomic_mut,1,1)
  the_position <- as.integer(substr(the_genomic_mut,2,nchar(the_genomic_mut)-1))
  new_nucl <- substr(the_genomic_mut,nchar(the_genomic_mut),nchar(the_genomic_mut))
  the_orf <- find_ORF_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    return(NA)
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
      mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
    }else if (pos_in_codon==2){
      ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
    }else if (pos_in_codon==3){
      ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  if (nchar(ref_codon)!=3){
    stop("codon length should be 3!")
  }
  pos_in_ORF_prot_seq <- unname(ceiling((the_position - v_start_orfs[the_orf] + 1)/3))
  ref_aa <- translate_seq(the_codon = ref_codon)
  new_aa <- translate_seq(the_codon = mut_codon)
  return(as.character(paste0(the_orf,":",ref_aa,pos_in_ORF_prot_seq,new_aa)))
}

#mutation prevalence data
df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- readRDS(file = paste0(output_workspace,"Table_df_all_mutations_prevalence_in_lineages.rds"))
#if the prevalence dataframe is incomplete, add the appropriate information
if ((!"label_mut_ORF_effect"%in%colnames(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages))|(!"B.1.617.X"%in%df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages$lineage)){
  df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- subset(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages, (!lineage%in%v_lineages_of_interest)&(prevalence>=0.9))
  df_sig_muts_prevalence_data_to_add <- read.csv2(file = paste0(output_workspace,"df_signature_muts_prevalence_VOCs_VUIs.csv"),sep = ",",header = T,stringsAsFactors = FALSE)
  df_sig_muts_prevalence_data_to_add$prevalence <- as.numeric(df_sig_muts_prevalence_data_to_add$prevalence)/100
  df_sig_muts_prevalence_data_to_add$genomic_position <- unname(vapply(X = df_sig_muts_prevalence_data_to_add$mutation_name,FUN = function(x) as.integer(substr(x,1,nchar(x)-1)),FUN.VALUE = c(0)))
  df_sig_muts_prevalence_data_to_add$mutation_name <- paste0(unname(vapply(X = df_sig_muts_prevalence_data_to_add$genomic_position, FUN = function(x) substr(genome_refseq,x,x),FUN.VALUE = c(""))),df_sig_muts_prevalence_data_to_add$mutation_name)
  df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- rbind(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages[,c("lineage","mutation_name","prevalence")],df_sig_muts_prevalence_data_to_add[,c("lineage","mutation_name","prevalence")])
  df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages$label_mut_ORF_effect <- NA
  for (i in 1:nrow(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages)){
    df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages$label_mut_ORF_effect[i] <- convert_genomic_mut_to_ORF_prot_mut(the_genomic_mut = df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages$mutation_name[i])
  }
  jsonlite::write_json(x = df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages,path = paste0(output_workspace,"Database_all_mutations_prevalence_in_SC2_lineages_consensus_sequences_as_of_2021_07_08.json"))
  saveRDS(object = df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages,file = paste0(output_workspace,"Table_df_all_mutations_prevalence_in_lineages.rds"))
}
mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- (as.matrix((reshape2::acast(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages, lineage~mutation_name, value.var="prevalence"))))
mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages[is.na(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages)] <- 0
#remove sub-lineages (level3+)
mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages[(unname(vapply(X = rownames(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages),FUN = function(x) sum(gregexpr(pattern = ".",text = x,fixed = T)[[1]]>=0),FUN.VALUE = c(0)))<3)|
                                                                                                                         (!(unname(vapply(X = rownames(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages),FUN = function(the_lin) any(unname(vapply(X = paste0(v_lineages_of_interest,"."),FUN = function(x) grepl(pattern = x,x = the_lin,fixed = T),FUN.VALUE = c(F)))),FUN.VALUE = c(F))))),]
rownames(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages)[which(rownames(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages)=="B.1.617")] <- "B.1.617.X"

saveRDS(object = mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages,file = paste0(output_workspace,"mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages.rds"))
mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- ifelse(test = mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages>=0.9,yes = 1,no=0)

v_lst_marker_mutations <- colnames(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)[colSums(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)==1]
v_lineage_marker_mutations <- NULL
for (current_marker_mut in v_lst_marker_mutations){
  v_lineage_marker_mutations <- c(v_lineage_marker_mutations,rownames(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)[which(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[,current_marker_mut]==1)])
  names(v_lineage_marker_mutations)[length(v_lineage_marker_mutations)] <- current_marker_mut
}
df_variants$pos_in_ORF_protein_seq <- vapply(X = 1:nrow(df_variants),FUN = function(i) ceiling((df_variants$Position[i] - v_start_orfs[df_variants$ORF[i]] + 1)/3),FUN.VALUE = c(0))
df_variants$label_mut_in_marker_fmt <- paste0(df_variants$Ref,df_variants$Position,df_variants$VarAllele) #paste0(df_variants$ORF,":",df_variants$old_aa,df_variants$pos_in_ORF_protein_seq,df_variants$new_aa)
v_label_marker_mut_to_mut_name <- unique(df_variants[,c("label_mut_in_marker_fmt","mutation_name")])$mutation_name
names(v_label_marker_mut_to_mut_name) <- unique(df_variants[,c("label_mut_in_marker_fmt","mutation_name")])$label_mut_in_marker_fmt
df_detected_marker_mutations_in_ww_samples <- subset(df_variants,label_mut_in_marker_fmt%in%names(v_lineage_marker_mutations))
df_detected_marker_mutations_in_ww_samples$PANGO_lineage <- v_lineage_marker_mutations[df_detected_marker_mutations_in_ww_samples$label_mut_in_marker_fmt]
df_detected_marker_mutations_in_ww_samples$date <- ifelse(test = (is.na(df_detected_marker_mutations_in_ww_samples$date))&(unname(vapply(X = df_detected_marker_mutations_in_ww_samples$Sample,FUN = function(x) grepl(pattern = "_202",x = x,fixed = T),FUN.VALUE = c(T)))),yes = unname(vapply(X = df_detected_marker_mutations_in_ww_samples$Sample,FUN = get_date_from_sample_name_with_date_en,FUN.VALUE = c(""))),no = df_detected_marker_mutations_in_ww_samples$date)
#function to get the position of a mutation from the mutation name
get_position_from_label_mut <- function(the_label_mut){
  pos_init <- 2
  pos_end <- nchar(the_label_mut) - 1
  return(as.integer(substr(the_label_mut,pos_init,pos_end)))
}
#function to get max confidence score of lineage inference
get_lst_lineage_mutations_and_conf_score <- function(the_lineage, the_sample){
  v_current_lineage_signature_muts <- colnames(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)[mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[the_lineage,]==1] 
  v_label_lineage_muts <- v_label_marker_mut_to_mut_name[v_current_lineage_signature_muts[v_current_lineage_signature_muts%in%unique(subset(df_variants,(Sample==the_sample))$label_mut_in_marker_fmt)]]
  the_stringent_confidence_score <- sum(v_current_lineage_signature_muts%in%unique(subset(df_variants,(Sample==the_sample))$label_mut_in_marker_fmt))/length(v_current_lineage_signature_muts)
  nb_signature_mut <- length(v_label_lineage_muts)
  v_genomic_positions_signature_muts <- unname(vapply(X = v_current_lineage_signature_muts[v_current_lineage_signature_muts%in%unique(subset(df_variants,(Sample==the_sample))$label_mut_in_marker_fmt)],FUN = get_position_from_label_mut,FUN.VALUE = c(0)))
  v_site_depth_genomic_positions_signature_muts <- unname(vapply(X = v_genomic_positions_signature_muts,FUN = function(x) subset(df_depth,(position==x)&(sample==the_sample))$depth[1],FUN.VALUE = c(0)))
  the_coverage_adjusted_confidence_score <- sum(v_current_lineage_signature_muts%in%unique(subset(df_variants,(Sample==the_sample))$label_mut_in_marker_fmt))/sum(v_site_depth_genomic_positions_signature_muts>=50)
  return(list(label_marker_muts = paste0(v_label_lineage_muts,collapse = "/"), stringent_confidence_score = the_stringent_confidence_score,nb_signature_mutations = nb_signature_mut, coverage_adjusted_confidence_score = the_coverage_adjusted_confidence_score))
}
df_detected_marker_mutations_in_ww_samples$str_signature_mutations <- NA
df_detected_marker_mutations_in_ww_samples$nb_signature_mutations <- NA
df_detected_marker_mutations_in_ww_samples$stringent_confidence_score <- NA
df_detected_marker_mutations_in_ww_samples$coverage_adjusted_confidence_score <- NA
for (i in 1:nrow(df_detected_marker_mutations_in_ww_samples)){
  list_lineage_mutations_and_conf_score <- get_lst_lineage_mutations_and_conf_score(the_lineage = df_detected_marker_mutations_in_ww_samples$PANGO_lineage[i],the_sample = df_detected_marker_mutations_in_ww_samples$Sample[i])
  df_detected_marker_mutations_in_ww_samples$str_signature_mutations[i] <- list_lineage_mutations_and_conf_score$label_marker_muts
  df_detected_marker_mutations_in_ww_samples$nb_signature_mutations[i] <- list_lineage_mutations_and_conf_score$nb_signature_mutations
  df_detected_marker_mutations_in_ww_samples$stringent_confidence_score[i] <- list_lineage_mutations_and_conf_score$stringent_confidence_score
  df_detected_marker_mutations_in_ww_samples$coverage_adjusted_confidence_score[i] <- list_lineage_mutations_and_conf_score$coverage_adjusted_confidence_score
}
#Remove Ottawa samples 
df_detected_marker_mutations_in_ww_samples <- subset(df_detected_marker_mutations_in_ww_samples,(!is.na(location))&(location!="OTTAWA"))

#Number of signature mutations per sample
mtx_is_signature_mutation_lineage <- mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[,colSums(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)>0]
df_plot_nb_sig_mutations_per_sample <- data.frame(Sample=lst_samples_original,x=unname(vapply(X = lst_samples_original,FUN = function(x) sum(subset(df_variants,Sample==x)$label_mut_in_marker_fmt%in%colnames(mtx_is_signature_mutation_lineage)),FUN.VALUE = c(0))),stringsAsFactors = F)
Nb_sig_mutations_per_sample_gg <- ggplot(data = subset(df_plot_nb_sig_mutations_per_sample,!grepl(pattern="Ctrl",x=Sample,fixed=T))) + geom_col(aes(x=reorder(Sample,-x),y=x),fill="black") + xlab("Sample") + ylab("Number of signature mutations")  + theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.title = element_text(size=14),legend.text = element_text(size=12),axis.text.x = element_blank()) 
Nb_sig_mutations_per_sample_gg
ggsave(filename = "Nb_sig_mutations_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")
summary(df_plot_nb_sig_mutations_per_sample$x);sd(df_plot_nb_sig_mutations_per_sample$x);

#Number of marker mutations per sample
df_plot_nb_marker_mutations_per_sample <- data.frame(Sample=lst_samples_original,x=unname(vapply(X = lst_samples_original,FUN = function(x) sum(subset(df_variants,Sample==x)$label_mut_in_marker_fmt%in%names(v_lineage_marker_mutations)),FUN.VALUE = c(0))),stringsAsFactors = F)
Nb_marker_mutations_per_sample_gg <- ggplot(data = subset(df_plot_nb_marker_mutations_per_sample,!grepl(pattern="Ctrl",x=Sample,fixed=T))) + geom_col(aes(x=reorder(Sample,-x),y=x),fill="black") + xlab("Sample") + ylab("Number of marker mutations")  + theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.title = element_text(size=14),legend.text = element_text(size=12),axis.text.x = element_blank())
Nb_marker_mutations_per_sample_gg
ggsave(filename = "Nb_marker_mutations_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")
summary(df_plot_nb_marker_mutations_per_sample$x);sd(df_plot_nb_marker_mutations_per_sample$x);

#Number of lineages per sample
df_plot_nb_lineages_per_sample <- data.frame(Sample=lst_samples_original,x=unname(vapply(X = lst_samples_original,FUN = function(x) length(unique(subset(df_detected_marker_mutations_in_ww_samples,Sample==x)$PANGO_lineage)),FUN.VALUE = c(0))),stringsAsFactors = F)
Nb_lineages_per_sample_gg <- ggplot(data = subset(df_plot_nb_lineages_per_sample,!grepl(pattern="Ctrl",x=Sample,fixed=T))) + geom_col(aes(x=reorder(Sample,-x),y=x),fill="black") + xlab("Sample") + ylab("Number of lineages")  + theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.title = element_text(size=14),legend.text = element_text(size=12),axis.text.x = element_blank()) 
Nb_lineages_per_sample_gg
ggsave(filename = "Nb_lineages_per_sample.png", path=output_workspace, width = 20, height = 15, units = "cm")
summary(df_plot_nb_lineages_per_sample$x);sd(df_plot_nb_lineages_per_sample$x);

#Sample metrics distributions
png(filename = paste0(output_workspace,"Sample_metrics_distributions.png"),width = 30,height=20, units = "cm",res = 1200)
grid.arrange(Nb_mutations_per_sample_gg,Nb_sig_mutations_per_sample_gg,Nb_marker_mutations_per_sample_gg,Nb_lineages_per_sample_gg,nrow=2,ncol=2)
dev.off()

#function to get the vector of marker mutations of a specific lineage
get_marker_mutations_of_PANGO_lineage <- function(the_pango_lin){
  return(names(v_lineage_marker_mutations[unname(v_lineage_marker_mutations)==the_pango_lin]))
}

df_plot_time_series_lineage_prevalence_per_location <- aggregate(x = df_detected_marker_mutations_in_ww_samples$Sample, by=list(PANGO_lineage=df_detected_marker_mutations_in_ww_samples$PANGO_lineage,Sampling_date=df_detected_marker_mutations_in_ww_samples$date,location=df_detected_marker_mutations_in_ww_samples$location),FUN=function(x) length(unique(x)))
df_plot_time_series_lineage_prevalence_per_location$nb_samples_in_which_lineage_is_detectable <- unname(vapply(X = 1:nrow(df_plot_time_series_lineage_prevalence_per_location),FUN = function(i) length(unique(intersect(subset(df_depth,(position%in%unname(vapply(X = get_marker_mutations_of_PANGO_lineage(the_pango_lin = df_plot_time_series_lineage_prevalence_per_location$PANGO_lineage[i]),FUN = get_position_from_label_mut,FUN.VALUE = c(0))))&(depth>=50))$sample,names(table(names((v_samples_date[v_samples_date==df_plot_time_series_lineage_prevalence_per_location$Sampling_date[i]]))))))),FUN.VALUE = c(0)))
df_plot_time_series_lineage_prevalence_per_location$prevalence_proportion <- df_plot_time_series_lineage_prevalence_per_location$x/df_plot_time_series_lineage_prevalence_per_location$nb_samples_in_which_lineage_is_detectable
ggplot(data = df_plot_time_series_lineage_prevalence_per_location,mapping=aes(x=Sampling_date,y=x)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Prevalence (Number of samples in which \nthe lineage is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_y_continuous(limits=c(0,max(df_plot_time_series_lineage_prevalence_per_location$x,na.rm=T)),breaks=seq(0,max(df_plot_time_series_lineage_prevalence_per_location$x,na.rm=T)+2,2))
ggsave(filename = "Time_series_PANGO_lineage_prevalence_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm")
ggplot(data = df_plot_time_series_lineage_prevalence_per_location,mapping=aes(x=Sampling_date,y=prevalence_proportion)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Normalized prevalence")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1.01,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages")
ggsave(filename = "Time_series_PANGO_lineage_normalized_prevalence_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm")

df_plot_confidence_scores_PANGO_lineage_per_location <- unique(df_detected_marker_mutations_in_ww_samples[,c("Sample","location","date","PANGO_lineage","nb_signature_mutations","stringent_confidence_score","coverage_adjusted_confidence_score")])
df_plot_confidence_scores_PANGO_lineage_per_location <- subset(df_plot_confidence_scores_PANGO_lineage_per_location,!is.na(date))
ggplot(data = df_plot_confidence_scores_PANGO_lineage_per_location,mapping=aes(x=date,y=stringent_confidence_score)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Confidence score\n (proportion of lineage signature mutations detected)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1.01,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") #+ scale_fill_manual(values=c("B.1.1.306"="skyblue","B.1.1.7"="red2","B.1.147"="blue3","B.1.160"="green3","B.1.36.17"="tomato", "B.40"="tan2","B.1.22"="orange","B.1.2"="brown"))
ggsave(filename = "Time_series_PANGO_lineage_stringent_confidence_score_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm")
ggplot(data = df_plot_confidence_scores_PANGO_lineage_per_location,mapping=aes(x=date,y=coverage_adjusted_confidence_score)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Confidence score\n (proportion of detectable signature mutations that are detected)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1.01,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") #+ scale_fill_manual(values=c("B.1.1.306"="skyblue","B.1.1.7"="red2","B.1.147"="blue3","B.1.160"="green3","B.1.36.17"="tomato", "B.40"="tan2","B.1.22"="orange","B.1.2"="brown"))
ggsave(filename = "Time_series_PANGO_lineage_coverage_adj_confidence_score_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm")
ggplot(data = df_plot_confidence_scores_PANGO_lineage_per_location,mapping=aes(x=date,y=nb_signature_mutations)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Number of detected signature mutations")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,max(df_plot_confidence_scores_PANGO_lineage_per_location$nb_signature_mutations)+3),breaks = seq(0,max(df_plot_confidence_scores_PANGO_lineage_per_location$nb_signature_mutations)+3,3)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") #+ scale_fill_manual(values=c("B.1.1.306"="skyblue","B.1.1.7"="red2","B.1.147"="blue3","B.1.160"="green3","B.1.36.17"="tomato", "B.40"="tan2","B.1.22"="orange","B.1.2"="brown"))
ggsave(filename = "Time_series_PANGO_lineage_nb_detected_signature_mutations_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm")

df_plot_prevalence_PANGO_lineage_per_location <- aggregate(x = df_detected_marker_mutations_in_ww_samples$Sample, by=list(PANGO_lineage=df_detected_marker_mutations_in_ww_samples$PANGO_lineage,location=df_detected_marker_mutations_in_ww_samples$location),FUN=function(x) length(unique(x)))
ggplot(data = df_plot_prevalence_PANGO_lineage_per_location,mapping=aes(x=location,y=x)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) + xlab("Location") + ylab(paste0("Prevalence (Number of samples in which \nthe lineage is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1)) + labs(fill="PANGO lineages") + scale_y_continuous(limits=c(0,max(df_plot_prevalence_PANGO_lineage_per_location$x,na.rm=T)),breaks=seq(0,max(df_plot_prevalence_PANGO_lineage_per_location$x,na.rm=T)+2,2)) #+ scale_fill_manual(values=c("B.1.1.306"="skyblue","B.1.1.7"="red2","B.1.147"="blue3","B.1.160"="green3","B.1.36.17"="tomato", "B.40"="tan2","B.1.22"="orange","B.1.2"="brown"))
ggsave(filename = "prevalence_PANGO_lineage_per_location.png", path=output_workspace, width = 30, height = 20, units = "cm")

df_plot_time_series_lineages_of_interest_prevalence_per_location <- aggregate(x = subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$Sample, by=list(PANGO_lineage=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage,Sampling_date=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$date,location=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$location),FUN=function(x) length(unique(x)))
df_plot_time_series_lineages_of_interest_prevalence_per_location$nb_samples_in_which_lineage_is_detectable <- unname(vapply(X = 1:nrow(df_plot_time_series_lineages_of_interest_prevalence_per_location),FUN = function(i) length(unique(intersect(subset(df_depth,(position%in%unname(vapply(X = get_marker_mutations_of_PANGO_lineage(the_pango_lin = df_plot_time_series_lineages_of_interest_prevalence_per_location$PANGO_lineage[i]),FUN = get_position_from_label_mut,FUN.VALUE = c(0))))&(depth>=50))$sample,names(table(names((v_samples_date[v_samples_date==df_plot_time_series_lineages_of_interest_prevalence_per_location$Sampling_date[i]]))))))),FUN.VALUE = c(0)))
df_plot_time_series_lineages_of_interest_prevalence_per_location$prevalence_proportion <- df_plot_time_series_lineages_of_interest_prevalence_per_location$x/df_plot_time_series_lineages_of_interest_prevalence_per_location$nb_samples_in_which_lineage_is_detectable
df_plot_time_series_lineages_of_interest_prevalence_per_location$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_plot_time_series_lineages_of_interest_prevalence_per_location$PANGO_lineage]
ggplot(data = df_plot_time_series_lineages_of_interest_prevalence_per_location,mapping=aes(x=Sampling_date,y=x)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Prevalence (Number of samples in which \nthe lineage is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_y_continuous(limits=c(0,max(df_plot_time_series_lineages_of_interest_prevalence_per_location$x,na.rm=T)),breaks=seq(0,max(df_plot_time_series_lineages_of_interest_prevalence_per_location$x,na.rm=T)+2,2)) + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_prevalence_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm")
ggplot(data = df_plot_time_series_lineages_of_interest_prevalence_per_location,mapping=aes(x=Sampling_date,y=prevalence_proportion)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Normalized prevalence")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1.01,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_normalized_prevalence_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm")

df_plot_confidence_scores_PANGO_lineages_of_interest_per_location <- unique(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)[,c("Sample","location","date","PANGO_lineage","nb_signature_mutations","stringent_confidence_score","coverage_adjusted_confidence_score")])
df_plot_confidence_scores_PANGO_lineages_of_interest_per_location <- subset(df_plot_confidence_scores_PANGO_lineages_of_interest_per_location,!is.na(date))
df_plot_confidence_scores_PANGO_lineages_of_interest_per_location$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_plot_confidence_scores_PANGO_lineages_of_interest_per_location$PANGO_lineage]
ggplot(data = df_plot_confidence_scores_PANGO_lineages_of_interest_per_location,mapping=aes(x=date,y=stringent_confidence_score)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Confidence score\n (proportion of lineage signature mutations detected)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1.01,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_stringent_confidence_score_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm")
ggplot(data = df_plot_confidence_scores_PANGO_lineages_of_interest_per_location,mapping=aes(x=date,y=coverage_adjusted_confidence_score)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Confidence score\n (proportion of detectable signature mutations that are detected)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1.01,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_coverage_adj_confidence_score_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm")
ggplot(data = df_plot_confidence_scores_PANGO_lineages_of_interest_per_location,mapping=aes(x=date,y=nb_signature_mutations)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Number of detected signature mutations")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,max(df_plot_confidence_scores_PANGO_lineages_of_interest_per_location$nb_signature_mutations)+3),breaks = seq(0,max(df_plot_confidence_scores_PANGO_lineages_of_interest_per_location$nb_signature_mutations)+3,3)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_nb_detected_signature_mutations_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm")

df_plot_prevalence_PANGO_lineage_of_interest_per_location <- aggregate(x = subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$Sample, by=list(PANGO_lineage=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage,location=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$location),FUN=function(x) length(unique(x)))
df_plot_prevalence_PANGO_lineage_of_interest_per_location$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_plot_prevalence_PANGO_lineage_of_interest_per_location$PANGO_lineage]
ggplot(data = df_plot_prevalence_PANGO_lineage_of_interest_per_location,mapping=aes(x=location,y=x)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) + xlab("Location") + ylab(paste0("Prevalence (Number of samples in which \nthe lineage is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1)) + labs(fill="PANGO lineages") + scale_y_continuous(limits=c(0,max(df_plot_prevalence_PANGO_lineage_of_interest_per_location$x,na.rm=T)),breaks=seq(0,max(df_plot_prevalence_PANGO_lineage_of_interest_per_location$x,na.rm=T)+2,2)) + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "prevalence_PANGO_lineage_of_interest_per_location.png", path=output_workspace, width = 40, height = 20, units = "cm")

#lineage presence inference and confidence score per sample site
df_detected_marker_mutations_in_ww_samples$site <- unname(vapply(X = df_detected_marker_mutations_in_ww_samples$Sample,FUN = get_site_of_sample,FUN.VALUE = c("")))
df_detected_marker_mutations_in_ww_samples$site <- ifelse(test = df_detected_marker_mutations_in_ww_samples$site=="NA",yes = NA,no=df_detected_marker_mutations_in_ww_samples$site)
df_plot_confidence_scores_PANGO_lineage_per_site <- unique(df_detected_marker_mutations_in_ww_samples[,c("Sample","site","date","PANGO_lineage","nb_signature_mutations","stringent_confidence_score","coverage_adjusted_confidence_score")])
df_plot_confidence_scores_PANGO_lineage_per_site <- subset(df_plot_confidence_scores_PANGO_lineage_per_site,!is.na(date))
ggplot(data = df_plot_confidence_scores_PANGO_lineage_per_site,mapping=aes(x=date,y=nb_signature_mutations)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Number of detected signature mutations")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,max(df_plot_confidence_scores_PANGO_lineage_per_site$nb_signature_mutations)+3),breaks = seq(0,max(df_plot_confidence_scores_PANGO_lineage_per_site$nb_signature_mutations)+3,3)) + facet_wrap(~site, ncol=3) + labs(fill="PANGO lineages") #+ scale_fill_manual(values=c("B.1.1.306"="skyblue","B.1.1.7"="red2","B.1.147"="blue3","B.1.160"="green3","B.1.36.17"="tomato", "B.40"="tan2","B.1.22"="orange","B.1.2"="brown"))
ggsave(filename = "Time_series_PANGO_lineage_nb_detected_signature_mutations_per_site.png", path=output_workspace, width = 40, height = 30, units = "cm")
df_plot_confidence_scores_PANGO_lineages_of_interest_per_site <- unique(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)[,c("Sample","site","date","PANGO_lineage","nb_signature_mutations","stringent_confidence_score","coverage_adjusted_confidence_score")])
df_plot_confidence_scores_PANGO_lineages_of_interest_per_site <- subset(df_plot_confidence_scores_PANGO_lineages_of_interest_per_site,!is.na(date))
df_plot_confidence_scores_PANGO_lineages_of_interest_per_site$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_plot_confidence_scores_PANGO_lineages_of_interest_per_site$PANGO_lineage]
ggplot(data = df_plot_confidence_scores_PANGO_lineages_of_interest_per_site,mapping=aes(x=date,y=nb_signature_mutations)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) +
  xlab("Sampling date") + ylab(paste0("Number of detected signature mutations")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,max(df_plot_confidence_scores_PANGO_lineages_of_interest_per_site$nb_signature_mutations)+3),breaks = seq(0,max(df_plot_confidence_scores_PANGO_lineages_of_interest_per_site$nb_signature_mutations)+3,3)) + facet_wrap(~site, ncol=3) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_nb_detected_signature_mutations_per_site.png", path=output_workspace, width = 40, height = 30, units = "cm")
df_plot_prevalence_PANGO_lineage_of_interest_per_site <- aggregate(x = subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$Sample, by=list(PANGO_lineage=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage,site=subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$site),FUN=function(x) length(unique(x)))
df_plot_prevalence_PANGO_lineage_of_interest_per_site$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_plot_prevalence_PANGO_lineage_of_interest_per_site$PANGO_lineage]
ggplot(data = df_plot_prevalence_PANGO_lineage_of_interest_per_site,mapping=aes(x=site,y=x)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = position_dodge(width = 0.9)) + xlab("site") + ylab(paste0("Prevalence (Number of samples in which \nthe lineage is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.text.x = element_text(angle = 60,hjust = 1,size=12),axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8)) + labs(fill="PANGO lineages") + scale_y_continuous(limits=c(0,max(df_plot_prevalence_PANGO_lineage_of_interest_per_site$x,na.rm=T)),breaks=seq(0,max(df_plot_prevalence_PANGO_lineage_of_interest_per_site$x,na.rm=T)+2,2)) + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "prevalence_PANGO_lineage_of_interest_per_site.png", path=output_workspace, width = 40, height = 20, units = "cm")
df_plot_prevalence_PANGO_lineage_per_site <- aggregate(x = df_detected_marker_mutations_in_ww_samples$Sample, by=list(PANGO_lineage=df_detected_marker_mutations_in_ww_samples$PANGO_lineage,site=df_detected_marker_mutations_in_ww_samples$site),FUN=function(x) length(unique(x)))
ggplot(data = df_plot_prevalence_PANGO_lineage_per_site,mapping=aes(x=site,y=x)) + geom_col(mapping=aes(fill=factor(PANGO_lineage,levels=names(sort(table(df_detected_marker_mutations_in_ww_samples$PANGO_lineage),decreasing = F)))),position = position_dodge(width = 0.9)) + xlab("site") + ylab(paste0("Prevalence (Number of samples in which \nthe lineage is detected;n=",nb_samples_original,")")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=12, angle = 60,hjust=1)) + labs(fill="PANGO lineages") + scale_y_continuous(limits=c(0,max(df_plot_prevalence_PANGO_lineage_per_site$x,na.rm=T)),breaks=seq(0,max(df_plot_prevalence_PANGO_lineage_per_site$x,na.rm=T)+2,2)) #+ scale_fill_manual(values=c("B.1.1.306"="skyblue","B.1.1.7"="red2","B.1.147"="blue3","B.1.160"="green3","B.1.36.17"="tomato", "B.40"="tan2","B.1.22"="orange","B.1.2"="brown"))
ggsave(filename = "prevalence_PANGO_lineage_per_site.png", path=output_workspace, width = 40, height = 20, units = "cm")


df_variants$month_year <- format(as.Date(df_variants$date,"%Y-%m-%d"),"%Y-%m")
v_lst_year_month <- c(paste0("2020-",c("01","02","03","04","05","06","07","08","09","10","11","12")),sort(intersect(format(as.Date(df_variants$date,"%Y-%m-%d"),"%Y-%m"), paste0("2021-",c("01","02","03","04","05","06","07","08","09","10","11","12")))))
v_lst_unique_locations <- sort(unique(df_variants$location))

#Data stratification analysis 
v_unique_site_type <- sort(unique(df_sample_stratifications_of_interest$Site_type))
v_unique_site_type <- v_unique_site_type[v_unique_site_type!=""]
v_unique_Location_in_network <- sort(unique(df_sample_stratifications_of_interest$Location_in_network))
v_unique_Location_in_network <- v_unique_Location_in_network[v_unique_Location_in_network!=""]
v_unique_Collector_type <- sort(unique(df_sample_stratifications_of_interest$Collector_type))
v_unique_Collector_type <- v_unique_Collector_type[v_unique_Collector_type!=""]
v_unique_site <- sort(unique(df_sample_stratifications_of_interest$site))
v_unique_site <- v_unique_site[v_unique_site!=""]
#Effect of site_type on WW detections
df_detected_marker_mutations_in_ww_samples$month_year <-  format(as.Date(df_detected_marker_mutations_in_ww_samples$date,"%Y-%m-%d"),"%Y-%m")
df_unique_detections_PANGO_lin <- unique(df_detected_marker_mutations_in_ww_samples[,c("Sample","date","month_year","location","PANGO_lineage","str_signature_mutations","nb_signature_mutations","stringent_confidence_score","coverage_adjusted_confidence_score")])
df_unique_detections_PANGO_lin <- subset(df_unique_detections_PANGO_lin, date != "2021-02-29")
  #Add site_type to dataframes
df_variants$site_type <- unname(vapply(X = df_variants$Sample,FUN = get_site_type_of_sample,FUN.VALUE = c("")))
df_unique_detections_PANGO_lin$site_type <- unname(vapply(X = df_unique_detections_PANGO_lin$Sample,FUN = get_site_type_of_sample,FUN.VALUE = c("")))
df_effect_of_site_type_on_WW_detections <- data.frame(month=rep(v_lst_year_month,times=length(v_unique_site_type)),site_type=rep(v_unique_site_type,each=length(v_lst_year_month)),nb_WW_samples=NA,nb_WW_detections=NA,WW_monthly_richness=NA,WW_monthly_Shannon_entropy=NA,WW_monthly_antilog_Shannon_entropy=NA,WW_monthly_Evenness=NA,stringsAsFactors = F)
for (i in 1:nrow(df_effect_of_site_type_on_WW_detections)){
  df_effect_of_site_type_on_WW_detections$nb_WW_samples[i] <- length(unique(subset(df_variants,(site_type==df_effect_of_site_type_on_WW_detections$site_type[i])&(!is.na(month_year))&(month_year==df_effect_of_site_type_on_WW_detections$month[i]))$Sample))
  df_effect_of_site_type_on_WW_detections$nb_WW_detections[i] <- sum((!is.na(df_unique_detections_PANGO_lin$month_year))&(df_unique_detections_PANGO_lin$month_year==df_effect_of_site_type_on_WW_detections$month[i])&(df_unique_detections_PANGO_lin$site_type==df_effect_of_site_type_on_WW_detections$site_type[i]),na.rm=T)
  df_effect_of_site_type_on_WW_detections$WW_monthly_richness[i] <- length(unique((subset(df_unique_detections_PANGO_lin,(site_type==df_effect_of_site_type_on_WW_detections$site_type[i])&(!is.na(month_year))&(month_year==df_effect_of_site_type_on_WW_detections$month[i])&(!is.na(PANGO_lineage))))$PANGO_lineage))
  df_effect_of_site_type_on_WW_detections$WW_monthly_Shannon_entropy[i] <- ifelse(test = df_effect_of_site_type_on_WW_detections$WW_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_unique_detections_PANGO_lin,(site_type==df_effect_of_site_type_on_WW_detections$site_type[i])&(!is.na(month_year))&(month_year==df_effect_of_site_type_on_WW_detections$month[i])))$PANGO_lineage),no=NA)
  df_effect_of_site_type_on_WW_detections$WW_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_effect_of_site_type_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_site_type_on_WW_detections$WW_monthly_Shannon_entropy[i],no=NA))
  df_effect_of_site_type_on_WW_detections$WW_monthly_Evenness[i] <- ifelse(test = df_effect_of_site_type_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_site_type_on_WW_detections$WW_monthly_Shannon_entropy[i]/(log2(df_effect_of_site_type_on_WW_detections$WW_monthly_richness[i])),no=NA)
}

Nb_detections_over_Seqrate_across_site_types_gg <- ggplot(data = df_effect_of_site_type_on_WW_detections,aes(x=as.factor(as.character(site_type)),y = nb_WW_detections/nb_WW_samples,fill=as.factor(as.character(site_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site type") + ylab("Number of lineage detections per month/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site type")
# Nb_detections_over_Seqrate_across_site_types_gg 
# ggsave(filename = "Nb_detections_over_Seqrate_per_month_across_site_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Richness_over_Seqrate_across_site_types_gg <- ggplot(data = df_effect_of_site_type_on_WW_detections,aes(x=as.factor(as.character(site_type)),y = WW_monthly_richness/nb_WW_samples,fill=as.factor(as.character(site_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site type") + ylab("Monthly richness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site type")
# Richness_over_Seqrate_across_site_types_gg
# ggsave(filename = "Richness_over_Seqrate_per_month_across_site_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Shannon_Entropy_over_Seqrate_across_site_types_gg <- ggplot(data = df_effect_of_site_type_on_WW_detections,aes(x=as.factor(as.character(site_type)),y = WW_monthly_Shannon_entropy/nb_WW_samples,fill=as.factor(as.character(site_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site type") + ylab("Monthly Shannon Entropy/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site type")
# Shannon_Entropy_over_Seqrate_across_site_types_gg
# ggsave(filename = "Shannon_Entropy_over_Seqrate_per_month_across_site_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Evenness_over_Seqrate_across_site_types_gg <- ggplot(data = df_effect_of_site_type_on_WW_detections,aes(x=as.factor(as.character(site_type)),y = WW_monthly_Evenness/nb_WW_samples,fill=as.factor(as.character(site_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site type") + ylab("Monthly Evenness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site type")
# Evenness_over_Seqrate_across_site_types_gg
# ggsave(filename = "Evenness_over_Seqrate_per_month_across_site_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

png(filename = paste0(output_workspace,"Effect_of_site_type_on_lineage_detections.png"),width = 40,height=30, units = "cm",res = 1200)
grid.arrange(Nb_detections_over_Seqrate_across_site_types_gg,Richness_over_Seqrate_across_site_types_gg,Shannon_Entropy_over_Seqrate_across_site_types_gg,Evenness_over_Seqrate_across_site_types_gg,nrow=2,ncol=2)
dev.off()

#Effect of Location_in_network on WW detections
  #Add Location_in_network to dataframes
df_variants$Location_in_network <- unname(vapply(X = df_variants$Sample,FUN = get_Location_in_network_of_sample,FUN.VALUE = c("")))
df_unique_detections_PANGO_lin$Location_in_network <- unname(vapply(X = df_unique_detections_PANGO_lin$Sample,FUN = get_Location_in_network_of_sample,FUN.VALUE = c("")))
df_effect_of_Location_in_network_on_WW_detections <- data.frame(month=rep(v_lst_year_month,times=length(v_unique_Location_in_network)),Location_in_network=rep(v_unique_Location_in_network,each=length(v_lst_year_month)),nb_WW_samples=NA,nb_WW_detections=NA,WW_monthly_richness=NA,WW_monthly_Shannon_entropy=NA,WW_monthly_antilog_Shannon_entropy=NA,WW_monthly_Evenness=NA,stringsAsFactors = F)
for (i in 1:nrow(df_effect_of_Location_in_network_on_WW_detections)){
  df_effect_of_Location_in_network_on_WW_detections$nb_WW_samples[i] <- length(unique(subset(df_variants,(Location_in_network==df_effect_of_Location_in_network_on_WW_detections$Location_in_network[i])&(!is.na(month_year))&(month_year==df_effect_of_Location_in_network_on_WW_detections$month[i]))$Sample))
  df_effect_of_Location_in_network_on_WW_detections$nb_WW_detections[i] <- sum((!is.na(df_unique_detections_PANGO_lin$month_year))&(df_unique_detections_PANGO_lin$month_year==df_effect_of_Location_in_network_on_WW_detections$month[i])&(df_unique_detections_PANGO_lin$Location_in_network==df_effect_of_Location_in_network_on_WW_detections$Location_in_network[i]),na.rm=T)
  df_effect_of_Location_in_network_on_WW_detections$WW_monthly_richness[i] <- length(unique((subset(df_unique_detections_PANGO_lin,(Location_in_network==df_effect_of_Location_in_network_on_WW_detections$Location_in_network[i])&(!is.na(month_year))&(month_year==df_effect_of_Location_in_network_on_WW_detections$month[i])&(!is.na(PANGO_lineage))))$PANGO_lineage))
  df_effect_of_Location_in_network_on_WW_detections$WW_monthly_Shannon_entropy[i] <- ifelse(test = df_effect_of_Location_in_network_on_WW_detections$WW_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_unique_detections_PANGO_lin,(Location_in_network==df_effect_of_Location_in_network_on_WW_detections$Location_in_network[i])&(!is.na(month_year))&(month_year==df_effect_of_Location_in_network_on_WW_detections$month[i])))$PANGO_lineage),no=NA)
  df_effect_of_Location_in_network_on_WW_detections$WW_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_effect_of_Location_in_network_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_Location_in_network_on_WW_detections$WW_monthly_Shannon_entropy[i],no=NA))
  df_effect_of_Location_in_network_on_WW_detections$WW_monthly_Evenness[i] <- ifelse(test = df_effect_of_Location_in_network_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_Location_in_network_on_WW_detections$WW_monthly_Shannon_entropy[i]/(log2(df_effect_of_Location_in_network_on_WW_detections$WW_monthly_richness[i])),no=NA)
}

Nb_detections_over_Seqrate_across_Locations_in_network_gg <- ggplot(data = df_effect_of_Location_in_network_on_WW_detections,aes(x=as.factor(as.character(Location_in_network)),y = nb_WW_detections/nb_WW_samples,fill=as.factor(as.character(Location_in_network)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Location in network") + ylab("Number of lineage detections per month/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Location in network")
# Nb_detections_over_Seqrate_across_Locations_in_network_gg 
# ggsave(filename = "Nb_detections_over_Seqrate_per_month_across_Locations_in_network.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Richness_over_Seqrate_across_Locations_in_network_gg <- ggplot(data = df_effect_of_Location_in_network_on_WW_detections,aes(x=as.factor(as.character(Location_in_network)),y = WW_monthly_richness/nb_WW_samples,fill=as.factor(as.character(Location_in_network)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Location in network") + ylab("Monthly richness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Location in network")
# Richness_over_Seqrate_across_Locations_in_network_gg
# ggsave(filename = "Richness_over_Seqrate_per_month_across_Locations_in_network.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Shannon_Entropy_over_Seqrate_across_Locations_in_network_gg <- ggplot(data = df_effect_of_Location_in_network_on_WW_detections,aes(x=as.factor(as.character(Location_in_network)),y = WW_monthly_Shannon_entropy/nb_WW_samples,fill=as.factor(as.character(Location_in_network)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Location in network") + ylab("Monthly Shannon Entropy/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Location in network")
# Shannon_Entropy_over_Seqrate_across_Locations_in_network_gg
# ggsave(filename = "Shannon_Entropy_over_Seqrate_per_month_across_Locations_in_network.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Evenness_over_Seqrate_across_Locations_in_network_gg <- ggplot(data = df_effect_of_Location_in_network_on_WW_detections,aes(x=as.factor(as.character(Location_in_network)),y = WW_monthly_Evenness/nb_WW_samples,fill=as.factor(as.character(Location_in_network)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Location in network") + ylab("Monthly Evenness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Location in network")
# Evenness_over_Seqrate_across_Locations_in_network_gg
# ggsave(filename = "Evenness_over_Seqrate_per_month_across_Locations_in_network.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

png(filename = paste0(output_workspace,"Effect_of_Location_in_network_on_lineage_detections.png"),width = 40,height=30, units = "cm",res = 1200)
grid.arrange(Nb_detections_over_Seqrate_across_Locations_in_network_gg,Richness_over_Seqrate_across_Locations_in_network_gg,Shannon_Entropy_over_Seqrate_across_Locations_in_network_gg,Evenness_over_Seqrate_across_Locations_in_network_gg,nrow=2,ncol=2)
dev.off()

#Effect of Collector_type on WW detections
  #Add Collector_type to dataframes
df_variants$Collector_type <- unname(vapply(X = df_variants$Sample,FUN = get_collector_type_of_sample,FUN.VALUE = c("")))
df_unique_detections_PANGO_lin$Collector_type <- unname(vapply(X = df_unique_detections_PANGO_lin$Sample,FUN = get_collector_type_of_sample,FUN.VALUE = c("")))
df_effect_of_Collector_type_on_WW_detections <- data.frame(month=rep(v_lst_year_month,times=length(v_unique_Collector_type)),Collector_type=rep(v_unique_Collector_type,each=length(v_lst_year_month)),nb_WW_samples=NA,nb_WW_detections=NA,WW_monthly_richness=NA,WW_monthly_Shannon_entropy=NA,WW_monthly_antilog_Shannon_entropy=NA,WW_monthly_Evenness=NA,stringsAsFactors = F)
for (i in 1:nrow(df_effect_of_Collector_type_on_WW_detections)){
  df_effect_of_Collector_type_on_WW_detections$nb_WW_samples[i] <- length(unique(subset(df_variants,(Collector_type==df_effect_of_Collector_type_on_WW_detections$Collector_type[i])&(!is.na(month_year))&(month_year==df_effect_of_Collector_type_on_WW_detections$month[i]))$Sample))
  df_effect_of_Collector_type_on_WW_detections$nb_WW_detections[i] <- sum((!is.na(df_unique_detections_PANGO_lin$month_year))&(df_unique_detections_PANGO_lin$month_year==df_effect_of_Collector_type_on_WW_detections$month[i])&(df_unique_detections_PANGO_lin$Collector_type==df_effect_of_Collector_type_on_WW_detections$Collector_type[i]),na.rm=T)
  df_effect_of_Collector_type_on_WW_detections$WW_monthly_richness[i] <- length(unique((subset(df_unique_detections_PANGO_lin,(Collector_type==df_effect_of_Collector_type_on_WW_detections$Collector_type[i])&(!is.na(month_year))&(month_year==df_effect_of_Collector_type_on_WW_detections$month[i])&(!is.na(PANGO_lineage))))$PANGO_lineage))
  df_effect_of_Collector_type_on_WW_detections$WW_monthly_Shannon_entropy[i] <- ifelse(test = df_effect_of_Collector_type_on_WW_detections$WW_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_unique_detections_PANGO_lin,(Collector_type==df_effect_of_Collector_type_on_WW_detections$Collector_type[i])&(!is.na(month_year))&(month_year==df_effect_of_Collector_type_on_WW_detections$month[i])))$PANGO_lineage),no=NA)
  df_effect_of_Collector_type_on_WW_detections$WW_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_effect_of_Collector_type_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_Collector_type_on_WW_detections$WW_monthly_Shannon_entropy[i],no=NA))
  df_effect_of_Collector_type_on_WW_detections$WW_monthly_Evenness[i] <- ifelse(test = df_effect_of_Collector_type_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_Collector_type_on_WW_detections$WW_monthly_Shannon_entropy[i]/(log2(df_effect_of_Collector_type_on_WW_detections$WW_monthly_richness[i])),no=NA)
}

Nb_detections_over_Seqrate_across_Collector_types_gg <- ggplot(data = df_effect_of_Collector_type_on_WW_detections,aes(x=as.factor(as.character(Collector_type)),y = nb_WW_detections/nb_WW_samples,fill=as.factor(as.character(Collector_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Collector type") + ylab("Number of lineage detections per month/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Collector type")
# Nb_detections_over_Seqrate_across_Collector_types_gg 
# ggsave(filename = "Nb_detections_over_Seqrate_per_month_across_Collector_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Richness_over_Seqrate_across_Collector_types_gg <- ggplot(data = df_effect_of_Collector_type_on_WW_detections,aes(x=as.factor(as.character(Collector_type)),y = WW_monthly_richness/nb_WW_samples,fill=as.factor(as.character(Collector_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Collector type") + ylab("Monthly richness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Collector type")
# Richness_over_Seqrate_across_Collector_types_gg
# ggsave(filename = "Richness_over_Seqrate_per_month_across_Collector_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Shannon_Entropy_over_Seqrate_across_Collector_types_gg <- ggplot(data = df_effect_of_Collector_type_on_WW_detections,aes(x=as.factor(as.character(Collector_type)),y = WW_monthly_Shannon_entropy/nb_WW_samples,fill=as.factor(as.character(Collector_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Collector type") + ylab("Monthly Shannon Entropy/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Collector type")
# Shannon_Entropy_over_Seqrate_across_Collector_types_gg
# ggsave(filename = "Shannon_Entropy_over_Seqrate_per_month_across_Collector_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Evenness_over_Seqrate_across_Collector_types_gg <- ggplot(data = df_effect_of_Collector_type_on_WW_detections,aes(x=as.factor(as.character(Collector_type)),y = WW_monthly_Evenness/nb_WW_samples,fill=as.factor(as.character(Collector_type)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Collector type") + ylab("Monthly Evenness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Collector type")
# Evenness_over_Seqrate_across_Collector_types_gg
# ggsave(filename = "Evenness_over_Seqrate_per_month_across_Collector_types.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

png(filename = paste0(output_workspace,"Effect_of_Collector_type_on_lineage_detections.png"),width = 40,height=30, units = "cm",res = 1200)
grid.arrange(Nb_detections_over_Seqrate_across_Collector_types_gg,Richness_over_Seqrate_across_Collector_types_gg,Shannon_Entropy_over_Seqrate_across_Collector_types_gg,Evenness_over_Seqrate_across_Collector_types_gg,nrow=2,ncol=2)
dev.off()

#Effect of site on WW detections
  #Add site to dataframes
df_variants$site <- unname(vapply(X = df_variants$Sample,FUN = get_site_of_sample,FUN.VALUE = c("")))
df_unique_detections_PANGO_lin$site <- unname(vapply(X = df_unique_detections_PANGO_lin$Sample,FUN = get_site_of_sample,FUN.VALUE = c("")))
df_effect_of_site_on_WW_detections <- data.frame(month=rep(v_lst_year_month,times=length(v_unique_site)),site=rep(v_unique_site,each=length(v_lst_year_month)),nb_WW_samples=NA,nb_WW_detections=NA,WW_monthly_richness=NA,WW_monthly_Shannon_entropy=NA,WW_monthly_antilog_Shannon_entropy=NA,WW_monthly_Evenness=NA,stringsAsFactors = F)
for (i in 1:nrow(df_effect_of_site_on_WW_detections)){
  df_effect_of_site_on_WW_detections$nb_WW_samples[i] <- length(unique(subset(df_variants,(site==df_effect_of_site_on_WW_detections$site[i])&(!is.na(month_year))&(month_year==df_effect_of_site_on_WW_detections$month[i]))$Sample))
  df_effect_of_site_on_WW_detections$nb_WW_detections[i] <- sum((!is.na(df_unique_detections_PANGO_lin$month_year))&(df_unique_detections_PANGO_lin$month_year==df_effect_of_site_on_WW_detections$month[i])&(df_unique_detections_PANGO_lin$site==df_effect_of_site_on_WW_detections$site[i]),na.rm=T)
  df_effect_of_site_on_WW_detections$WW_monthly_richness[i] <- length(unique((subset(df_unique_detections_PANGO_lin,(site==df_effect_of_site_on_WW_detections$site[i])&(!is.na(month_year))&(month_year==df_effect_of_site_on_WW_detections$month[i])&(!is.na(PANGO_lineage))))$PANGO_lineage))
  df_effect_of_site_on_WW_detections$WW_monthly_Shannon_entropy[i] <- ifelse(test = df_effect_of_site_on_WW_detections$WW_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_unique_detections_PANGO_lin,(site==df_effect_of_site_on_WW_detections$site[i])&(!is.na(month_year))&(month_year==df_effect_of_site_on_WW_detections$month[i])))$PANGO_lineage),no=NA)
  df_effect_of_site_on_WW_detections$WW_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_effect_of_site_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_site_on_WW_detections$WW_monthly_Shannon_entropy[i],no=NA))
  df_effect_of_site_on_WW_detections$WW_monthly_Evenness[i] <- ifelse(test = df_effect_of_site_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_site_on_WW_detections$WW_monthly_Shannon_entropy[i]/(log2(df_effect_of_site_on_WW_detections$WW_monthly_richness[i])),no=NA)
}

Nb_detections_over_Seqrate_across_sites_gg <- ggplot(data = df_effect_of_site_on_WW_detections,aes(x=as.factor(as.character(site)),y = nb_WW_detections/nb_WW_samples,fill=as.factor(as.character(site)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site") + ylab("Number of lineage detections per month/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site")
# Nb_detections_over_Seqrate_across_sites_gg 
# ggsave(filename = "Nb_detections_over_Seqrate_per_month_across_sites.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Richness_over_Seqrate_across_sites_gg <- ggplot(data = df_effect_of_site_on_WW_detections,aes(x=as.factor(as.character(site)),y = WW_monthly_richness/nb_WW_samples,fill=as.factor(as.character(site)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site") + ylab("Monthly richness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site")
# Richness_over_Seqrate_across_sites_gg
# ggsave(filename = "Richness_over_Seqrate_per_month_across_sites.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Shannon_Entropy_over_Seqrate_across_sites_gg <- ggplot(data = df_effect_of_site_on_WW_detections,aes(x=as.factor(as.character(site)),y = WW_monthly_Shannon_entropy/nb_WW_samples,fill=as.factor(as.character(site)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site") + ylab("Monthly Shannon Entropy/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site")
# Shannon_Entropy_over_Seqrate_across_sites_gg
# ggsave(filename = "Shannon_Entropy_over_Seqrate_per_month_across_sites.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

Evenness_over_Seqrate_across_sites_gg <- ggplot(data = df_effect_of_site_on_WW_detections,aes(x=as.factor(as.character(site)),y = WW_monthly_Evenness/nb_WW_samples,fill=as.factor(as.character(site)))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("Site") + ylab("Monthly Evenness/\nNumber of sequenced samples per month") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_text(angle=60,hjust=1)) + stat_compare_means(method = "kruskal") + labs(fill="Site")
# Evenness_over_Seqrate_across_sites_gg
# ggsave(filename = "Evenness_over_Seqrate_per_month_across_sites.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

png(filename = paste0(output_workspace,"Effect_of_site_on_lineage_detections.png"),width = 40,height=30, units = "cm",res = 1200)
grid.arrange(Nb_detections_over_Seqrate_across_sites_gg,Richness_over_Seqrate_across_sites_gg,Shannon_Entropy_over_Seqrate_across_sites_gg,Evenness_over_Seqrate_across_sites_gg,nrow=2,ncol=2)
dev.off()

#Varition partitioning analysis
df_effect_of_5factors_on_WW_detections <- expand.grid(v_lst_year_month,v_unique_Collector_type,v_unique_Location_in_network,v_unique_site,v_unique_site_type)
colnames(df_effect_of_5factors_on_WW_detections) <- c("month","Collector_type","Location_in_network","site","site_type")
df_effect_of_5factors_on_WW_detections$nb_WW_samples <- NA
df_effect_of_5factors_on_WW_detections$nb_WW_detections <- NA
df_effect_of_5factors_on_WW_detections$WW_monthly_richness <- NA
df_effect_of_5factors_on_WW_detections$WW_monthly_Shannon_entropy <- NA
df_effect_of_5factors_on_WW_detections$WW_monthly_antilog_Shannon_entropy <- NA
df_effect_of_5factors_on_WW_detections$WW_monthly_Evenness <- NA
df_effect_of_5factors_on_WW_detections <- unique(df_effect_of_5factors_on_WW_detections)

for (i in 1:nrow(df_effect_of_5factors_on_WW_detections)){
  df_effect_of_5factors_on_WW_detections$nb_WW_samples[i] <- length(unique(subset(df_variants,(Collector_type==df_effect_of_5factors_on_WW_detections$Collector_type[i])&(Location_in_network==df_effect_of_5factors_on_WW_detections$Location_in_network[i])&(site==df_effect_of_5factors_on_WW_detections$site[i])&(site_type==df_effect_of_5factors_on_WW_detections$site_type[i])&(!is.na(month_year))&(month_year==df_effect_of_5factors_on_WW_detections$month[i]))$Sample))
  df_effect_of_5factors_on_WW_detections$nb_WW_detections[i] <- sum((!is.na(df_unique_detections_PANGO_lin$month_year))&(df_unique_detections_PANGO_lin$month_year==df_effect_of_5factors_on_WW_detections$month[i])&(df_unique_detections_PANGO_lin$Collector_type==df_effect_of_5factors_on_WW_detections$Collector_type[i])&(df_unique_detections_PANGO_lin$Location_in_network==df_effect_of_5factors_on_WW_detections$Location_in_network[i])&(df_unique_detections_PANGO_lin$site==df_effect_of_5factors_on_WW_detections$site[i])&(df_unique_detections_PANGO_lin$site_type==df_effect_of_5factors_on_WW_detections$site_type[i]),na.rm=T)
  df_effect_of_5factors_on_WW_detections$WW_monthly_richness[i] <- length(unique((subset(df_unique_detections_PANGO_lin,(Collector_type==df_effect_of_5factors_on_WW_detections$Collector_type[i])&(Location_in_network==df_effect_of_5factors_on_WW_detections$Location_in_network[i])&(site==df_effect_of_5factors_on_WW_detections$site[i])&(site_type==df_effect_of_5factors_on_WW_detections$site_type[i])&(!is.na(month_year))&(month_year==df_effect_of_5factors_on_WW_detections$month[i])&(!is.na(PANGO_lineage))))$PANGO_lineage))
  df_effect_of_5factors_on_WW_detections$WW_monthly_Shannon_entropy[i] <- ifelse(test = df_effect_of_5factors_on_WW_detections$WW_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_unique_detections_PANGO_lin,(Collector_type==df_effect_of_5factors_on_WW_detections$Collector_type[i])&(Location_in_network==df_effect_of_5factors_on_WW_detections$Location_in_network[i])&(site==df_effect_of_5factors_on_WW_detections$site[i])&(site_type==df_effect_of_5factors_on_WW_detections$site_type[i])&(!is.na(month_year))&(month_year==df_effect_of_5factors_on_WW_detections$month[i])))$PANGO_lineage),no=NA)
  df_effect_of_5factors_on_WW_detections$WW_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_effect_of_5factors_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_5factors_on_WW_detections$WW_monthly_Shannon_entropy[i],no=NA))
  df_effect_of_5factors_on_WW_detections$WW_monthly_Evenness[i] <- ifelse(test = df_effect_of_5factors_on_WW_detections$WW_monthly_richness[i]>0,yes=df_effect_of_5factors_on_WW_detections$WW_monthly_Shannon_entropy[i]/(log2(df_effect_of_5factors_on_WW_detections$WW_monthly_richness[i])),no=NA)
}

# df_effect_of_5factors_on_WW_detections <- NULL
# i <- 1
# for (current_mont_year in v_lst_year_month){
#   for (current_Collector_type in v_unique_Collector_type){
#     for (current_Location_in_network in v_unique_Location_in_network){
#       for (current_site in v_unique_site){
#         for (current_site_type in v_unique_site_type){
#           df_effect_of_5factors_on_WW_detections <- rbind(df_effect_of_5factors_on_WW_detections,data.frame(month=current_mont_year,Collector_type=current_Collector_type,Location_in_network=current_Location_in_network,site=current_site,site_type=current_site_type,nb_WW_samples=NA,nb_WW_detections=NA,WW_monthly_richness=NA,WW_monthly_Shannon_entropy=NA,WW_monthly_antilog_Shannon_entropy=NA,WW_monthly_Evenness=NA,stringsAsFactors = F))
#           print(i)
#           i <- i+1
#         }
#       }
#     }
#   }
# }
  #variation partitioning of Number of Monthly detections
varp_nb_monthly_detections <- varpart (Y = df_effect_of_5factors_on_WW_detections$nb_WW_detections,~ month, ~ nb_WW_samples, ~ Collector_type, ~ site_type, data = df_effect_of_5factors_on_WW_detections)
png(filename = paste0(output_workspace,"Variation_partitioning_nb_monthly_detections.png"),width = 30,height=25, units = "cm",res = 1200)
plot(varp_nb_monthly_detections, digits = 3, Xnames = c("month","Sequencing rate","Collector type","Site type"))
dev.off()

#variation partitioning of Richness
varp_Richness <- varpart (Y = df_effect_of_5factors_on_WW_detections$WW_monthly_richness,~ month, ~ nb_WW_samples, ~ Collector_type, ~ site_type, data = df_effect_of_5factors_on_WW_detections)
png(filename = paste0(output_workspace,"Variation_partitioning_Richness.png"),width = 30,height=25, units = "cm",res = 1200)
plot(varp_Richness, digits = 3, Xnames = c("month","Sequencing rate","Collector type","Site type"))
dev.off()

# #variation partitioning of Shannon Entropy
# varp_Shannon_Entropy <- varpart (Y = df_effect_of_5factors_on_WW_detections$WW_monthly_Shannon_entropy,~ month, ~ nb_WW_samples, ~ Collector_type, ~ site_type, data = df_effect_of_5factors_on_WW_detections)
# png(filename = paste0(output_workspace,"Variation_partitioning_Shannon_Entropy.png"),width = 30,height=25, units = "cm",res = 1200)
# plot(varp_Shannon_Entropy, digits = 3, Xnames = c("month","Sequencing rate","Collector type","Site type"))
# dev.off()
# 
# #variation partitioning of Evenness
# varp_Evenness <- varpart (Y = df_effect_of_5factors_on_WW_detections$WW_monthly_Evenness,~ month, ~ nb_WW_samples, ~ Collector_type, ~ site_type, data = df_effect_of_5factors_on_WW_detections)
# png(filename = paste0(output_workspace,"Variation_partitioning_Evenness.png"),width = 30,height=25, units = "cm",res = 1200)
# plot(varp_Evenness, digits = 3, Xnames = c("month","Sequencing rate","Collector type","Site type"))
# dev.off()


#Stacked barplot of the number of marker mutations per ORF colored by lineage
v_orfs_of_marker_mutations <- unname(vapply(as.integer(substr(v_lst_marker_mutations,2,nchar(v_lst_marker_mutations)-1)),FUN = find_ORF_of_mutation,FUN.VALUE = c("")))
df_freq_orfs_of_marker_mutations <- NULL
for (current_orf in v_orfs){
  v_lineages_of_current_orf <- unname(v_lineage_marker_mutations[v_orfs_of_marker_mutations==current_orf])
  v_lineages_of_current_orf <- v_lineages_of_current_orf[!is.na(v_lineages_of_current_orf)]
  for (current_lineage in unique(v_lineages_of_current_orf)){
    df_freq_orfs_of_marker_mutations <- rbind(df_freq_orfs_of_marker_mutations,data.frame(ORF=current_orf,lineage=current_lineage,number_of_marker_mutations=sum(v_lineages_of_current_orf==current_lineage)))
  }
}
ggplot(data = df_freq_orfs_of_marker_mutations,mapping=aes(x=factor(ORF,levels=v_orfs),y=number_of_marker_mutations,fill=factor(lineage,levels=sort(unique(unname(v_lineage_marker_mutations)))) )) + geom_col() + xlab("Genomic region") + ylab("Number of marker mutations") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none",legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=10,angle = 60,hjust=1)) + labs(fill="Lineage")
ggsave(filename = "Number_of_marker_mutations_per_ORF_colored_by_lineage.png", path=output_workspace, width = 20, height = 15, units = "cm")
df_freq_orfs_of_marker_mutations$density_of_marker_mutations <- df_freq_orfs_of_marker_mutations$number_of_marker_mutations/v_orfs_length[as.character(df_freq_orfs_of_marker_mutations$ORF)]
ggplot(data = df_freq_orfs_of_marker_mutations,mapping=aes(x=factor(ORF,levels=v_orfs),y=density_of_marker_mutations,fill=factor(lineage,levels=sort(unique(unname(v_lineage_marker_mutations)))) )) + geom_col() + xlab("Genomic region") + ylab("Density of marker mutations\n(number of marker mutations/ORF length)") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none",legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=10,angle = 60,hjust=1)) + labs(fill="Lineage")
ggsave(filename = "Density_of_marker_mutations_per_ORF_colored_by_lineage.png", path=output_workspace, width = 20, height = 15, units = "cm")

df_nb_samples_per_location_across_time <- data.frame(location= rep(v_lst_unique_locations,each=length(v_lst_year_month)), time = rep(v_lst_year_month,length(v_lst_unique_locations)), Nb_samples=0, stringsAsFactors = F)
for(i in 1:nrow(df_nb_samples_per_location_across_time)){
  df_nb_samples_per_location_across_time$Nb_samples[i] <- sum(unname(vapply(X = lst_samples_original,FUN = function(x) (!is.na(v_samples_date[x]))&(!is.na(v_samples_location[x]))&grepl(pattern = df_nb_samples_per_location_across_time$time[i],x = v_samples_date[x],fixed = T)&grepl(pattern = df_nb_samples_per_location_across_time$location[i],x = v_samples_location[x],fixed = T),FUN.VALUE = c(F))))
}
df_nb_samples_per_location_across_time <- subset(df_nb_samples_per_location_across_time,(!is.na(location))&(location!="OTTAWA")&(location!="ON"))
ggplot(data = df_nb_samples_per_location_across_time,mapping=aes(x=factor(time,levels=v_lst_year_month),y=Nb_samples)) + geom_col(mapping=aes(fill=factor(location,levels=v_lst_unique_locations))) + xlab("Time") + ylab(paste0("Number of Wastewater samples with SNV(s)")) + theme_bw() + theme(axis.title = element_text(size=14),axis.text = element_text(size=14),legend.title = element_text(size=12),legend.text = element_text(size=12),axis.text.x = element_text(size=12, angle = 60,hjust=1),axis.text.y = element_text(size=12)) + labs(fill="Location") + scale_y_continuous(limits=c(0,max(df_nb_samples_per_location_across_time$Nb_samples*1.5,na.rm=T)),breaks=seq(0,max(df_nb_samples_per_location_across_time$Nb_samples*1.5,na.rm=T)+20,20)) + scale_fill_manual(values=c("QUEBEC"="blue3","MONTREAL"="green3","OTTAWA"="red","ON"="grey","LAVAL"="violet"))
ggsave(filename = "Nb_samples_per_location_across_time.png", path=output_workspace, width = 20, height = 15, units = "cm",dpi=1200)

#Concordance between WW and clinical sampling
#Add sampling site information
df_metadata_LSPQ_samples <- read.csv2(file = paste0(output_workspace,"LATEST_REPORT.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_metadata_LSPQ_samples <- subset(df_metadata_LSPQ_samples,(run.type%in%c("production"))&(sample.type=="clinical"))
df_HospCenter <- read.csv2(file = paste0(output_workspace,"ListeCH_prefix.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
v_site_to_RSS <- df_HospCenter$RSS
names(v_site_to_RSS) <- df_HospCenter$Prefix
df_metadata_LSPQ_samples$RSS <- ifelse(test=is.na(df_metadata_LSPQ_samples$site),yes=NA,no=v_site_to_RSS[df_metadata_LSPQ_samples$site])
df_metadata_LSPQ_samples$RSS <- toupper(df_metadata_LSPQ_samples$RSS)
for (i in 1:nrow(df_metadata_LSPQ_samples)){
  if (grepl(pattern = "B.1.617",x = df_metadata_LSPQ_samples$PANGOLIN[i],fixed = T)){
    df_metadata_LSPQ_samples$PANGOLIN[i] <- "B.1.617.X"
  }
}
#Determining the thresholds confidence score and time gaps
v_tested_min_stringent_confidence <- seq(0,1,0.05)
v_tested_max_time_gap <- 0:8
df_concordance_score_PANGO_lineages_inference <- NULL
for (current_min_stringent_confidence in v_tested_min_stringent_confidence){
  #only consider lineages for which we have both WW and clinical sampling data != outbreak
  current_df_inferred_pango_lineages <- subset(df_unique_detections_PANGO_lin,(stringent_confidence_score >= current_min_stringent_confidence)&(!PANGO_lineage%in%setdiff(subset(df_metadata_LSPQ_samples,run.type=="outbreak")$PANGOLIN,subset(df_metadata_LSPQ_samples,run.type!="outbreak")$PANGOLIN)))
  current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling <- NA
  for (i in 1:nrow(current_df_inferred_pango_lineages)){
    tryCatch(expr = {current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling[i] <- min(abs(as.integer(as.Date(subset(df_metadata_LSPQ_samples,(PANGOLIN==current_df_inferred_pango_lineages$PANGO_lineage[i])&(RSS==current_df_inferred_pango_lineages$location[i]))$DATE_PRELEV)-as.Date(current_df_inferred_pango_lineages$date[i]))),na.rm=T)},error = function(e) print(e))
  }
  for (current_max_time_gap in v_tested_max_time_gap){
    df_concordance_score_PANGO_lineages_inference <- rbind(df_concordance_score_PANGO_lineages_inference,data.frame(min_stringent_confidence=current_min_stringent_confidence,max_time_gap_tolerated=current_max_time_gap,concordance_score=((sum((!is.na(current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling))&(current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling<=current_max_time_gap)))/(nrow(current_df_inferred_pango_lineages))),nb_detections=nrow(current_df_inferred_pango_lineages)))
  }
}
df_concordance_score_PANGO_lineages_inference$scaled_number_of_detection <- scale(df_concordance_score_PANGO_lineages_inference$nb_detections,T,T)
df_concordance_score_PANGO_lineages_inference$label_facet <- paste0("Max. time difference = ",df_concordance_score_PANGO_lineages_inference$max_time_gap_tolerated)
#plot sensitivity analysis (how does the concordance between WW and clinical sampling depends of the maximum time gap tolerated and the minimum confidence score?)
ggplot(data = df_concordance_score_PANGO_lineages_inference,mapping=aes(x=min_stringent_confidence,y=concordance_score)) + geom_point(aes(size=nb_detections)) + xlab("Minimum confidence score") + ylab(paste0("Concordance between Wastewater\nand clinical sampling data")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1),legend.position = "right") + scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.1)) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.05)) + facet_wrap(~label_facet,nrow=3) + labs(size="Number of PANGO lineage\n detections in Wastewater")
ggsave(filename = "Sensitvity_analysis_for_concordance_WW_and_clinical_sampling.png", path=output_workspace, width = 30, height = 30, units = "cm",dpi=1200)

v_min_nb_signature_muts_detected <- 1:10
df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts <- NULL
for (current_min_nb_signature_muts_detected in v_min_nb_signature_muts_detected){
  #only consider lineages for which we have both WW and clinical sampling data != outbreak
  current_df_inferred_pango_lineages <- subset(df_unique_detections_PANGO_lin,(nb_signature_mutations >= current_min_nb_signature_muts_detected)&(!PANGO_lineage%in%setdiff(subset(df_metadata_LSPQ_samples,run.type=="outbreak")$PANGOLIN,subset(df_metadata_LSPQ_samples,run.type!="outbreak")$PANGOLIN)))
  current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling <- NA
  for (i in 1:nrow(current_df_inferred_pango_lineages)){
    tryCatch(expr = {current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling[i] <- min(abs(as.integer(as.Date(subset(df_metadata_LSPQ_samples,(PANGOLIN==current_df_inferred_pango_lineages$PANGO_lineage[i])&(RSS==current_df_inferred_pango_lineages$location[i]))$DATE_PRELEV)-as.Date(current_df_inferred_pango_lineages$date[i]))),na.rm=T)},error = function(e) print(e))
  }
  for (current_max_time_gap in v_tested_max_time_gap){
    df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts <- rbind(df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts,data.frame(min_nb_signature_muts_detected=current_min_nb_signature_muts_detected,max_time_gap_tolerated=current_max_time_gap,concordance_score=((sum((!is.na(current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling))&(current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling<=current_max_time_gap)))/(nrow(current_df_inferred_pango_lineages))),nb_detections=nrow(current_df_inferred_pango_lineages)))
  }
}
df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts$scaled_number_of_detection <- scale(df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts$nb_detections,T,T)
df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts$label_facet <- paste0("Max. time difference = ",df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts$max_time_gap_tolerated)
#plot sensitivity analysis (how does the concordance between WW and clinical sampling depends of the maximum time gap tolerated and the minimum number of signature mutations?)
ggplot(data = df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts,mapping=aes(x=min_nb_signature_muts_detected,y=concordance_score)) + geom_point(aes(size=nb_detections)) + xlab("Minimum number of\nsignature mutations detected") + ylab(paste0("Concordance between Wastewater\nand clinical sampling data")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1),legend.position = "right") +scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.1)) + scale_x_continuous(limits = c(0,10),breaks=0:10) + facet_wrap(~label_facet,nrow=3) + labs(size="Number of PANGO lineage\n detections in Wastewater")
ggsave(filename = "Sensitvity_analysis_for_concordance_WW_and_clinical_sampling_WITH_abs_nb_sig_muts_detected.png", path=output_workspace, width = 30, height = 30, units = "cm",dpi=1200)

df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts$label_facet2 <- factor(paste0("Min. number of detected signature SNVs = ",df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts$min_nb_signature_muts_detected),levels=paste0("Min. number of detected signature SNVs = ",1:10))
#v2
ggplot(data = df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts,mapping=aes(x=max_time_gap_tolerated,y=concordance_score)) + geom_point(aes(size=nb_detections)) + xlab("Max. time difference") + ylab(paste0("Concordance between Wastewater\nand clinical sampling data")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1),legend.position = "right") +scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.1)) + scale_x_continuous(limits = c(0,8),breaks=0:8) + facet_wrap(~label_facet2,nrow=5) + labs(size="Number of PANGO lineage\n detections in Wastewater")
ggsave(filename = "Sensitvity_analysis_for_concordance_WW_and_clinical_sampling_WITH_abs_nb_sig_muts_detected_v_2_0.png", path=output_workspace, width = 30, height = 30, units = "cm",dpi=1200)

#Permutation tests to get p-value 
df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts_and_pvalues <- subset(df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts,(min_nb_signature_muts_detected%in%(3:8))&(max_time_gap_tolerated%in%(5:7))) #df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts
df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts_and_pvalues$p.value <- NA
for (indx_row in 1:nrow(df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts_and_pvalues)){
  v_the_concordance_score_at_selected_thresholds <- subset(df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts,(min_nb_signature_muts_detected==current_min_nb_signature_muts_detected)&(max_time_gap_tolerated==current_max_time_gap))$concordance_score
  nb_perms <- 999
  lst_splits <- split(1:nb_perms, ceiling(seq_along(1:nb_perms)/(nb_perms/nb_cores)))
  the_f_parallel <- function(i_cl){
    the_vec<- lst_splits[[i_cl]]
    current_df_concordance_score_PANGO_lineages_inference_after_permutation <- NULL
    z <- 1
    for (k in the_vec){
      permuted_dates_in_df_unique_detections_PANGO_lin <- df_unique_detections_PANGO_lin
      permuted_dates_in_df_unique_detections_PANGO_lin$date <- sample(x=df_unique_detections_PANGO_lin$date,size = nrow(df_unique_detections_PANGO_lin),replace = F)
      permuted_dates_in_df_unique_detections_PANGO_lin$location <- sample(x=df_unique_detections_PANGO_lin$location,size = nrow(df_unique_detections_PANGO_lin),replace = F)
      #only consider lineages for which we have both WW and clinical sampling data != outbreak
      current_df_inferred_pango_lineages <- subset(permuted_dates_in_df_unique_detections_PANGO_lin,(nb_signature_mutations >= current_min_nb_signature_muts_detected)&(!PANGO_lineage%in%setdiff(subset(df_metadata_LSPQ_samples,run.type=="outbreak")$PANGOLIN,subset(df_metadata_LSPQ_samples,run.type!="outbreak")$PANGOLIN)))
      current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling <- NA
      for (i in 1:nrow(current_df_inferred_pango_lineages)){
        tryCatch(expr = {current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling[i] <- min(abs(as.integer(as.Date(subset(df_metadata_LSPQ_samples,(PANGOLIN==current_df_inferred_pango_lineages$PANGO_lineage[i])&(RSS==current_df_inferred_pango_lineages$location[i]))$DATE_PRELEV)-as.Date(current_df_inferred_pango_lineages$date[i]))),na.rm=T)},error = function(e) print(e))
      }
      current_df_concordance_score_PANGO_lineages_inference_after_permutation <- rbind(current_df_concordance_score_PANGO_lineages_inference_after_permutation,data.frame(min_nb_signature_muts_detected=current_min_nb_signature_muts_detected,max_time_gap_tolerated=current_max_time_gap,concordance_score=((sum((!is.na(current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling))&(current_df_inferred_pango_lineages$best_time_gap_with_clinical_sampling<=current_max_time_gap)))/(nrow(current_df_inferred_pango_lineages))),nb_detections=nrow(current_df_inferred_pango_lineages),ID_permutation=k))
      if (z%%5==0){
        print(paste0(z," permutations done out of ",length(the_vec),"!"))
      }
      z <- z + 1
    }
    
    return(current_df_concordance_score_PANGO_lineages_inference_after_permutation)
  }
  cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_Concordance_analysis.txt"))
  registerDoParallel(cl)
  current_threshold_set_df_concordance_score_PANGO_lineages_inference_after_permutation <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel(i_cl)
  stopCluster(cl)
  
  df_concordance_score_PANGO_lineages_inference_with_abs_nb_sig_muts_and_pvalues$p.value[indx_row] <- sum(current_threshold_set_df_concordance_score_PANGO_lineages_inference_after_permutation$concordance_score >= v_the_concordance_score_at_selected_thresholds)/nrow(current_threshold_set_df_concordance_score_PANGO_lineages_inference_after_permutation)
}

#Determine which method detected the first arrival of each lineage in Qc (DESPITE HIGHER SEQUENCING RATE)
v_lineages_to_check <- union(df_detected_marker_mutations_in_ww_samples$PANGO_lineage,df_metadata_LSPQ_samples$PANGOLIN)[(unname(vapply(X = union(df_detected_marker_mutations_in_ww_samples$PANGO_lineage,df_metadata_LSPQ_samples$PANGOLIN),FUN = function(x) grepl(pattern = ".",x = x,fixed = T),FUN.VALUE = F)))&(!is.na(union(df_detected_marker_mutations_in_ww_samples$PANGO_lineage,df_metadata_LSPQ_samples$PANGOLIN)))]
v_WW_first_arrival <- rep(NA,length(v_lineages_to_check))
v_clinical_sampling_first_arrival <- rep(NA,length(v_lineages_to_check))
for (i in 1:length(v_lineages_to_check)){
  v_WW_first_arrival[i] <- min(subset(df_detected_marker_mutations_in_ww_samples,(nb_signature_mutations>=3)&(PANGO_lineage==v_lineages_to_check[i]))$date,na.rm=T)
  v_clinical_sampling_first_arrival[i] <- min(subset(df_metadata_LSPQ_samples,(run.type!="outbreak")&(PANGOLIN==v_lineages_to_check[i]))$DATE_PRELEV,na.rm=T)
}

percentage_first_lineage_occurences_detected_in_WW <- 100*sum(((v_WW_first_arrival<v_clinical_sampling_first_arrival)&(!is.na(v_WW_first_arrival))&(!is.na(v_clinical_sampling_first_arrival)))|((!is.na(v_WW_first_arrival))&(is.na(v_clinical_sampling_first_arrival))),na.rm = T)/length(v_lineages_to_check)
v_which_lineages_detected_first_in_WW <- v_lineages_to_check[((v_WW_first_arrival<v_clinical_sampling_first_arrival)&(!is.na(v_WW_first_arrival))&(!is.na(v_clinical_sampling_first_arrival)))|((!is.na(v_WW_first_arrival))&(is.na(v_clinical_sampling_first_arrival)))] 
v_date_detection_lineages_detected_first_in_WW <- v_WW_first_arrival[((v_WW_first_arrival<v_clinical_sampling_first_arrival)&(!is.na(v_WW_first_arrival))&(!is.na(v_clinical_sampling_first_arrival)))|((!is.na(v_WW_first_arrival))&(is.na(v_clinical_sampling_first_arrival)))] 

v_which_lineages_detected_first_in_clinical_sampling <- v_lineages_to_check[((v_clinical_sampling_first_arrival<v_WW_first_arrival)&(!is.na(v_clinical_sampling_first_arrival))&(!is.na(v_WW_first_arrival)))|((!is.na(v_clinical_sampling_first_arrival))&(is.na(v_WW_first_arrival)))] 
v_date_detection_lineages_detected_first_in_clinical_sampling <- v_clinical_sampling_first_arrival[((v_clinical_sampling_first_arrival<v_WW_first_arrival)&(!is.na(v_clinical_sampling_first_arrival))&(!is.na(v_WW_first_arrival)))|((!is.na(v_clinical_sampling_first_arrival))&(is.na(v_WW_first_arrival)))]

#Venn Diagram detected lineages in WW vs Clinical sampling
venn.diagram(
  x = list(v_lineages_to_check[(!is.na(v_clinical_sampling_first_arrival))],v_lineages_to_check[(!is.na(v_WW_first_arrival))]),
  category.names = c("Clinical" , "WW") ,
  filename =  paste0(output_workspace,"Venn_diagram_lineages_detections_Clinical_vs_WW.png"),
  output = TRUE ,
  imagetype="png" ,
  fill =  c(alpha("tomato",0.6), alpha('deepskyblue',0.6))
)

#compute Shannon entropy
get_entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

#Relationship between number of detections and sequencing rate (nb_samples_per_month) 
df_variants$month_year <- format(as.Date(df_variants$date,"%Y-%m-%d"),"%Y-%m")
df_unique_detections_PANGO_lin$month_year <- format(as.Date(df_unique_detections_PANGO_lin$date,"%Y-%m-%d"),"%Y-%m")
df_clinical_detections <- subset(df_metadata_LSPQ_samples,(run.type!="outbreak")&(DATE_PRELEV>="2020-01-01"))
df_clinical_detections$month_year <- format(as.Date(df_clinical_detections$DATE_PRELEV,"%Y-%m-%d"),"%Y-%m")
df_nb_detections_PANGO_lineage_vs_seqrate <- data.frame(month=c(paste0("2020-",c("01","02","03","04","05","06","07","08","09","10","11","12")),paste0("2021-",c("01","02","03","04","05","06"))),nb_WW_samples=NA,nb_clinical_samples=NA,nb_WW_detections=NA,nb_clinical_detections=NA,WW_monthly_richness=NA,WW_monthly_Shannon_entropy=NA,WW_monthly_antilog_Shannon_entropy=NA,WW_monthly_Evenness=NA,Clinical_sampling_monthly_richness=NA,Clinical_sampling_monthly_Shannon_entropy=NA,Clinical_sampling_monthly_antilog_Shannon_entropy=NA,Clinical_sampling_monthly_Evenness=NA,stringsAsFactors = F)
for (i in 1:nrow(df_nb_detections_PANGO_lineage_vs_seqrate)){
  df_nb_detections_PANGO_lineage_vs_seqrate$nb_WW_samples[i] <- length(unique(subset(df_variants,(!is.na(month_year))&(month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i]))$Sample))
  df_nb_detections_PANGO_lineage_vs_seqrate$nb_WW_detections[i] <- sum((!is.na(df_unique_detections_PANGO_lin$month_year))&(df_unique_detections_PANGO_lin$month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i]),na.rm=T)
  df_nb_detections_PANGO_lineage_vs_seqrate$nb_clinical_samples[i] <- length(unique(subset(df_clinical_detections,(!is.na(month_year))&(month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i]))$Sample.Name))
  df_nb_detections_PANGO_lineage_vs_seqrate$nb_clinical_detections[i] <- sum((!is.na(df_clinical_detections$month_year))&(df_clinical_detections$month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i]),na.rm=T)
  df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_richness[i] <- length(unique((subset(df_unique_detections_PANGO_lin,(!is.na(month_year))&(month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i])&(!is.na(PANGO_lineage))))$PANGO_lineage))
  df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_Shannon_entropy[i] <- ifelse(test = df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_unique_detections_PANGO_lin,(!is.na(month_year))&(month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i])))$PANGO_lineage),no=NA)
  df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_richness[i]>0,yes=df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_Shannon_entropy[i],no=NA))
  df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_Evenness[i] <- ifelse(test = df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_richness[i]>0,yes=df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_Shannon_entropy[i]/(log2(df_nb_detections_PANGO_lineage_vs_seqrate$WW_monthly_richness[i])),no=NA)
  df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_richness[i] <- length(unique((subset(df_clinical_detections,(!is.na(month_year))&(month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i])))$PANGOLIN))
  df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_Shannon_entropy[i] <- ifelse(test = df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_richness[i]>0,yes=get_entropy(target = (subset(df_clinical_detections,(!is.na(month_year))&(month_year==df_nb_detections_PANGO_lineage_vs_seqrate$month[i])))$PANGOLIN),no=NA)
  df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_antilog_Shannon_entropy[i] <- 2^(ifelse(test = df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_richness[i]>0,yes=df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_Shannon_entropy[i],no=NA))
  df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_Evenness[i] <- ifelse(test = df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_richness[i]>0,yes=df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_Shannon_entropy[i]/(log2(df_nb_detections_PANGO_lineage_vs_seqrate$Clinical_sampling_monthly_richness[i])),no=NA)
}
#Nb detections vs nb samples (Clinical sampling LSPQ)
ggplotRegression(fit = lmp(formula = nb_clinical_detections~nb_clinical_samples,data =df_nb_detections_PANGO_lineage_vs_seqrate,center = FALSE,perm = 9999),ggsave_path = output_workspace,the_filename = paste0("Nb_clinical_PANGO_lineages_detections_per_month_vs_sequencing_rate.png"),xlabl ="Sequencing rate\n(Number of sequenced samples per month)",ylabl = "Number of occurences of PANGO lineages\ndetected by clinical sampling")
#Nb detections vs nb samples (Wastewater)
ggplotRegression(fit = lmp(formula = nb_WW_detections~nb_WW_samples,data =df_nb_detections_PANGO_lineage_vs_seqrate,center = FALSE,perm = 9999),ggsave_path = output_workspace,the_filename = paste0("Nb_WW_PANGO_lineages_detections_per_month_vs_sequencing_rate.png"),xlabl ="Sequencing rate\n(Number of sequenced samples per month)",ylabl = "Number of occurences of PANGO lineages\ndetected in Wastewater")

df_nb_detections_PANGO_lineage_vs_seqrate$log10_p1_nb_WW_samples <- log10(df_nb_detections_PANGO_lineage_vs_seqrate$nb_WW_samples+1)
df_nb_detections_PANGO_lineage_vs_seqrate$log10_p1_nb_WW_detections <- log10(df_nb_detections_PANGO_lineage_vs_seqrate$nb_WW_detections+1)
df_nb_detections_PANGO_lineage_vs_seqrate$log10_p1_nb_clinical_samples <- log10(df_nb_detections_PANGO_lineage_vs_seqrate$nb_clinical_samples+1)
df_nb_detections_PANGO_lineage_vs_seqrate$log10_p1_nb_clinical_detections <- log10(df_nb_detections_PANGO_lineage_vs_seqrate$nb_clinical_detections+1)
#Nb detections vs nb samples per month (Clinical sampling LSPQ vs Wastewater)
formula <- y ~ poly(x,1,raw=T)
WW_vs_Clinical_samplings_Nb_detections_gg <- ggplot(df_nb_detections_PANGO_lineage_vs_seqrate) +
  geom_point(aes(x=log10_p1_nb_clinical_samples, y=log10_p1_nb_clinical_detections,col="Clinical")) +
  stat_smooth(aes(x=log10_p1_nb_clinical_samples, y=log10_p1_nb_clinical_detections),method = "lm", formula = formula,lty=2,col='blue') +
  stat_regline_equation(aes(x=log10_p1_nb_clinical_samples, y=log10_p1_nb_clinical_detections,label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 2, label.y = 1.5,col="blue") +
  stat_cor(aes(x=log10_p1_nb_clinical_samples, y=log10_p1_nb_clinical_detections,label =  ..p.label..),label.x = 2.5, label.y = 1.25,col="blue") +
  geom_point(aes(x=log10_p1_nb_WW_samples, y=log10_p1_nb_WW_detections,col="Wastewater")) +
  stat_smooth(aes(x=log10_p1_nb_WW_samples, y=log10_p1_nb_WW_detections),method = "lm", formula = formula,lty=2,col='red') +
  stat_regline_equation(aes(x=log10_p1_nb_WW_samples, y=log10_p1_nb_WW_detections,label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 1, label.y = 0,col="red") +
  stat_cor(aes(x=log10_p1_nb_WW_samples, y=log10_p1_nb_WW_detections,label =  ..p.label..),label.x = 1.5, label.y = -0.25,col="red") +
  theme_bw()+ ylab("Number of detections of\nPANGO lineages (log10 +1)")+ xlab("Sequencing rate\nlog10(Number of sequenced samples per month +1)") + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + scale_color_manual(values = c("Clinical"="blue","Wastewater"="red")) + labs(col="Sampling")
WW_vs_Clinical_samplings_Nb_detections_gg
ggsave(filename = paste0("Nb_WW_PANGO_lineages_detections_per_month_vs_sequencing_rate.png"), path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#WW_monthly_richness vs nb samples per month (Wastewater)
formula <- y ~ poly(x,1,raw=T)
WW_vs_Clinical_samplings_Richness_gg <- ggplot(df_nb_detections_PANGO_lineage_vs_seqrate) +
  geom_point(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_richness+1),col="Clinical")) +
  stat_smooth(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_richness+1)),method = "lm", formula = formula,lty=2,col='blue') +
  stat_regline_equation(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_richness+1),label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 2, label.y = 0.95,col="blue") +
  stat_cor(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_richness+1),label =  ..p.label..),label.x = 2.5, label.y = 0.75,col="blue") +
  geom_point(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_richness+1),col="Wastewater")) +
  stat_smooth(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_richness+1)),method = "lm", formula = formula,lty=2,col='red') +
  stat_regline_equation(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_richness+1),label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 0, label.y = 2,col="red") +
  stat_cor(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_richness+1),label =  ..p.label..),label.x = 0.5, label.y = 1.8,col="red") +
  theme_bw()+ ylab("Richness log10(number of unique lineages\ndetected per month+1)")+ xlab("Sequencing rate\nlog10(Number of sequenced samples per month+1)") + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + scale_color_manual(values = c("Clinical"="blue","Wastewater"="red")) + labs(col="Sampling")
WW_vs_Clinical_samplings_Richness_gg
ggsave(filename = paste0("WW_Monthly_richness_vs_sequencing_rate.png"), path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#Shannon entropy vs log10(nb samples per month+1) (Wastewater)
formula <- y ~ poly(x,1,raw=T)
WW_vs_Clinical_samplings_Shannon_gg <- ggplot(df_nb_detections_PANGO_lineage_vs_seqrate) +
  geom_point(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_Shannon_entropy+1),col="Clinical")) +
  stat_smooth(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_Shannon_entropy+1)),method = "lm", formula = formula,lty=2,col='blue') +
  stat_regline_equation(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_Shannon_entropy+1),label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 1.9, label.y = 0.25,col="blue") +
  stat_cor(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_Shannon_entropy+1),label =  ..p.label..),label.x = 2.4, label.y = 0.15,col="blue") +
  geom_point(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_Shannon_entropy+1),col="Wastewater")) +
  stat_smooth(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_Shannon_entropy+1)),method = "lm", formula = formula,lty=2,col='red') +
  stat_regline_equation(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_Shannon_entropy+1),label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 0, label.y = 1,col="red") +
  stat_cor(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_Shannon_entropy+1),label =  ..p.label..),label.x = 0.5, label.y = 0.9,col="red") +
  theme_bw()+ ylab("Monthly Shannon entropy (log10(y+1))")+ xlab("Sequencing rate\nlog10(Number of sequenced samples per month+1)") + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + scale_color_manual(values = c("Clinical"="blue","Wastewater"="red")) + labs(col="Sampling")
WW_vs_Clinical_samplings_Shannon_gg
ggsave(filename = paste0("WW_Monthly_Shannon_entropy_vs_sequencing_rate.png"), path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#WW_monthly_Evenness vs log10(nb samples per month+1) (Wastewater)
formula <- y ~ poly(x,1,raw=T)
WW_vs_Clinical_samplings_Evenness_gg <- ggplot(df_nb_detections_PANGO_lineage_vs_seqrate) +
  geom_point(aes(x=log10_p1_nb_clinical_samples, y=Clinical_sampling_monthly_Evenness,col="Clinical")) +
  stat_smooth(aes(x=log10_p1_nb_clinical_samples, y=Clinical_sampling_monthly_Evenness),method = "lm", formula = formula,lty=2,col='blue') +
  stat_regline_equation(aes(x=log10_p1_nb_clinical_samples, y=Clinical_sampling_monthly_Evenness,label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 1.8, label.y = 0.2,col="blue") +
  stat_cor(aes(x=log10_p1_nb_clinical_samples, y=Clinical_sampling_monthly_Evenness,label =  ..p.label..),label.x = 2.5, label.y = 0.15,col="blue") +
  geom_point(aes(x=log10_p1_nb_WW_samples, y=WW_monthly_Evenness,col="Wastewater")) +
  stat_smooth(aes(x=log10_p1_nb_WW_samples, y=WW_monthly_Evenness),method = "lm", formula = formula,lty=2,col='red') +
  stat_regline_equation(aes(x=log10_p1_nb_WW_samples, y=WW_monthly_Evenness,label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 1, label.y = 1.05,col="red") +
  stat_cor(aes(x=log10_p1_nb_WW_samples, y=WW_monthly_Evenness,label =  ..p.label..),label.x = 1.5, label.y = 1,col="red") +
  theme_bw()+ ylab("Monthly Evenness\n(Entropy/Max Entropy)")+ xlab("Sequencing rate\nlog10(Number of sequenced samples per month+1)") + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + scale_y_continuous(limits = c(0,1.2),breaks = seq(0,1,0.1)) + scale_color_manual(values = c("Clinical"="blue","Wastewater"="red")) + labs(col="Sampling")
WW_vs_Clinical_samplings_Evenness_gg
ggsave(filename = paste0("WW_Monthly_Evenness_vs_sequencing_rate.png"), path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

# #antilog2 Shannon entropy vs log10(nb samples per month+1) (Wastewater)
# formula <- y ~ poly(x,1,raw=T)
# WW_vs_Clinical_samplings_antilog_Shannon_gg <- ggplot(df_nb_detections_PANGO_lineage_vs_seqrate) +
#   geom_point(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_antilog_Shannon_entropy+1),col="Clinical")) +
#   stat_smooth(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_antilog_Shannon_entropy+1)),method = "lm", formula = formula,lty=2,col='blue') +
#   stat_regline_equation(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_antilog_Shannon_entropy+1),label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 0.2, label.y = 15,col="blue") +
#   stat_cor(aes(x=log10_p1_nb_clinical_samples, y=log10(Clinical_sampling_monthly_antilog_Shannon_entropy+1),label =  ..p.label..),label.x = 0.41, label.y = 13,col="blue") +
#   geom_point(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_antilog_Shannon_entropy+1),col="Wastewater")) +
#   stat_smooth(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_antilog_Shannon_entropy+1)),method = "lm", formula = formula,lty=2,col='red') +
#   stat_regline_equation(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_antilog_Shannon_entropy+1),label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),formula = formula,label.x = 0.2, label.y = 15,col="red") +
#   stat_cor(aes(x=log10_p1_nb_WW_samples, y=log10(WW_monthly_antilog_Shannon_entropy+1),label =  ..p.label..),label.x = 0.41, label.y = 13,col="red") +
#   theme_bw()+ ylab("Monthly antilog of Shannon entropy (log10(y+1))")+ xlab("Sequencing rate\nlog10(Number of sequenced samples per month+1)") + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + scale_color_manual(values = c("Clinical"="blue","Wastewater"="red")) + labs(col="Sampling")
# WW_vs_Clinical_samplings_antilog_Shannon_gg
# ggsave(filename = paste0("WW_Monthly_antilog_Shannon_entropy_vs_sequencing_rate.png"), path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#Effect of sequencing rate on lineage_detections 
png(filename = paste0(output_workspace,"Effect_of_sequencing_rate_on_lineage_detections.png"),width = 40,height=30, units = "cm",res = 1200)
grid.arrange(WW_vs_Clinical_samplings_Nb_detections_gg,WW_vs_Clinical_samplings_Richness_gg,WW_vs_Clinical_samplings_Shannon_gg,WW_vs_Clinical_samplings_Evenness_gg,nrow=2,ncol=2)
dev.off()

#Difference in Sequencing rate distributions
df_seqrate_clinical_only <- df_nb_detections_PANGO_lineage_vs_seqrate[,c("month","nb_clinical_samples")]
df_seqrate_clinical_only$Sampling <- "Clinical"
names(df_seqrate_clinical_only)[names(df_seqrate_clinical_only)=="nb_clinical_samples"] <- "metric"
df_seqrate_WW_only <- df_nb_detections_PANGO_lineage_vs_seqrate[,c("month","nb_WW_samples")]
df_seqrate_WW_only$Sampling <- "Wastewater"
names(df_seqrate_WW_only)[names(df_seqrate_WW_only)=="nb_WW_samples"] <- "metric"

ggplot(data = rbind(df_seqrate_clinical_only,df_seqrate_WW_only),aes(x=factor(as.character(Sampling),levels=c("Clinical","Wastewater")),y = log10(metric + 1),fill=factor(as.character(Sampling),levels=c("Clinical","Wastewater")))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("") + ylab("log10(Sequencing rate + 1)\n(Rate = Number of sequenced samples per month)") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) + scale_fill_manual(values = c("Clinical"="Blue","Wastewater"="red2")) + stat_compare_means(method = "wilcox") +scale_y_continuous(limits = c(0,ceiling(max(log10(rbind(df_seqrate_clinical_only,df_seqrate_WW_only)$metric+1),na.rm=T))),breaks=0:ceiling(max(log10(rbind(df_seqrate_clinical_only,df_seqrate_WW_only)$metric+1),na.rm=T))) + labs(fill="Sampling type")
ggsave(filename = "Sequencing_rates_distr_WW_vs_Clinical_sampling.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

#Difference in the number of detections distributions
df_nb_detections_PANGO_lineage_clinical_only <- df_nb_detections_PANGO_lineage_vs_seqrate[,c("month","nb_clinical_detections")]
df_nb_detections_PANGO_lineage_clinical_only$Sampling <- "Clinical"
names(df_nb_detections_PANGO_lineage_clinical_only)[names(df_nb_detections_PANGO_lineage_clinical_only)=="nb_clinical_detections"] <- "metric"
df_nb_detections_PANGO_lineage_WW_only <- df_nb_detections_PANGO_lineage_vs_seqrate[,c("month","nb_WW_detections")]
df_nb_detections_PANGO_lineage_WW_only$Sampling <- "Wastewater"
names(df_nb_detections_PANGO_lineage_WW_only)[names(df_nb_detections_PANGO_lineage_WW_only)=="nb_WW_detections"] <- "metric"

ggplot(data = rbind(df_nb_detections_PANGO_lineage_clinical_only,df_nb_detections_PANGO_lineage_WW_only),aes(x=factor(as.character(Sampling),levels=c("Clinical","Wastewater")),y = log10(metric + 1),fill=factor(as.character(Sampling),levels=c("Clinical","Wastewater")))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_jitter() + xlab("") + ylab("log10(Number of inferred occurences of\n a PANGO lineage per month+ 1)") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) + scale_fill_manual(values = c("Clinical"="Blue","Wastewater"="red2")) + stat_compare_means(method = "wilcox") +scale_y_continuous(limits = c(0,ceiling(max(log10(rbind(df_nb_detections_PANGO_lineage_clinical_only,df_nb_detections_PANGO_lineage_WW_only)$metric+1),na.rm=T))),breaks=0:ceiling(max(log10(rbind(df_nb_detections_PANGO_lineage_clinical_only,df_nb_detections_PANGO_lineage_WW_only)$metric+1),na.rm=T))) + labs(fill="Sampling type")
ggsave(filename = "Nb_monthly_detections_PANGO_lineage_distr_WW_vs_Clinical_sampling.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

#function to get mutation in genomic format
find_candidate_genomic_mutation_causing_aa_change_in_orf <- function(current_mutation_name){
  the_orf <- strsplit(x = current_mutation_name,split = ":")[[1]][1]
  the_mut <-strsplit(x = current_mutation_name,split = ":")[[1]][2]
  old_aa <- substr(the_mut,1,1)
  new_aa <- substr(the_mut,nchar(the_mut),nchar(the_mut))
  pos_in_orf_prot_seq <- as.integer(substr(the_mut,2,nchar(the_mut)-1))
  pos_start_codon_in_orf <- ((pos_in_orf_prot_seq-1)*3)+1
  pos_middle_codon_in_orf <- pos_start_codon_in_orf + 1
  pos_stop_codon_in_orf <- pos_start_codon_in_orf + 2
  pos_start_codon_in_genome <- v_start_orfs[the_orf] + pos_start_codon_in_orf - 1
  pos_middle_codon_in_genome <- v_start_orfs[the_orf] + pos_middle_codon_in_orf - 1
  pos_stop_codon_in_genome <- v_start_orfs[the_orf] + pos_stop_codon_in_orf - 1
  ref_codon <- substr(genome_refseq,pos_start_codon_in_genome,pos_stop_codon_in_genome)
  if (old_aa!=translate_seq(the_codon = ref_codon)){
    stop("Logical error: extracted reference codon does not translate into old aa!")
  }
  v_candidate_genomic_mutations <- NULL
  for (pos_in_codon in 1:3){
    for (current_new_nucl in setdiff(c("A","T","C","G"),substr(ref_codon,pos_in_codon,pos_in_codon))){
      if (pos_in_codon ==1){
        new_codon <- paste0(current_new_nucl,substr(ref_codon,2,3))
      }else if (pos_in_codon ==2){
        new_codon <- paste0(substr(ref_codon,1,1),current_new_nucl,substr(ref_codon,3,3))
      }else if (pos_in_codon ==3){
        new_codon <- paste0(substr(ref_codon,1,1),substr(ref_codon,2,2),current_new_nucl)
      }else{
        stop("Logical error: position in codon cannot be different than 1, 2 or 3!")
      }
      if (translate_seq(new_codon)==new_aa){
        v_candidate_genomic_mutations <- c(v_candidate_genomic_mutations,paste0(substr(ref_codon,pos_in_codon,pos_in_codon),pos_start_codon_in_genome+pos_in_codon-1,current_new_nucl))
      }
    }
  }
  return(v_candidate_genomic_mutations)
}

#Controls analysis 
#import and concatenate samples variants calling file
df_variants_Controls <- read.csv2(file = paste0(output_workspace,"variants_",lst_samples[1],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
if (nrow(df_variants_Controls)!=0){
  df_variants_Controls$Sample <- lst_samples[1]
}
for (i in 2:nb_samples){
  df_to_add_to_variants_df <- read.csv2(file = paste0(output_workspace,"variants_",lst_samples[i],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  if (nrow(df_to_add_to_variants_df)==0){
    next()
  }
  df_to_add_to_variants_df$Sample <- lst_samples[i]
  df_variants_Controls <- rbind(df_variants_Controls,df_to_add_to_variants_df)
  print(i)
}
#Do not confuse nucleotide T with the alias for the boolean value TRUE
df_variants_Controls$Ref <- gsub(pattern = "TRUE", replacement = "T",x = df_variants_Controls$Ref,fixed=TRUE) 
df_variants_Controls$VarAllele <- gsub(pattern = "TRUE", replacement = "T",x = df_variants_Controls$VarAllele,fixed=TRUE) 
df_variants_Controls$VarAllele <- gsub(pattern = "+", replacement = "",x = df_variants_Controls$VarAllele,fixed=TRUE) 
df_variants_Controls$VarAllele <- gsub(pattern = "-", replacement = "",x = df_variants_Controls$VarAllele,fixed=TRUE)
df_variants_Controls <- subset(df_variants_Controls,nchar(VarAllele)==1)

df_variants_Controls <- subset(df_variants_Controls,subset = !duplicated(paste0(Position,VarAllele,Sample,sep="")))
#make sure that the variant frequency filter is applied (freq variant >=0.05 on at least one strand)
df_variants_Controls$VarFreq <- (as.numeric(gsub(pattern = "%",replacement = "",x = df_variants_Controls$VarFreq,fixed = TRUE))/100)
df_variants_Controls <- subset(df_variants_Controls, (Reads1+Reads2>=100)&(VarFreq>=0.05))
#apply strand bias filter on variants (freq variant >=0.02 on each strand)
df_variants_Controls <- df_variants_Controls[(((df_variants_Controls$Reads2Plus/(df_variants_Controls$Reads1Plus+df_variants_Controls$Reads2Plus))>=0.02)&((df_variants_Controls$Reads2Minus/(df_variants_Controls$Reads1Minus+df_variants_Controls$Reads2Minus))>=0.02)),]
#Eliminate site with no variants
df_variants_Controls <- subset(df_variants_Controls,!is.na(VarAllele))
#sekect control samples 
df_variants_Controls <- subset(df_variants_Controls,unname(vapply(X = toupper(df_variants_Controls$Sample),FUN = function(x) grepl(pattern = "CTRL",x = x,fixed = T),FUN.VALUE = c(F))))
df_variants_Controls$ORF <- vapply(X = df_variants_Controls$Position,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_variants_Controls$gene <- vapply(X = df_variants_Controls$Position,FUN = function(x) return(find_gene_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
#duplicate variants of orf3a if they occur also in orf3b or orf3c
df_subset_orf3b_orf3c <- subset(df_variants_Controls,subset=(Position>=v_start_orfs["ORF3a"])&(Position<=v_end_orfs["ORF3a"]))
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"])),"ORF"] <- "ORF3b"
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"])),"ORF"] <- "ORF3c"
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"])),"gene"] <- "ORF3b"
df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"])),"gene"] <- "ORF3c"
df_subset_orf3b_orf3c <- subset(df_subset_orf3b_orf3c,ORF!="ORF3a")
df_variants_Controls <- rbind(df_variants_Controls,df_subset_orf3b_orf3c)
#duplicate variants of ORF7a if they occur also in ORF7b
df_subset_orf7b <- subset(df_variants_Controls,subset=(Position>=v_start_orfs["ORF7a"])&(Position<=v_end_orfs["ORF7a"]))
df_subset_orf7b[which((df_subset_orf7b$Position>=v_start_orfs["ORF7b"])&(df_subset_orf7b$Position<=v_end_orfs["ORF7b"])),"ORF"] <- "ORF7b"
df_subset_orf7b <- subset(df_subset_orf7b,ORF!="ORF7a")
df_variants_Controls <- rbind(df_variants_Controls,df_subset_orf7b)
#duplicate variants of N if they occur also in orf9c
df_subset_ORF9c <- subset(df_variants_Controls,subset=(Position>=v_start_orfs["N"])&(Position<=v_end_orfs["N"]))
df_subset_ORF9c[which((df_subset_ORF9c$Position>=v_start_orfs["ORF9c"])&(df_subset_ORF9c$Position<=v_end_orfs["ORF9c"])),"ORF"] <- "ORF9c"
df_subset_ORF9c <- subset(df_subset_ORF9c,ORF!="N")
df_variants_Controls <- rbind(df_variants_Controls,df_subset_ORF9c)
#Original codon and mutated codon
df_variants_Controls$ref_codon <- NA
df_variants_Controls$mut_codon <- NA
df_variants_Controls$pos_in_ORF <- NA
df_variants_Controls$pos_in_gene <- NA
df_variants_Controls$pos_in_protein <- NA
for (i in 1:nrow(df_variants_Controls)){
  df_variants_Controls$ref_codon[i] <-(get_ref_and_mutated_codon(the_position = df_variants_Controls$Position[i],ref_nucl = df_variants_Controls$Ref[i],new_nucl = df_variants_Controls$VarAllele[i]))$ref_codon
  df_variants_Controls$mut_codon[i] <-(get_ref_and_mutated_codon(the_position = df_variants_Controls$Position[i],ref_nucl = df_variants_Controls$Ref[i],new_nucl = df_variants_Controls$VarAllele[i]))$mutated_codon
  df_variants_Controls$old_aa[i] <- translate_seq(the_codon = df_variants_Controls$ref_codon[i] )
  df_variants_Controls$new_aa[i] <- translate_seq(the_codon = df_variants_Controls$mut_codon[i] )
  df_variants_Controls$pos_in_ORF[i] <- df_variants_Controls$Position[i] - v_start_orfs[df_variants_Controls$ORF[i]] + 1
  df_variants_Controls$pos_in_gene[i] <- df_variants_Controls$Position[i] - v_start_genes[df_variants_Controls$gene[i]] + 1
  df_variants_Controls$pos_in_protein[i] <- ceiling(df_variants_Controls$pos_in_gene[i]/3)
  
}
df_variants_Controls$mutation_name <- paste0(paste0(df_variants_Controls$Ref,df_variants_Controls$Position,df_variants_Controls$VarAllele,""),";",paste0(df_variants_Controls$old_aa,df_variants_Controls$pos_in_protein,df_variants_Controls$new_aa),";",df_variants_Controls$ORF,";",df_variants_Controls$gene)
#Define Nonsense and non-coding mutations
df_variants_Controls$is_nonsense <- (df_variants_Controls$new_aa=="Stop")
df_variants_Controls$is_UTR <- (is.na(df_variants_Controls$new_aa))
lst_control_samples <- unique(df_variants_Controls$Sample)
v_nb_SNVs_per_Control_sample <- as.vector(table(df_variants_Controls$Sample)[lst_control_samples])
names(v_nb_SNVs_per_Control_sample) <- lst_control_samples
v_nb_SNVs_per_Control_sample <- ifelse(test = is.na(v_nb_SNVs_per_Control_sample),yes=0,no=v_nb_SNVs_per_Control_sample)

df_plot_nb_mutations_per_Ctrl_sample <- data.frame(Sample=names(v_nb_SNVs_per_Control_sample),x=v_nb_SNVs_per_Control_sample,Control_type=ifelse(test = unname(vapply(X = toupper(lst_control_samples),FUN = function(x) grepl(pattern = "POSCTRL",x = x,fixed = T),FUN.VALUE = c(F))),yes="Positive controls", no="Negative controls"),stringsAsFactors = F)
ggplot(data = df_plot_nb_mutations_per_Ctrl_sample,aes(x=reorder(Sample,-x),y=x)) + geom_col(fill="black") + geom_text(aes(label=x), position=position_dodge(width=0.9), vjust=-0.25) + xlab("Sample") + ylab("Number of mutations")  + theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.title = element_text(size=16),legend.text = element_text(size=12),axis.text.x = element_text(size=7,angle=60,hjust=1),text=element_text(size=16)) + facet_wrap(~Control_type,ncol=1,scales = "free_x") + scale_y_continuous(limits = c(0,max(df_plot_nb_mutations_per_Ctrl_sample$x,na.rm=T)+50))
ggsave(filename = "Nb_mutations_per_Control_sample.png", path=output_workspace, width = 30, height = 25, units = "cm")

#Estimate lineages within-sample frequency
#list of samples in which lineages are detected
lst_unique_samples_for_PANGOlin_detection <- sort(unique(df_detected_marker_mutations_in_ww_samples$Sample))
#define parallelization parameters
lst_splits <- split(1:length(lst_unique_samples_for_PANGOlin_detection), ceiling(seq_along(1:length(lst_unique_samples_for_PANGOlin_detection))/(length(lst_unique_samples_for_PANGOlin_detection)/nb_cores)))
#Lineages pre-selection based on number of signature mutations + X=prevalence
the_f_parallel_preselect_based_on_nb_sig_muts_and_X_is_prevalence <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  current_df_sample_PANGO_lineages_frequencies <- NULL
  #find best linear model for the sample with grid search
  iii <- 1
  for (z in the_vec){
    current_sample <- lst_unique_samples_for_PANGOlin_detection[z]
    is_the_model_mtx_singular <- F
    current_lst_candidate_lineages <- rownames(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)[rowSums(as.matrix(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[,subset(df_detected_marker_mutations_in_ww_samples,(Sample==current_sample)&(nb_signature_mutations>=3))$label_mut_in_marker_fmt]))>0]
    #Signature mutations VAF in sample
    subset_df_current_sample <- subset(df_variants,(Sample==current_sample)&(label_mut_in_marker_fmt%in%colnames(mtx_is_signature_mutation_lineage)))
    subset_df_current_sample <- unique(subset_df_current_sample[,c("Sample","label_mut_in_marker_fmt","VarFreq")])
    #Unsolvable case: 1 marker mutation and 1 lineage (no regression possible)
    if (length(unique(subset_df_current_sample$label_mut_in_marker_fmt))==1){
      current_df_sample_PANGO_lineages_frequencies <- rbind(current_df_sample_PANGO_lineages_frequencies,data.frame(Sample=current_sample,lineage=ifelse(length(current_lst_candidate_lineages)>0, paste0(current_lst_candidate_lineages,collapse="/"),NA),frequency=NA,p.value=NA,rsq_sample_model=NA,adj_rsq_sample_model=NA,stringsAsFactors = F))
      next()
    }
    #prevalence data for mutations to fit in the model
    if (length(current_lst_candidate_lineages)==0){
      print(paste0("Cannot detect candidate lineages for Sample \'",current_sample,"\'!"))
      next()
    }else if (length(current_lst_candidate_lineages)==1){
      X <- matrix(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages[current_lst_candidate_lineages,subset_df_current_sample$label_mut_in_marker_fmt],ncol=1)
      X <- ifelse(test=X<0.5,yes=0,no=X)#remove the contribution of a lineage to the mutation frequency if the mutation is not present in most of the lineage sequences
      X <- matrix(X,ncol=1)
      vec_bool_rownames <- rowSums(X)>0
      X <- X[rowSums(X)>0,] # remove mutations that cannot be explained by the presence of any candidate lineages
      X <- matrix(X,ncol=1)
      rownames(X) <- subset_df_current_sample$label_mut_in_marker_fmt[vec_bool_rownames]
    }else{
      X <- t(mtx_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages[current_lst_candidate_lineages,subset_df_current_sample$label_mut_in_marker_fmt])
      X <- ifelse(test=X<0.5,yes=0,no=X)#remove the contribution of a lineage to the mutation frequency if the mutation is not present in most of the lineage sequences
      X <- X[rowSums(X)>0,] # remove mutations that cannot be explained by the presence of any candidate lineages
    }
    Y <- unique(subset_df_current_sample[,c("label_mut_in_marker_fmt","VarFreq")])$VarFreq
    names(Y) <- subset_df_current_sample$label_mut_in_marker_fmt
    Y <- unname(Y[rownames(X)])
    #grid search definition for coefficients
    eval(parse(text=paste0("df_grid_search_combinations <- ",paste0("expand.grid(",paste0(rep("c(0.1,0.5,0.9)", ncol(X)),collapse = ","),")"))))
    #make sure grid search values respect initial constraints
    df_grid_search_combinations <- subset(df_grid_search_combinations,(rowSums(df_grid_search_combinations) <= 1))
    # MCMC grid search
    current_sample_best_model <- NULL #initialization
    current_sample_best_model_loglik <- -Inf #initialization
    for (i in 1:nrow(df_grid_search_combinations)){
      if (length(current_lst_candidate_lineages)==1){
        is_the_model_mtx_singular <- T
        eval(parse(text=paste0("current_fit <- ConsReg(formula = Y~0+X,family = \'gaussian\',constraints = \'(X) <= 1,(X) > 0\',optimizer = \'mcmc\',LOWER = 0, UPPER = 1,ini.pars.coef = c(",paste0(unname(unlist(df_grid_search_combinations[i,])),collapse=","),"),penalty = 1E3)")))
        tryCatch(expr = {summary(current_fit);is_the_model_mtx_singular<-F},error=function(e) print(e))
        if (is_the_model_mtx_singular){
          break()
        }
        current_sample_signif_model_likelihood <- unname(unlist(summary(current_fit)$metrics["LogLik"]))
        if (current_sample_signif_model_likelihood > current_sample_best_model_loglik){
          current_sample_best_model <- current_fit 
          current_sample_best_model_loglik <- current_sample_signif_model_likelihood
        }
      }else{
        is_the_model_mtx_singular <- T
        eval(parse(text= paste0("current_fit <- ConsReg(formula = Y~0+X,family = \'gaussian\',constraints = \'(",paste0("X",colnames(X),collapse = " + "),") <= 1,(",paste0("X",colnames(X),collapse = " + "),") > 0\',optimizer = \'mcmc\',LOWER = 0, UPPER = 1,ini.pars.coef = c(",paste0(unname(unlist(df_grid_search_combinations[i,])),collapse=","),"),penalty = 1E3)")))
        tryCatch(expr = {summary(current_fit);is_the_model_mtx_singular<-F},error=function(e) print(e))
        if (is_the_model_mtx_singular){
          break()
        }
        if ((all(summary(current_fit)$coefficients[,"p.value"]<0.05))&(all(!is.na(summary(current_fit)$coefficients[,"p.value"])))){
          current_sample_signif_model_likelihood <- unname(unlist(summary(current_fit)$metrics["LogLik"]))
          if (current_sample_signif_model_likelihood > current_sample_best_model_loglik){
            current_sample_best_model <- current_fit 
            current_sample_best_model_loglik <- current_sample_signif_model_likelihood
          }
        }
      }
    }
    
    if (is_the_model_mtx_singular){ 
      current_df_sample_PANGO_lineages_frequencies <- rbind(current_df_sample_PANGO_lineages_frequencies,data.frame(Sample=current_sample,lineage=ifelse(length(current_lst_candidate_lineages)>0, paste0(current_lst_candidate_lineages,collapse="/"),NA),frequency=NA,p.value=NA,rsq_sample_model=NA,adj_rsq_sample_model=NA,stringsAsFactors = F))
      next()
    }
    if (length(current_sample_best_model)>0){
      Y_predicted <- predict(current_sample_best_model, newdata = data.frame(X))
      
      #find SST and SSE
      sst <- sum((Y - mean(Y))^2)
      sse <- sum((Y_predicted - Y)^2)
      
      #find R-Squared
      rsq <- 1 - sse/sst
      adj_rsq <- 1-(((1-rsq)*(nrow(X)-1))/(nrow(X)-ncol(X)-1))
      current_df_sample_PANGO_lineages_frequencies <- rbind(current_df_sample_PANGO_lineages_frequencies,data.frame(Sample=current_sample,lineage=current_lst_candidate_lineages,frequency=unname(unlist(summary(current_sample_best_model)$coefficients[1:(length(current_lst_candidate_lineages)),"Estimate"])),p.value=unname(unlist(summary(current_sample_best_model)$coefficients[1:(length(current_lst_candidate_lineages)),"p.value"])),rsq_sample_model=rsq,adj_rsq_sample_model=adj_rsq,stringsAsFactors = F))
    }
    if (iii%%5==0){
      print(paste0("Core ",i_cl,": ",iii," samples analyzed out of ",length(the_vec),"!"))
    }
    iii <- iii + 1
  }
  
  return(current_df_sample_PANGO_lineages_frequencies)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_PANGO_lineages_Freq_Estimation_analysis_preselect_based_on_nb_sig_muts_and_X_is_prevalence.txt"))
registerDoParallel(cl)
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings","glmnet","FD","vegan","ConsReg"))  %dopar% the_f_parallel_preselect_based_on_nb_sig_muts_and_X_is_prevalence(i_cl)
stopCluster(cl)

#Result of Lineage Frequencies estimation analysis
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$date <- v_samples_date[df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$Sample]
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$location <- v_samples_location[df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$Sample]
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence,!is.na(date))
#Time series estimated frequencies of lineages of interest
df_filtered_sample_PANGO_lineages_of_interest_frequencies <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence,(lineage%in%v_lineages_of_interest)&(!is.na(p.value))&(p.value<0.05)&(frequency>=0)&(frequency<=1))
df_filtered_sample_PANGO_lineages_of_interest_frequencies$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_filtered_sample_PANGO_lineages_of_interest_frequencies$lineage]
#Time series lineages frequencies
ggplot(data = df_filtered_sample_PANGO_lineages_of_interest_frequencies,mapping=aes(x=date,y=frequency)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = "stack") +
  xlab("Sampling date") + ylab(paste0("Frequency (t-test p-value < 0.05)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_Frequency_per_location_per_location_preselect_based_on_nb_sig_muts_and_X_is_prevalence.png", path=output_workspace, width = 40, height = 20, units = "cm")
#Time series estimated frequencies of lineages of interest
df_filtered_sample_PANGO_lineages_frequencies <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence,(!is.na(p.value))&(p.value<0.05)&(frequency>=0)&(frequency<=1))
#Time series lineages frequencies
ggplot(data = df_filtered_sample_PANGO_lineages_frequencies,mapping=aes(x=date,y=frequency)) + geom_col(mapping=aes(fill=factor(lineage,levels=names(sort(table(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$lineage),decreasing = F)))),position = "stack") +
  xlab("Sampling date") + ylab(paste0("Frequency (t-test p-value < 0.05)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") 
ggsave(filename = "Time_series_all_PANGO_lineages_Frequency_per_location_per_location_preselect_based_on_nb_sig_muts_and_X_is_prevalence.png", path=output_workspace, width = 40, height = 20, units = "cm")
#Barplot Model results
df_sample_frequencies_results <- aggregate(x = df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$p.value,by=list(Sample=df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_prevalence$Sample),FUN=function(x) ifelse(test=all((x<0.05)&(!is.na(x))),yes="Significant",no=ifelse(test=all((x>=0.05)&(!is.na(x))),yes="Non-significant",no=ifelse(test=all(is.na(x)),yes="Unsolvable",no="Mix"))))
df_percent_model_results <- data.frame(model_result=c("Significant","Non-significant","Unsolvable"),frequency=c(nrow(subset(df_sample_frequencies_results,x=="Significant"))/nrow(df_sample_frequencies_results),nrow(subset(df_sample_frequencies_results,x=="Non-significant"))/nrow(df_sample_frequencies_results),nrow(subset(df_sample_frequencies_results,x=="Unsolvable"))/nrow(df_sample_frequencies_results)))
ggplot(data = df_percent_model_results,mapping=aes(x=factor(model_result,levels=c("Significant","Non-significant","Unsolvable")),y=frequency)) + geom_col(mapping=aes(fill=factor(model_result,levels=c("Significant","Non-significant","Unsolvable"))),position = "stack") +
  xlab("") + ylab("Frequency") + theme_bw() + theme(axis.title = element_text(size=12), axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.1)) + labs(fill="Sample lineage frequency\nmodel result") + scale_fill_manual(values=c("Significant"="blue","Non-significant"="black","Unsolvable"="red"))
ggsave(filename = "WW_samples_lineage_freq_model_result_per_location_preselect_based_on_nb_sig_muts_and_X_is_prevalence.png", path=output_workspace, width = 15, height = 15, units = "cm")

#Lineages pre-selection based on number of signature mutations + X is binary (signature mutation or not)
the_f_parallel_preselect_based_on_nb_sig_muts_and_X_is_binary <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  current_df_sample_PANGO_lineages_frequencies <- NULL
  #find best linear model for the sample with grid search
  iii <- 1
  for (z in the_vec){
    current_sample <- lst_unique_samples_for_PANGOlin_detection[z]
    is_the_model_mtx_singular <- F
    current_lst_candidate_lineages <- rownames(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages)[rowSums(as.matrix(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[,subset(df_detected_marker_mutations_in_ww_samples,(Sample==current_sample)&(nb_signature_mutations>=3))$label_mut_in_marker_fmt]))>0]
    #Signature mutations VAF in sample
    subset_df_current_sample <- subset(df_variants,(Sample==current_sample)&(label_mut_in_marker_fmt%in%colnames(mtx_is_signature_mutation_lineage)))
    subset_df_current_sample <- unique(subset_df_current_sample[,c("Sample","label_mut_in_marker_fmt","VarFreq")])
    #Unsolvable case: 1 marker mutation and 1 lineage (no regression possible)
    if (length(unique(subset_df_current_sample$label_mut_in_marker_fmt))==1){
      current_df_sample_PANGO_lineages_frequencies <- rbind(current_df_sample_PANGO_lineages_frequencies,data.frame(Sample=current_sample,lineage=ifelse(length(current_lst_candidate_lineages)>0, paste0(current_lst_candidate_lineages,collapse="/"),NA),frequency=NA,p.value=NA,rsq_sample_model=NA,adj_rsq_sample_model=NA,stringsAsFactors = F))
      next()
    }
    #prevalence data for mutations to fit in the model
    if (length(current_lst_candidate_lineages)==0){
      print(paste0("Cannot detect candidate lineages for Sample \'",current_sample,"\'!"))
      next()
    }else if (length(current_lst_candidate_lineages)==1){
      X <- matrix(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[current_lst_candidate_lineages,subset_df_current_sample$label_mut_in_marker_fmt],ncol=1)
      X <- ifelse(test=X<0.5,yes=0,no=X)#remove the contribution of a lineage to the mutation frequency if the mutation is not present in most of the lineage sequences
      X <- matrix(X,ncol=1)
      vec_bool_rownames <- rowSums(X)>0
      X <- X[rowSums(X)>0,] # remove mutations that cannot be explained by the presence of any candidate lineages
      X <- matrix(X,ncol=1)
      rownames(X) <- subset_df_current_sample$label_mut_in_marker_fmt[vec_bool_rownames]
    }else{
      X <- t(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[current_lst_candidate_lineages,subset_df_current_sample$label_mut_in_marker_fmt])
      X <- ifelse(test=X<0.5,yes=0,no=X)#remove the contribution of a lineage to the mutation frequency if the mutation is not present in most of the lineage sequences
      X <- X[rowSums(X)>0,] # remove mutations that cannot be explained by the presence of any candidate lineages
    }
    Y <- unique(subset_df_current_sample[,c("label_mut_in_marker_fmt","VarFreq")])$VarFreq
    names(Y) <- subset_df_current_sample$label_mut_in_marker_fmt
    Y <- unname(Y[rownames(X)])
    #grid search definition for coefficients
    eval(parse(text=paste0("df_grid_search_combinations <- ",paste0("expand.grid(",paste0(rep("c(0.1,0.5,0.9)", ncol(X)),collapse = ","),")"))))
    #make sure grid search values respect initial constraints
    df_grid_search_combinations <- subset(df_grid_search_combinations,(rowSums(df_grid_search_combinations) <= 1))
    # MCMC optimization with grid search over initial parameters
    current_sample_best_model <- NULL #initialization
    current_sample_best_model_loglik <- -Inf #initialization
    for (i in 1:nrow(df_grid_search_combinations)){
      if (length(current_lst_candidate_lineages)==1){
        is_the_model_mtx_singular <- T
        eval(parse(text=paste0("current_fit <- ConsReg(formula = Y~0+X,family = \'gaussian\',constraints = \'(X) <= 1,(X) > 0\',optimizer = \'mcmc\',LOWER = 0, UPPER = 1,ini.pars.coef = c(",paste0(unname(unlist(df_grid_search_combinations[i,])),collapse=","),"),penalty = 1E3)")))
        tryCatch(expr = {summary(current_fit);is_the_model_mtx_singular<-F},error=function(e) print(e))
        if (is_the_model_mtx_singular){
          break()
        }
        current_sample_signif_model_likelihood <- unname(unlist(summary(current_fit)$metrics["LogLik"]))
        if (current_sample_signif_model_likelihood > current_sample_best_model_loglik){
          current_sample_best_model <- current_fit 
          current_sample_best_model_loglik <- current_sample_signif_model_likelihood
        }
      }else{
        is_the_model_mtx_singular <- T
        eval(parse(text= paste0("current_fit <- ConsReg(formula = Y~0+X,family = \'gaussian\',constraints = \'(",paste0("X",colnames(X),collapse = " + "),") <= 1,(",paste0("X",colnames(X),collapse = " + "),") > 0\',optimizer = \'mcmc\',LOWER = 0, UPPER = 1,ini.pars.coef = c(",paste0(unname(unlist(df_grid_search_combinations[i,])),collapse=","),"),penalty = 1E3)")))
        tryCatch(expr = {summary(current_fit);is_the_model_mtx_singular<-F},error=function(e) print(e))
        if (is_the_model_mtx_singular){
          break()
        }
        if ((all(summary(current_fit)$coefficients[,"p.value"]<0.05))&(all(!is.na(summary(current_fit)$coefficients[,"p.value"])))){
          current_sample_signif_model_likelihood <- unname(unlist(summary(current_fit)$metrics["LogLik"]))
          if (current_sample_signif_model_likelihood > current_sample_best_model_loglik){
            current_sample_best_model <- current_fit 
            current_sample_best_model_loglik <- current_sample_signif_model_likelihood
          }
        }
      }
    }
    
    if (is_the_model_mtx_singular){ 
      current_df_sample_PANGO_lineages_frequencies <- rbind(current_df_sample_PANGO_lineages_frequencies,data.frame(Sample=current_sample,lineage=ifelse(length(current_lst_candidate_lineages)>0, paste0(current_lst_candidate_lineages,collapse="/"),NA),frequency=NA,p.value=NA,rsq_sample_model=NA,adj_rsq_sample_model=NA,stringsAsFactors = F))
      next()
    }
    if (length(current_sample_best_model)>0){
      Y_predicted <- predict(current_sample_best_model, newdata = data.frame(X))
      
      #find SST and SSE
      sst <- sum((Y - mean(Y))^2)
      sse <- sum((Y_predicted - Y)^2)
      
      #find R-Squared
      rsq <- 1 - sse/sst
      adj_rsq <- 1-(((1-rsq)*(nrow(X)-1))/(nrow(X)-ncol(X)-1))
      current_df_sample_PANGO_lineages_frequencies <- rbind(current_df_sample_PANGO_lineages_frequencies,data.frame(Sample=current_sample,lineage=current_lst_candidate_lineages,frequency=unname(unlist(summary(current_sample_best_model)$coefficients[1:(length(current_lst_candidate_lineages)),"Estimate"])),p.value=unname(unlist(summary(current_sample_best_model)$coefficients[1:(length(current_lst_candidate_lineages)),"p.value"])),rsq_sample_model=rsq,adj_rsq_sample_model=adj_rsq,stringsAsFactors = F))
    }
    if (iii%%5==0){
      print(paste0("Core ",i_cl,": ",iii," samples analyzed out of ",length(the_vec),"!"))
    }
    iii <- iii + 1
  }
  
  return(current_df_sample_PANGO_lineages_frequencies)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_PANGO_lineages_Freq_Estimation_analysis_preselect_based_on_nb_sig_muts_and_X_is_binary.txt"))
registerDoParallel(cl)
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings","glmnet","FD","vegan","ConsReg"))  %dopar% the_f_parallel_preselect_based_on_nb_sig_muts_and_X_is_binary(i_cl)
stopCluster(cl)

#Result of Lineage Frequencies estimation analysis
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$date <- v_samples_date[df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$Sample]
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$location <- v_samples_location[df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$Sample]
df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary,!is.na(date))
#Time series estimated frequencies of lineages of interest
df_filtered_sample_PANGO_lineages_of_interest_frequencies <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary,(lineage%in%v_lineages_of_interest)&(!is.na(p.value))&(p.value<0.05)&(frequency>=0)&(frequency<=1))
df_filtered_sample_PANGO_lineages_of_interest_frequencies$label_lineage <- v_lineages_of_interest_with_who_desgnation[df_filtered_sample_PANGO_lineages_of_interest_frequencies$lineage]
#Time series lineages frequencies
ggplot(data = df_filtered_sample_PANGO_lineages_of_interest_frequencies,mapping=aes(x=date,y=frequency)) + geom_col(mapping=aes(fill=factor(label_lineage,levels=v_lineages_of_interest_with_who_desgnation[names(sort(table(subset(df_detected_marker_mutations_in_ww_samples,PANGO_lineage%in%v_lineages_of_interest)$PANGO_lineage),decreasing = F))])),position = "stack") +
  xlab("Sampling date") + ylab(paste0("Frequency (t-test p-value < 0.05)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") + scale_fill_manual(values=palette_PANGO_lineages_of_interest_with_who_designation)
ggsave(filename = "Time_series_PANGO_lineage_of_interest_Frequency_per_location_per_location_preselect_based_on_nb_sig_muts_and_X_is_binary.png", path=output_workspace, width = 40, height = 20, units = "cm")
#Time series estimated frequencies of lineages of interest
df_filtered_sample_PANGO_lineages_frequencies <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary,(!is.na(p.value))&(p.value<0.05)&(frequency>=0)&(frequency<=1))
#Time series lineages frequencies
ggplot(data = df_filtered_sample_PANGO_lineages_frequencies,mapping=aes(x=date,y=frequency)) + geom_col(mapping=aes(fill=factor(lineage,levels=names(sort(table(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$lineage),decreasing = F)))),position = "stack") +
  xlab("Sampling date") + ylab(paste0("Frequency (t-test p-value < 0.05)")) + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=8, angle = 60,hjust=1))+ scale_y_continuous(limit=c(0,1),breaks = seq(0,1,0.1)) + facet_wrap(~location, ncol=1) + labs(fill="PANGO lineages") 
ggsave(filename = "Time_series_all_PANGO_lineages_Frequency_per_location_per_location_preselect_based_on_nb_sig_muts_and_X_is_binary.png", path=output_workspace, width = 40, height = 20, units = "cm")
#Barplot Model results
df_sample_frequencies_results <- aggregate(x = df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$p.value,by=list(Sample=df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$Sample),FUN=function(x) ifelse(test=all((x<0.05)&(!is.na(x))),yes="Significant",no=ifelse(test=all((x>=0.05)&(!is.na(x))),yes="Non-significant",no=ifelse(test=all(is.na(x)),yes="Unsolvable",no="Mix"))))
df_percent_model_results <- data.frame(model_result=c("Significant","Non-significant","Unsolvable"),frequency=c(nrow(subset(df_sample_frequencies_results,x=="Significant"))/nrow(df_sample_frequencies_results),nrow(subset(df_sample_frequencies_results,x=="Non-significant"))/nrow(df_sample_frequencies_results),nrow(subset(df_sample_frequencies_results,x=="Unsolvable"))/nrow(df_sample_frequencies_results)))
ggplot(data = df_percent_model_results,mapping=aes(x=factor(model_result,levels=c("Significant","Non-significant","Unsolvable")),y=frequency)) + geom_col(mapping=aes(fill=factor(model_result,levels=c("Significant","Non-significant","Unsolvable"))),position = "stack") +
  xlab("") + ylab("Frequency") + theme_bw() + theme(axis.title = element_text(size=12), axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.1)) + labs(fill="Sample lineage frequency\nmodel result") + scale_fill_manual(values=c("Significant"="blue","Non-significant"="black","Unsolvable"="red"))
ggsave(filename = "WW_samples_lineage_freq_model_result_preselect_based_on_nb_sig_muts_and_X_is_binary.png", path=output_workspace, width = 15, height = 15, units = "cm")

#Determine weights for agnostic lineage pre-selection 
df_residuals_lineage_freq_model_result_preselect_based_on_nb_sig_muts_and_X_is_binary <- NULL
for (current_sample in sort(unique(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary$Sample))){
  current_lst_lineages <- subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary,(Sample==current_sample)&(frequency>=0)&(frequency<=1))$lineage
  #Signature mutations VAF in sample
  subset_df_current_sample <- subset(df_variants,(Sample==current_sample)&(label_mut_in_marker_fmt%in%colnames(mtx_is_signature_mutation_lineage)))
  subset_df_current_sample <- unique(subset_df_current_sample[,c("Sample","label_mut_in_marker_fmt","VarFreq")])
  Y <- unique(subset_df_current_sample[,c("label_mut_in_marker_fmt","VarFreq")])$VarFreq
  names(Y) <- subset_df_current_sample$label_mut_in_marker_fmt
  mtx_freq_current_lineages <- matrix(subset(df_sample_PANGO_lineages_frequencies_preselect_based_on_nb_sig_muts_and_X_is_binary,(Sample==current_sample)&(frequency>=0)&(frequency<=1)&(!is.na(frequency)))$frequency,ncol=1)
  if (nrow(mtx_freq_current_lineages)==0){
    next()
  }
  Y_predicted <- as.vector(t(matrix(mtx_pres_abs_mut_of_interest_in_NCBI_WORLDWIDE_lineages[current_lst_lineages,names(Y)],nrow=length(current_lst_lineages)))%*%mtx_freq_current_lineages)
  current_muts_abs_residuals <- abs(Y - Y_predicted)
  df_residuals_lineage_freq_model_result_preselect_based_on_nb_sig_muts_and_X_is_binary <- rbind(df_residuals_lineage_freq_model_result_preselect_based_on_nb_sig_muts_and_X_is_binary,data.frame(Sample=current_sample,label_mut_in_marker_fmt=names(current_muts_abs_residuals),residual=unname(current_muts_abs_residuals),model_mut_type=ifelse(test=names(current_muts_abs_residuals)%in%names(v_lineage_marker_mutations),yes="Marker mutation",no=ifelse(test=names(current_muts_abs_residuals)%in%colnames(mtx_is_signature_mutation_lineage),yes="Signature mutation (not marker)",no= "Other?")),stringsAsFactors = F))
}

#Compare distributions of residuals between signature only and marker mutations
ggplot(data = df_residuals_lineage_freq_model_result_preselect_based_on_nb_sig_muts_and_X_is_binary,aes(x=factor(as.character(model_mut_type),levels=c("Signature mutation (not marker)","Marker mutation")),y = residual,fill=factor(as.character(model_mut_type),levels=c("Signature mutation (not marker)","Marker mutation")))) + geom_violin() + geom_jitter() + geom_boxplot(width=0.075,fill="white") + xlab("") + ylab("Residuals (in absolute value)") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "right",axis.text.x=element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank()) + scale_fill_manual(values = c("Signature mutation (not marker)"="Blue","Marker mutation"="red2")) + stat_compare_means(method = "wilcox") + labs(fill="Model mutation type") + scale_y_continuous(limits=c(0,1.1))
ggsave(filename = "Residuals_Signature_olny_vs_Marker_mutations.png", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200)

#Table for data stratification analysis and sample metrics
df_sample_detection_metrics <- df_sample_stratifications_of_interest
df_sample_detection_metrics$avg_cov <- NA
df_sample_detection_metrics$nb_mutations <- NA
df_sample_detection_metrics$nb_detected_lineages <- NA
names(v_avg_depth_samples) <- lst_samples_original
for (i in 1:nrow(df_sample_detection_metrics)){
  df_sample_detection_metrics$avg_cov[i] <- v_avg_depth_samples[df_sample_detection_metrics$Sample[i]]
  df_sample_detection_metrics$nb_mutations[i] <- v_nb_SNVs_per_sample[df_sample_detection_metrics$Sample[i]]
  df_sample_detection_metrics$nb_detected_lineages[i] <- length(unique(subset(df_detected_marker_mutations_in_ww_samples,Sample==df_sample_detection_metrics$Sample[i])$PANGO_lineage))
}
#save tables
write.table(x=df_sample_detection_metrics,file = paste0(output_workspace,"Table_sample_metrics.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)
write.table(x=df_variants,file = paste0(output_workspace,"Table_sample_metrics.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)
write.table(x=subset(df_detected_marker_mutations_in_ww_samples$nb_signature_mutations>=3),file = paste0(output_workspace,"Table_detected_marker_mutations_in_ww_samples.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

library("session")
save.session(file = paste0(output_workspace,"ALL_ILLUMINA_",gsub(pattern = ":",replacement = "_",x = gsub(pattern = " ",replacement = "_",x = date())),"_RSession.Rda"))
