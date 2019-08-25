###########################################################################################
###########################################################################################
##  Script for RNA-Seq pipeline using Deseq2/fcros for RNASeq data of Aline Dubos        ##
##  Only DEA Deseq2 and fcros to be followed by pathways gage analysis. PLots QC apart 	 ##
###########################################################################################
###########################################################################################
##README
#Done by Mar Muniz. YHerault team. Last update 041018.
#Deseq2: Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and #
# dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 					      #
#Firsst part of the script read all samples together and normalized together. I already 
# probed during my previous trnascriptomic analysis that is better if the normalization is 
# split by sample type(miceline or cell type) and then combined together to see the common
#genes and pathayws disregulated on the DEA and DFA. And for PCA, after normalization
# splited by cell type we combine the results to create  a matrix and apply the prcomp
#  function to calculate the PCA
# other methods: http://core.sam.pitt.edu/node/7552
#info
#https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/rnaseq_diff_Snf2/rnaseq_diff_Snf2.html
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
#http://www.nathalievialaneix.eu/doc/pdf/tutorial-rnaseq.pdf
############################################################################################
############# Packages needed ###############


#in case of errors:
######################################
#unlink("/home/me/src/Rlibs/00LOCK-Rcpp", recursive = TRUE)
#install.packages("biomaRt", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#biocLite("biomaRt", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

######################################
#   "gplots",
#   repos = c("http://rstudio.org/_packages",
#   "http://cran.rstudio.com", dependencies = TRUE)
#)
source("http://bioconductor.org/biocLite.R")
library("PoiClaClu")
library("DESeq2")
library("apeglm")
library("genefilter")
library("GenomicRanges")
library("tidyr")
library("dplyr")
library("gtools")
library("gdata")
library("fcros")
library("biomaRt")
#Using Deseq2 docs
#############################################

#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22) # need to set it to get always the same random results and plots
setwd("/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_fetalFib/") 
#sessionInfo()
#wd:
#####################################################

#setwd("/shared/space2/herault/OmicsAnalyses/omicsAnalyses/ongoing/projects/mus_musculus/180618_BrV0_S16138_Dp1Yey/") 
wd <- getwd()

set.seed(22) # need to set it to get always the same random results and plots
setwd("/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_adultFib/") 
wd <- getwd()

setwd("/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_timeCourse/") 
wd <- getwd()


#save.image()
#load(".RData")
#set output directories
f_counts<- "/HTSeq/"
f_input <- "/input/"
f_results <- "/results/"
f_norm <- "normalization/"
f_qc <- "QC/"
f_plots <- "plots/"
f_DEGs <- "DEGs/"
f_w_low <- "with_Low_expressed/"
f_wo_low <- "wo_Low_expressed/"
f_DEGs <- "DEGs/"
f_EGs <- "EGs/"
f_corr <- "correlation/"
f_dist <- "SamplesDistance/"
f_heatmaps <- "heatmaps/"
f_dyrk  <- "Dyrk1a/"
f_data <- "data/"
f_edgeR <- "edgeR/"
f_Rdata <- "RData/"
f_deseq2 <-"deseq2/"

##########################
#folder for each analysis
##############################
f_name<- "monozygotic_twins_fibroblasts/"
folder <- f_name
name <- "monozygotic_twins_fibroblasts"

f_name<- "adult_fibroblasts/"
folder <- f_name
name <- "adult_fibroblasts"

f_name<- "fibroblast_time_course/"
folder <- f_name
name <- "fibroblast_time_course"


##########################


#creating the folders
#######################################################
#dir.create(file.path(getwd (), f_results), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_Rdata), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm, f_plots), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_edgeR), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_deseq2), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_deseq2,f_name), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_deseq2,f_name,f_Rdata), showWarnings = F)


#################################################################################################################
######################################## CODE ###################################################################
#functions:

mergeCounts = function(path){
filenames=mixedsort(list.files(path=path, full.names=TRUE))
#datalist = lapply(filenames, function(x){read.csv(file=x,header=T, stringsAsFactors=F,sep="\t",colClasses=c("character","numeric"))})
#Reduce(function(x,y) {left_join(x,y,by="ID")}, datalist)
datalist = lapply(filenames, function(x){read.table(file=x,header=T,sep="\t",stringsAsFactors=F)})
Reduce(function(x,y) {merge(x,y)}, datalist) #for whateever reason the option merge/reduce is not working ...
}


#################################################################################################################
######################################### Reading the files  ##############################################
###########################################################################################################

# A)From hisat2 HTSEQ counts run these lines to generate the final counts files to import in R as dataCounts
#IMPORTANT NOTE: run line by line! if not errors are incorpored in the counts! and mess with the reduce function!!!!!
# NOTE2: only change the name of the file VQBT15_S1.htseq.txt for one of your experiment to get the relationship 
# between ENSEMBL ID and gene_symbol
##the direct reading of the htseq counts with the dese2 specific function does not wwork because it 
# expect the second column to be the counts and not the gene symbols and because there are 5 more lines with 
#stats at the end of each count file generated by htseq count that should not be imported
#so after i forated the files in linux, now i will import all of them and combine them 

#######################################################################################################################
#######################################################################################################################
#cd /home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_adultFib
# cd /home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_fetalFib
#mkdir formated
#for i in *.txt;  do head -n -5 $i > formated/$i; done
#cd formated
#
#mkdir last_format

#for i in *.txt;  do awk '{OFS="\t"}{print $1,$3}'  $i > last_format/$i; done
#awk '{OFS="\t"} {print $1,$2} ' SRR1182244.htseq.txt > last_format/rownames_names.txt #fetal
#awk '{OFS="\t"} {print $1,$2} ' SRR1182262.htseq.txt > last_format/rownames_names.txt #adult
#
#cd last_format
##
#for file in *.txt; do echo -e ""ID"\t$file"$'\n'"$(cat -- "$file")" > "$file"; done
#nano rownames_names.txt   # name of second colummn to Gene_symbol

#B) importing into R
##############################

directory <-"/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_fetalFib/input/HTSeq/formated/last_format"
metafile='/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_fetalFib/input/metadata/metafile_letourneau2014fetalfib.xls'


directory <-"/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_adultFib/input/HTSeq/formated/last_format"
metafile='/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_adultFib/input/metadata/metafile_letourneau_adultFib.xls'

directory <-"/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_timeCourse/input/HTSeq/formated/last_format"
metafile='/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/letourneau_2014/letourneau2014_timeCourse/input/metadata/metafile_letourneau_timeCourse.xls'

sampleFiles <- grep("htseq",list.files(directory),value=TRUE)
#the htseq counts generated have 5 lines in the end with stats that we need to get rid of 
# so when trying to use it does not give an error because there are not 3 columns on these 
RawdataCounts <-mergeCounts(directory)




###### READING THE ANNOTATION FILE #########
annotation <- read.table(file=paste0(directory,"/rownames_names.txt" ), header=T,stringsAsFactors=F, sep="\t",quote='')

############################################
#######################################################################################
#
# WORKING WITH THE METAFILE: ADDING info of the  txt files names
#
#######################################################################################
library("xlsx")
#colData <- read.table(metafile, header=T,stringsAsFactors=F, sep="\t",quote='')
colData <- xlsx::read.xlsx(metafile,1, header=T,stringsAsFactors=F)
#colsConserv <- c("sampleName", "fileName" ,"condition", "sample","model","genotype","batch") #instead of batch that in this case is relevant we can use type or other optional
colData$sampleName <- mixedsort(sub("(.*)\\.htseq.*","\\1",sampleFiles))
colData$fileName <- paste0(colData$sampleName, ".htseq.txt")
colsConserv <- c("sampleName", "fileName" ,"condition", "sample","model","genotype","age","family") #instead of batch that in this case is relevant we can use type or other optional
#colsConserv <- c("sampleName", "fileName" ,"condition", "sample","model","genotype","batch") #instead of batch that in this case is relevant we can use type or other optional
colData <- colData[,colsConserv]
colData$fileName <- paste0(colData$sampleName, ".htseq.txt")
rownames(colData) <-colData[,4]
#colData$batch <- as.factor(colData$batch)

######################################### Reading the files  ##############################################
###########################################################################################################

rawCounts.ensemblNameAnnot.function <- function(RawCounts,metafile) {
	#This function produces the counts and metadata file for the DEA analysis
	#It also produces one of the input files needed to run gage data.all.deseq.fcros
	#If there is an error is because of save, just run twice the function to be able 
	#to save it like that
	countData <- as.matrix(RawCounts[,sampleFiles])
	rownames(countData) <- RawCounts[,1]#if i used the column with geneSymbol is ok but 92 duplicted names
	countData <- countData[,c(match(metafile[,2],colnames(countData)))]
	colnames(countData) <- metafile$sample
	assign("RawcountData",countData,.GlobalEnv)
	group <- as.factor(as.character(metafile$condition))
	assign("group",group,.GlobalEnv)
	save(RawcountData,colData,file=paste0(wd, f_input, "counts_and_metafile_use_as_input_",name,".RData")) #file input GAGE

}

#head(RawcountData)
rawCounts.ensemblNameAnnot.function(RawdataCounts,colData)
dim(RawcountData) #58825    GRCm38.p6.v93 54146, Rno.v96 32623, human GRCh38.p13.v96 58825


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
###  DEA using edgeR  ######
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

#set on each experiment
nConditions<-length(unique(colData[,3]))
nReplicates <-length(colData[,3])/nConditions
##DEA: Built of the design matrix
################################
design <- model.matrix(~0+group)
colnames(design) <- gsub("group","", colnames(design))



#
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
###  DEA using Deseq2  ######
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
#NOTE:
##DEseq analysis:
###################################################################################################
#DEseq2 analysis:
#not recomended prefiltering of low counts genes, by M love developer of deseq2
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

#vignette("DESeq2")
#"1.3.5 Pre-filtering
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, 
# there are two reasons which make pre-filtering useful: by removing rows in which there are no 
# reads or nearly no reads, we reduce the memory size of the dds data object and we increase the 
# speed of the transformation and testing functions within DESeq2. "

# You don't have to filter at all though. The safest threshold would be to not filter anything 
# above row sum of 1, and just let the data-driven software (which lives in the genefilter package,
# outside of DESeq2) choose the threshold that maximizes power. For more details you can read the 
# citation for the genefilter package which is also referenced in the DESeq2 paper section on 
# independent filtering.
###################################################################################################

##########################################
#for  adullt and fetal analysis split
##########################################

deseq2_DEA.function <- function(countData, metafile,organism,refName){
(countData, metafile,organism,refName){
#new norm, deseq2  lfcShrink fails if the names in the contrast have _ or capital letters so not use them!
#metafile$condition <- gsub(".*_","",metafile$condition)
#metafile$genotype <- tolower(gsub(".*_","",metafile$condition))
#metafile$samples <- tolower(gsub("_","",metafile$samples))
# 1.function Deseq: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
#dds = estimateSizeFactors(dds) #these 3 functions ar eintegrated inside DESeq function and are sequential
#dds = estimateDispersions(dds) #these 3 functions ar eintegrated inside DESeq function and are sequential
#dds = nbinomWaldTest(dds) #these 3 functions ar eintegarted inside DESeq function and are sequential. the estimated 
# standard error of a log2 fold change  is used to test if it is equal to zero, there is another way the likelihood ratio test (LRT)
# DO we have technical replicates? NO
# Not needed step for these data: we dont have technical replicates so is not neccesary to collapse them  collapseReplicates. Technical replicate implies multiple sequencing runs of the same library.
#
####################################
# To know the organism dataset
####################################
#ensembl=useMart("ensembl") #usemart gives probs becasue the servers are unstable! do not use!
#datasets <- listDatasets(ensembl)
#head(datasets)
#creating the annotation file:
ensembl_mo =useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = organism, mirror = "useast")
attributes_m = listAttributes(ensembl_mo)
inputNames <- rownames(countData)
#interesting attributes
attrib <-c(as.numeric(rownames(attributes_m[which(attributes_m[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributes_m[which(attributes_m[,1]=='external_gene_name'),][1,])))
attrib <-attributes_m[attrib,1] # "ensembl_gene_id"    "external_gene_name" "entrezgene" 
annotation = getBM(attributes=attrib,filters="ensembl_gene_id",values=inputNames, mart=ensembl_mo)
colnames(annotation) <- c("ID","Gene_symbol")
###############################

metafile$condition <- gsub(".*_","",metafile$condition)
reference  <- as.character(unique(metafile$condition[grep(refName, metafile$condition)]))
metafile$condition <- factor(metafile$condition)
metafile$condition <- relevel(metafile$condition, ref = reference)
colnames(countData) <- gsub("_",".",colnames(countData))
rownames(metafile) <- gsub("_",".",rownames(metafile))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metafile, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]  #minimal pre-filtering to remove rows that have only 0 or 1 read
dds <- DESeq(dds)
#resultsNames(dds) #to know the name of the contrasts

# 1.function Deseq: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
#dds = estimateSizeFactors(dds) #these 3 functions ar reintegrated inside DESeq function and are sequential
#dds = estimateDispersions(dds) #these 3 functions ar reintegrated inside DESeq function and are sequential
#dds = nbinomWaldTest(dds) #these 3 functions ar eintegrated inside DESeq function and are sequential. the estimated 
# standard error of a log2 fold change  is used to test if it is equal to zero, there is another way the likelihood ratio test (LRT)
# DO we have technical replicates? NO
# Not needed step for these data: we dont have technical replicates so is not neccesary to collapse them  collapseReplicates. 
# Technical replicate implies multiple sequencing runs of the same library.

res.deseq2 <- results(dds,alpha = 0.1,cooksCutoff=FALSE) 
resALL.deseq2 <- results(dds,alpha = 0.999,cooksCutoff=FALSE) 
res.deseq2$Gene_symbol <- annotation$Gene_symbol[match(rownames(res.deseq2), annotation$ID)]
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2],type="apeglm")
assign("res.deseq2",res.deseq2,.GlobalEnv)
resOrdered.deseq2 <- res.deseq2[order(res.deseq2$padj,decreasing = TRUE),]
tab005 <- table(resOrdered.deseq2$padj < 0.05)
tab01 <-table(resOrdered.deseq2$padj < 0.1)
tab05 <-table(resOrdered.deseq2$padj < 0.5)
resSig.005.deseq2 <- subset(resOrdered.deseq2, padj < 0.05)
resSig.01.deseq2 <- subset(resOrdered.deseq2, padj < 0.1)
resSig.05.deseq2 <- subset(resOrdered.deseq2, padj < 0.5)

#joint pipeline with gage: using the DESeq2 norm counts for gage
deseq2.fc <- resOrdered.deseq2$log2FoldChange
names(deseq2.fc) <- rownames(resOrdered.deseq2)
exp.fc <- deseq2.fc
##################################################################

data.DEA.results<- data.frame(data=name)
data.DEA.results$`DEA_Deseq2_0.05` <- as.numeric(tab005[['TRUE']])
data.DEA.results$`Genes_Deseq2_0.05` <-paste(as.character(rownames(res.deseq2[which(res.deseq2$padj < 0.05),])),collapse=" ")
data.DEA.results$`GeneSymbol_Deseq2_0.05` <-paste(as.character(res.deseq2[which(res.deseq2$padj < 0.05),'Gene_symbol']),collapse=" ")
data.DEA.results$`Number_upregulated_Deseq2_0.05` <-length(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange>0),7])
data.DEA.results$`Upregulated_Deseq2_0.05` <-paste(as.character(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange>0),7]),collapse=" ")
data.DEA.results$`Number_downregulated_Deseq2_0.05` <-length(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange<0),7])
data.DEA.results$`Downregulated_Deseq2_0.05` <-paste(as.character(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange<0),7]),collapse=" ")

data.DEA.results$`DEA_Deseq2_0.1` <- as.numeric(tab01[['TRUE']])
data.DEA.results$`Genes_Deseq2_0.1` <-paste(as.character(rownames(res.deseq2[which(res.deseq2$padj < 0.1),])),collapse=" ")
data.DEA.results$`GeneSymbol_Deseq2_0.1` <-paste(as.character(res.deseq2[which(res.deseq2$padj < 0.1),'Gene_symbol']),collapse=" ")
data.DEA.results$`Number_upregulated_Deseq2_0.1` <-length(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange>0),7])
data.DEA.results$`Upregulated_Deseq2_0.1` <-paste(as.character(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange>0),7]),collapse=" ")
data.DEA.results$`Number_downregulated_Deseq2_0.1` <-length(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange<0),7])
data.DEA.results$`Downregulated_Deseq2_0.1` <-paste(as.character(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange<0),7]),collapse=" ")

data.DEA.results$`DEA_Deseq2_0.5` <- as.numeric(tab05[['TRUE']])
data.DEA.results$`Genes_Deseq2_0.5` <-paste(as.character(rownames(res.deseq2[which(res.deseq2$padj < 0.5),])),collapse=" ")
data.DEA.results$`GeneSymbol_Deseq2_0.5` <-paste(as.character(res.deseq2[which(res.deseq2$padj < 0.5),'Gene_symbol']),collapse=" ")
data.DEA.results$`Number_upregulated_Deseq2_0.5` <-length(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange>0),7])
data.DEA.results$`Upregulated_Deseq2_0.5` <-paste(as.character(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange>0),7]),collapse=" ")
data.DEA.results$`Number_downregulated_Deseq2_0.5` <-length(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange<0),7])
data.DEA.results$`Downregulated_Deseq2_0.5` <-paste(as.character(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange<0),7]),collapse=" ")

resTable005 <-  res.deseq2[which(rownames(res.deseq2) %in% rownames(res.deseq2[which(res.deseq2$padj < 0.05),])),]
resTable01 <-  res.deseq2[which(rownames(res.deseq2) %in% rownames(res.deseq2[which(res.deseq2$padj < 0.1),])),]
resTable05 <-  res.deseq2[which(rownames(res.deseq2) %in% rownames(res.deseq2[which(res.deseq2$padj < 0.5),])),]
write.table(resTable005, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_005_",name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
write.table(resTable01, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_01_",name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
write.table(resTable05, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_05_",name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
write.table(res.deseq2, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_ALL_",name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
#New lines
write.table(data.DEA.results, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "data.DEA.results_", name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
save(resALL.deseq2,exp.fc,res.deseq2,data.DEA.results,resTable005,resTable01,resTable05,file=paste0(wd,f_results, f_norm,f_deseq2,f_name,f_Rdata,name,"_tmp_Input_for_GAGE.RData"))

assign("resLFC",resLFC,.GlobalEnv)
assign("resTable005",resTable005,.GlobalEnv)
assign("resTable05",resTable05,.GlobalEnv)
assign("resTable01",resTable01,.GlobalEnv)
assign("dds",dds,.GlobalEnv)
assign("data.DEA.results",data.DEA.results ,.GlobalEnv)
assign("exp.fc",exp.fc,.GlobalEnv)
assign("resALL.deseq2",resALL.deseq2,.GlobalEnv)#needed for volcano and correlation plot is all EGs 
};

#organism<-"mmusculus_gene_ensembl"
#organism<-"rnorvegicus_gene_ensembl"
organism<-"hsapiens_gene_ensembl"
refName <- "Euploid"
deseq2_DEA.function(RawcountData,colData,organism,refName)



head(data.DEA.results[,c(1,2,5,7,9,12,14,16,19,21)]) 


###################################################################################################

###################################################################################################

##########################################
#for  adullt and fetal analysis normalised
#  together as a time series
##########################################



deseq2_DEA.function <- function(countData, metafile,organism,level1,level2,level3,level4,deseq2_contrast2wt,deseq2_contrast3wt,deseq2_contrast4wt,annotation){
	#new norm, deseq2  lfcShrink fails if the names in the contrast have _ or capital letters so not use them!
	#metafile$condition <- gsub(".*_","",metafile$condition)
	#metafile$genotype <- tolower(gsub(".*_","",metafile$condition))
	#metafile$samples <- tolower(gsub("_","",metafile$samples))
	# 1.function Deseq: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
	#dds = estimateSizeFactors(dds) #these 3 functions ar eintegrated inside DESeq function and are sequential
	#dds = estimateDispersions(dds) #these 3 functions ar eintegrated inside DESeq function and are sequential
	#dds = nbinomWaldTest(dds) #these 3 functions ar eintegarted inside DESeq function and are sequential. the estimated 
	# standard error of a log2 fold change  is used to test if it is equal to zero, there is another way the likelihood ratio test (LRT)
	# DO we have technical replicates? NO
	# Not needed step for these data: we dont have technical replicates so is not neccesary to collapse them  collapseReplicates.
	# Technical replicate implies multiple sequencing runs of the same library.
	#
	####################################
	# To know the organism dataset
	####################################
	#ensembl=useMart("ensembl") #usemart gives probs becasue the servers are unstable! do not use!
	#datasets <- listDatasets(ensembl)
	#head(datasets)
	#creating the annotation file:


	#ensembl_mo =useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = organism, mirror = "useast") 
	##change to "uswest" if the ENSEMBL server is down and u get a query error	
	#attributes_m = listAttributes(ensembl_mo)
	#inputNames <- rownames(countData)
	#interesting attributes
	#attrib <-c(as.numeric(rownames(attributes_m[which(attributes_m[,1]=='ensembl_gene_id'),][1,])),\
		#as.numeric(rownames(attributes_m[which(attributes_m[,1]=='external_gene_name'),][1,])))
	#attrib <-attributes_m[attrib,1] # "ensembl_gene_id"    "external_gene_name" "entrezgene" 
	#annotation = getBM(attributes=attrib,filters="ensembl_gene_id",values=inputNames, mart=ensembl_mo);
	#colnames(annotation) <- c("ID","Gene_symbol")
	###############################

	metafile$conditionAge <- factor(metafile$conditionAge)
	metafile$conditionAge <- factor(metafile$conditionAge, levels=c(level1,level2,level3,level4)) #"conditionAge_mut (level2) vs wt (level1)" 
	colnames(countData) <- gsub("_",".",colnames(countData))
	rownames(metafile) <- gsub("_",".",rownames(metafile))

	dds <- DESeqDataSetFromMatrix(countData = countData, colData = metafile, design = ~ conditionAge)
	dds <- dds[ rowSums(counts(dds)) > 1, ]  #minimal pre-filtering to remove rows that have only 0 or 1 read
	dds <- DESeq(dds)
	print(resultsNames(dds)) #to know the name of the contrasts

	# 1.function Deseq: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
	#	dds = estimateSizeFactors(dds) #these 3 functions ar reintegrated inside DESeq function and are sequential
	#	dds = estimateDispersions(dds) #these 3 functions ar reintegrated inside DESeq function and are sequential
	#	dds = nbinomWaldTest(dds) #these 3 functions ar eintegrated inside DESeq function and are sequential. the estimated 
	# standard error of a log2 fold change  is used to test if it is equal to zero, there is another way the likelihood ratio test (LRT)
	# DO we have technical replicates? NO
	# Not needed step for these data: we dont have technical replicates so is not neccesary to collapse them  collapseReplicates. 
	# Technical replicate implies multiple sequencing runs of the same library.
			
	dds.one <- dds
	res.deseq2.one <- results(dds.one ,name=deseq2_contrast2wt, alpha = 0.1,cooksCutoff=FALSE) 
	res.deseq2.one$Gene_symbol <- annotation$Gene_symbol[match(rownames(res.deseq2.one), annotation$ID)]

	resOrdered.deseq2.one <- res.deseq2.one[order(res.deseq2.one$padj,decreasing = TRUE),]
	resALL.deseq2.one <- results(dds.one,name=deseq2_contrast2wt, alpha = 0.999,cooksCutoff=FALSE) 
	resLFC.one <- lfcShrink(dds.one, coef=resultsNames(dds.one)[2],type="apeglm")
	
	
	dds.two <- dds
	res.deseq2.two <- results(dds.two,name=deseq2_contrast3wt,  alpha = 0.1,cooksCutoff=FALSE) 
	res.deseq2.two$Gene_symbol <- annotation$Gene_symbol[match(rownames(res.deseq2.two), annotation$ID)]
	resOrdered.deseq2.two <- res.deseq2.two[order(res.deseq2.two$padj,decreasing = TRUE),]
	resLFC.two <- lfcShrink(dds.two, coef=resultsNames(dds.two)[2],type="apeglm")
	resALL.deseq2.two <- results(dds.two,name=deseq2_contrast3wt, alpha = 0.999,cooksCutoff=FALSE) 

	metafile$conditionAge <- factor(metafile$conditionAge, levels=c(level2,level4,level1,level3)) #"conditionAge_mut (level2) vs wt (level1)" 
	dds <- DESeqDataSetFromMatrix(countData = countData, colData = metafile, design = ~ conditionAge)
	dds <- dds[ rowSums(counts(dds)) > 1, ]  #minimal pre-filtering to remove rows that have only 0 or 1 read
	dds <- DESeq(dds)
	print(resultsNames(dds)) #to know the name of the contrasts

	dds.three <- dds
	res.deseq2.three <- results(dds.three ,name=deseq2_contrast4wt, alpha = 0.1,cooksCutoff=FALSE) 
	res.deseq2.three$Gene_symbol <- annotation$Gene_symbol[match(rownames(res.deseq2.three), annotation$ID)]
	resOrdered.deseq2.three <- res.deseq2.three[order(res.deseq2.three$padj,decreasing = TRUE),]
	resALL.deseq2.three <- results(dds.three,name=deseq2_contrast4wt, alpha = 0.999,cooksCutoff=FALSE) 
	resLFC.three <- lfcShrink(dds.three, coef=resultsNames(dds.three)[2],type="apeglm")
	
	assign("dds",dds,.GlobalEnv)
	assign("res.deseq2.adult",res.deseq2.one ,.GlobalEnv)
	assign("res.deseq2.euploid_fetal_vs_adult",res.deseq2.two ,.GlobalEnv)
	assign("res.deseq2.fetal",res.deseq2.three ,.GlobalEnv)
	assign("resLFC.adult",resLFC.one ,.GlobalEnv)
	assign("resLFC.euploid_fetal_vs_adult",resLFC.two ,.GlobalEnv)
	assign("resLFC.fetal",resLFC.three ,.GlobalEnv)
	assign("resALL.deseq2.adult",resALL.deseq2.one ,.GlobalEnv)#needed for volcano and correlation plots all EGs
	assign("resALL.deseq2.euploid_fetal_vs_adult",resALL.deseq2.two ,.GlobalEnv)#needed for volcano and correlation plots all EGs
	assign("resALL.deseq2.fetal",resALL.deseq2.three ,.GlobalEnv)#needed for volcano and correlation plots all EGs
	res.list<- list(resOrdered.deseq2.adult=resOrdered.deseq2.one,resOrdered.deseq2.euploid_fetal_vs_adult=resOrdered.deseq2.two,resOrdered.deseq2.fetal=resOrdered.deseq2.three)
	assign("res.list",res.list ,.GlobalEnv)

	for (i in 1:length(res.list)){ #run one by one
		contrastName <- gsub("resOrdered.deseq2.","",names(res.list)[i])
		resOrdered.deseq2 <- res.list[[i]]
		resSig.005.deseq2 <- subset(resOrdered.deseq2, padj < 0.05)
		resSig.01.deseq2 <- subset(resOrdered.deseq2, padj < 0.1)
		resSig.05.deseq2 <- subset(resOrdered.deseq2, padj < 0.5)
	
		tab005 <- table(resOrdered.deseq2$padj < 0.05)
		tab01 <-table(resOrdered.deseq2$padj < 0.1)
		tab05 <-table(resOrdered.deseq2$padj < 0.5)

		#joint pipeline with gage: using the DESeq2 norm counts for gage
		deseq2.fc <- resOrdered.deseq2$log2FoldChange
		names(deseq2.fc) <- rownames(resOrdered.deseq2)
		exp.fc <- deseq2.fc
		##################################################################

		data.DEA.results<- data.frame(data=gsub("resOrdered.deseq2.","",names(res.list)[i]))
		data.DEA.results$`DEA_Deseq2_0.05` <- as.numeric(tab005[['TRUE']])
		data.DEA.results$`Genes_Deseq2_0.05` <-paste(as.character(rownames(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.05),])),collapse=" ")
		data.DEA.results$`GeneSymbol_Deseq2_0.05` <-paste(as.character(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.05),'Gene_symbol']),collapse=" ")
		data.DEA.results$`Number_upregulated_Deseq2_0.05` <-length(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange>0),7])
		data.DEA.results$`Upregulated_Deseq2_0.05` <-paste(as.character(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange>0),7]),collapse=" ")
		data.DEA.results$`Number_downregulated_Deseq2_0.05` <-length(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange<0),7])
		data.DEA.results$`Downregulated_Deseq2_0.05` <-paste(as.character(resSig.005.deseq2[which(resSig.005.deseq2$log2FoldChange<0),7]),collapse=" ")

		data.DEA.results$`DEA_Deseq2_0.1` <- as.numeric(tab01[['TRUE']])
		data.DEA.results$`Genes_Deseq2_0.1` <-paste(as.character(rownames(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.1),])),collapse=" ")
		data.DEA.results$`GeneSymbol_Deseq2_0.1` <-paste(as.character(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.1),'Gene_symbol']),collapse=" ")
		data.DEA.results$`Number_upregulated_Deseq2_0.1` <-length(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange>0),7])
		data.DEA.results$`Upregulated_Deseq2_0.1` <-paste(as.character(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange>0),7]),collapse=" ")
		data.DEA.results$`Number_downregulated_Deseq2_0.1` <-length(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange<0),7])
		data.DEA.results$`Downregulated_Deseq2_0.1` <-paste(as.character(resSig.01.deseq2[which(resSig.01.deseq2$log2FoldChange<0),7]),collapse=" ")
		
		data.DEA.results$`DEA_Deseq2_0.5` <- as.numeric(tab05[['TRUE']])
		data.DEA.results$`Genes_Deseq2_0.5` <-paste(as.character(rownames(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.5),])),collapse=" ")
		data.DEA.results$`GeneSymbol_Deseq2_0.5` <-paste(as.character(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.5),'Gene_symbol']),collapse=" ")
		data.DEA.results$`Number_upregulated_Deseq2_0.5` <-length(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange>0),7])
		data.DEA.results$`Upregulated_Deseq2_0.5` <-paste(as.character(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange>0),7]),collapse=" ")
		data.DEA.results$`Number_downregulated_Deseq2_0.5` <-length(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange<0),7])
		data.DEA.results$`Downregulated_Deseq2_0.5` <-paste(as.character(resSig.05.deseq2[which(resSig.05.deseq2$log2FoldChange<0),7]),collapse=" ")
		
		resTable005 <-  resOrdered.deseq2[which(rownames(resOrdered.deseq2) %in% rownames(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.05),])),]
		resTable01 <-  resOrdered.deseq2[which(rownames(resOrdered.deseq2) %in% rownames(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.1),])),]
		resTable05 <-  resOrdered.deseq2[which(rownames(resOrdered.deseq2) %in% rownames(resOrdered.deseq2[which(resOrdered.deseq2$padj < 0.5),])),]
		
		write.table(resTable005, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_005_",contrastName,"_", name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
		write.table(resTable01, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_01_",contrastName,"_", name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
		write.table(resTable05, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "DEA_deseq2_05_",contrastName,"_", name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
		write.table(data.DEA.results, file=paste0(wd,f_results, f_norm,f_deseq2,f_name, "data.DEA.results_",contrastName,"_", name,".txt"), sep = "\t", col.names=TRUE, row.names=T,na = " ",)
		save(exp.fc,data.DEA.results,resTable005,resTable01,resTable05,file=paste0(wd,f_results, f_norm,f_deseq2,f_name,f_Rdata,name,"_",contrastName,"_tmp_Input_for_GAGE.RData")) #data.all.deseq.fcros.two,data.all.deseq.fcros.hippo,

		assign(paste0("resTable005.",contrastName),resTable005,.GlobalEnv)
		assign(paste0("resTable05.",contrastName),resTable05,.GlobalEnv)
		assign(paste0("resTable01.",contrastName),resTable01,.GlobalEnv)
		assign(paste0("data.DEA.results.",contrastName),data.DEA.results ,.GlobalEnv)
		assign(paste0("exp.fc.",contrastName),exp.fc,.GlobalEnv)
		i=i+1				
	};
						

};


organism<-"hsapiens_gene_ensembl"
#model all together
colData$conditionAge <- paste0(colData$genotype ,"_",colData$age)
colnames(colData)[9] <- "conditionAge"

level1 <- "Euploid_Adult"
level2 <- "Euploid_fetal"
level3 <- "T21_Adult"
level4 <- "T21_fetal"  
deseq2_contrast2wt <- "conditionAge_T21_Adult_vs_Euploid_Adult"
deseq2_contrast3wt <- "conditionAge_Euploid_fetal_vs_Euploid_Adult"
deseq2_contrast4wt <- "conditionAge_T21_fetal_vs_Euploid_fetal"
deseq2_DEA.function(RawcountData,colData,organism,level1,level2,level3,level4,deseq2_contrast2wt,deseq2_contrast3wt,deseq2_contrast4wt,annotation)
 
#see results  
data.DEA.results.adult[,c(1,2,5,7,9,12,14,16,19,21)]     
data.DEA.results.euploid_fetal_vs_adult[,c(1,2,5,7,9,12,14,16,19,21)]
data.DEA.results.fetal[,c(1,2,5,7,9,12,14,16,19,21)]  

#data.DEA.results.adult[,c(1,2,5,7,9,12,14,16,19,21)]     
#   data DEA_Deseq2_0.05 Number_upregulated_Deseq2_0.05
# adult           12673                           6875
#  Number_downregulated_Deseq2_0.05 DEA_Deseq2_0.1 Number_upregulated_Deseq2_0.1
#                             5798          14724                          8082
#  Number_downregulated_Deseq2_0.1 DEA_Deseq2_0.5 Number_upregulated_Deseq2_0.5
#                            6642          22850                         12450
#  Number_downregulated_Deseq2_0.5
#                           10400
# data.DEA.results.euploid_fetal_vs_adult[,c(1,2,5,7,9,12,14,16,19,21)]
#                    data DEA_Deseq2_0.05 Number_upregulated_Deseq2_0.05
# euploid_fetal_vs_adult            9353                           4663
#  Number_downregulated_Deseq2_0.05 DEA_Deseq2_0.1 Number_upregulated_Deseq2_0.1
#                             4690          11387                          5773
#  Number_downregulated_Deseq2_0.1 DEA_Deseq2_0.5 Number_upregulated_Deseq2_0.5
#                            5614          20108                         10203
#  Number_downregulated_Deseq2_0.5
#                            9905
# data.DEA.results.fetal[,c(1,2,5,7,9,12,14,16,19,21)]  
#  data DEA_Deseq2_0.05 Number_upregulated_Deseq2_0.05
# fetal               8                              8
# Number_downregulated_Deseq2_0.05 DEA_Deseq2_0.1 Number_upregulated_Deseq2_0.1
#                                0             17                            11
# Number_downregulated_Deseq2_0.1 DEA_Deseq2_0.5 Number_upregulated_Deseq2_0.5
#                               6             36                            22
#  Number_downregulated_Deseq2_0.5
#                              14



######################
#
# fcros DEGs analysis  ####
#
#######################

#saving the files in the following folders
##################################################
f_fcros<- "fcros/"
f_down <- "downregulated/"
f_up <- "upregulated/"
f_all <- "all/"


######
dir.create(file.path(getwd (), f_results, f_norm,f_fcros), showWarnings = FALSE) 
dir.create(file.path(getwd (), f_results, f_norm,f_fcros, f_name), showWarnings = FALSE) 
dir.create(file.path(getwd (), f_results, f_norm,f_fcros, f_name,f_down), showWarnings = FALSE) 
dir.create(file.path(getwd (), f_results, f_norm,f_fcros, f_name,f_all), showWarnings = FALSE) 
dir.create(file.path(getwd (), f_results, f_norm,f_fcros, f_name,f_up), showWarnings = FALSE) 
######
######

fcros_DEA.function<- function(dds,genoControlName,genoMutName, trim.opt,res.deseq2,annotation,nameContrast,name,nConditions){
# This function perform the DEA using fcros package 
# Doc: https://cran.r-project.org/web/packages/fcros/index.html
#    https://www.ncbi.nlm.nih.gov/pubmed/24423217
#  The input counts for fcros, should be the normalised counts, after adjusting by library size, possibly gene length,
#  batch correction if needed. 
#1. Using the normalised counts by Deseq2
#2. We identify the DEGs by fcros, the up and regulated ones, and merge into one df. 
#3. We also annotate these fcrosDEGs with the DESeq2 stats(FC, basemean, pval, fdr...)
#4. This function also produces one of the input files needed to run gage data.all.deseq.fcros
if (name == "Allsamples") {
data <- log2(data.frame(counts(dds, normalized=TRUE)+1))  #input data the normalised counts by deseq2 in log2
#ncolData <- length(colnames(data)[grep(nameContrast, colnames(data))])
cont <- colnames(data[,grep(genoControlName,colnames(data))])
test <- colnames(data[,grep(genoMutName,colnames(data))])
log2.opt <- 0 #data in the matrix "data" are expressed in a log2 scale
idnames <- rownames(data)
index <- 1:nrow(data)
data$EnsembleGeneName <- idnames
colnames(annotation) <- c("EnsembleGeneName","swissprot")
data <- left_join(data,annotation,by="EnsembleGeneName") 
data$index <- index
data <- data[, c("EnsembleGeneName",cont,test,"swissprot","index")]
af <- fcros(data, cont, test, log2.opt, trim.opt) #The first column of the matrix "xdata"  contain the gene IDs or their names.
af$fdr <- p.adjust(af$p.value, method="BH")
#Ex
#idnames       FC      FC2       ri  f.value  p.value   bounds   params     params_t      fdr     
 #  36627    36627    36627    36627    36627    36627        2        3            3    36627      
dataStats <- data.frame(minPval=min(af$p.value),maxPval=max(af$p.value),minFDR=min(af$fdr),maxFDR=max(af$fdr))
##############################################################
cuack <-af
cuack$bounds <- NULL
cuack$params <- NULL
cuack$params_t <- NULL
afDf <- do.call(cbind.data.frame, cuack)
colnames(afDf)[1] <- "EnsembleGeneName"
afDf$bounds <- paste0(unlist(af[[7]][1]),":",unlist(af[[7]][2]))
afDf$params <- paste0(unlist(af[[8]][1]),":",unlist(af[[8]][2]),":",unlist(af[[8]][3]))
afDf$params_t <- paste0(unlist(af[[9]][1]),":",unlist(af[[9]][2]),":",unlist(af[[9]][3]))
afDf <- left_join(afDf,annotation,by="EnsembleGeneName") 
afDf <- afDf[,c(1,11,2:10)]
totalETs <- length(unique(afDf[,1]))  #   #total expressed transcripts
totalEGs <- length(unique(afDf[,2]))  #   #total expressed genes
dupGenes <- unique(afDf[duplicated(afDf[,2]),2])
dataDup <- list()
for (i in 1:length(dupGenes)){
dataDup[[i]] <-afDf[which(afDf[,2]==dupGenes[i]),]
names(dataDup)[i] <- dupGenes[i]
i=i+1
}
#############################################################
#B. IDENTIFYING UP AND DOWN REGULATED DEGs
#############################################################
#############################################################
#ading the alpha values manually, as it cannot be computed by the fcrostopN functionwith so many DEGs i our list.
# So as this was the threshold defined by Dembele D. to be used equal to get a fdr <0,05

f.value <- afDf[,6]
id.up <- matrix(0, 1)
id.down <- matrix(0, 1)
alpha_up <- 0.975
alpha_dn <- 0.025

down.all <- 1;
up.all <- 1;
for (i in 1:totalETs) {
if (f.value[i] <= alpha_dn) { id.down[down.all] <- i; down.all <- down.all + 1; }
if (f.value[i] >= alpha_up) { id.up[up.all] <- i; up.all <- up.all + 1; }
}
data.down_all <- afDf[id.down[1:(down.all-1)], ];
ndown_all <- nrow(data.down_all);
data.up_all <- afDf[id.up[1:(up.all-1)], ];
nup_all <- nrow(data.up_all) 
fcrosDEGsStats <- data.frame(totalDEGs=ndown_all + nup_all,upregulatedDEGs=nup_all,downregulatedDEGs=ndown_all)
data.down_all$Regulation <- rep("DownRegulated", times=nrow(data.down_all))
data.up_all$Regulation <- rep("UpRegulated", times=nrow(data.up_all))
data.down_all$FoldChange <- as.numeric(2^data.down_all$FC2)*-1
data.up_all$FoldChange <- as.numeric(2^data.up_all$FC2)
options(scipen = 999) #turn off scientific notation

data.all <- rbind(data.up_all,data.down_all)
data.all <- data.all[,c(1:2,13,3:4,12,5:11)]
doubleRegGenes.tmp <- data.down_all[which(data.down_all[,2] %in% data.up_all[,2]),]
doubleRegGenes.tmp2 <- data.up_all[which(data.up_all[,2] %in% data.down_all[,2]),]
doubleRegGenes <- rbind(doubleRegGenes.tmp,doubleRegGenes.tmp2)
if (length(doubleRegGenes[,1]>0)){
data.all <- data.all[ !(data.all[,1] %in% doubleRegGenes[,1]), ]
double <-unique(doubleRegGenes[,2])
for (i in 1:length(double)){
doubleRegGenes[,2] <-paste0(doubleRegGenes[,2],"_",doubleRegGenes[,1])
i=i+1 
}
data.all <- rbind(data.all,doubleRegGenes)
}
res.deseq2Df <-data.frame(res.deseq2)
colnames(res.deseq2Df)<- paste0(colnames(res.deseq2Df),"_DESeq2")
colnames(res.deseq2Df)
res.deseq2Df$EnsembleGeneName <- rownames(res.deseq2Df) 
colnames(data.all)<- ifelse(colnames(data.all)=="EnsembleGeneName","EnsembleGeneName",ifelse(colnames(data.all)=="swissprot","swissprot",paste0(colnames(data.all),"_fcros")))
data.all.deseq.fcros <- left_join(data.all,res.deseq2Df,by="EnsembleGeneName")[,-20]
save(data.all.deseq.fcros,file=paste0(f_input,"tgdyrk_tmp_Input_for_GAGE.RData")) #file input GAGE

#############################################################
#creating output files
#############################################################
write.table(data.down_all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_down, "DEA_fcros_downregulated_",name,"_",nameContrast,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
write.table(data.up_all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_up, "DEA_fcros_upregulated_",name,"_",nameContrast,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
write.table(data.all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_all, "DEA_fcros_",name,"_",nameContrast,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
#write.table(data.all.deseq.fcros, file=paste0(wd, f_results, f_norm, "DEA_fcros_deseq2_comparison",name,"_",nameContrast,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
assign(paste0("data.up_all_",name,"_",nameContrast),data.up_all,.GlobalEnv)
assign(paste0("data.down_all_",name,"_",nameContrast),data.down_all,.GlobalEnv)
assign(paste0("res.deseq2Df_",name,"_",nameContrast),res.deseq2Df,.GlobalEnv)
assign(paste0("afDf_",name,"_",nameContrast),afDf,.GlobalEnv)
assign(paste0("dataStats_",name,"_",nameContrast),dataStats,.GlobalEnv)
assign(paste0("totalETs_",name,"_",nameContrast),totalETs,.GlobalEnv)
assign(paste0("totalEGs_",name,"_",nameContrast),totalEGs,.GlobalEnv)
assign(paste0("fcrosDEGsStats_",name,"_",nameContrast),fcrosDEGsStats,.GlobalEnv)
assign(paste0("doubleRegGenes_",name,"_",nameContrast),doubleRegGenes,.GlobalEnv)
assign(paste0("data.all.deseq.fcros_",name,"_",nameContrast),data.all.deseq.fcros,.GlobalEnv)


}else{

data <- log2(data.frame(counts(dds, normalized=TRUE)+1))  #input data the normalised counts by deseq2 in log2
ncolData <- ncol(data)
cont <- colnames(data[,grep(genoControlName,colnames(data))])
print(cont)
test <- colnames(data[,grep(genoMutName,colnames(data))])
print(test)
log2.opt <- 0 #data in the matrix "data" are expressed in a log2 scale
idnames <- rownames(data)
index <- 1:nrow(data)
data$EnsembleGeneName <- idnames
colnames(annotation) <- c("EnsembleGeneName","swissprot")
data <- left_join(data,annotation,by="EnsembleGeneName") 
data$index <- index
data <- data[, c("EnsembleGeneName",cont,test,"swissprot","index")]
af <- fcros(data, cont, test, log2.opt, trim.opt) #The first column of the matrix "xdata" should contain the gene IDs or their names.
af$fdr <- p.adjust(af$p.value, method="BH")
dataStats <- data.frame(minPval=min(af$p.value),maxPval=max(af$p.value),minFDR=min(af$fdr),maxFDR=max(af$fdr))
##############################################################
cuack <-af
cuack$bounds <- NULL
cuack$params <- NULL
cuack$params_t <- NULL
afDf <- do.call(cbind.data.frame, cuack)
colnames(afDf)[1] <- "EnsembleGeneName"
afDf$bounds <- paste0(unlist(af[[7]][1]),":",unlist(af[[7]][2]))
afDf$params <- paste0(unlist(af[[8]][1]),":",unlist(af[[8]][2]),":",unlist(af[[8]][3]))
afDf$params_t <- paste0(unlist(af[[9]][1]),":",unlist(af[[9]][2]),":",unlist(af[[9]][3]))
afDf <- left_join(afDf,annotation,by="EnsembleGeneName") 
afDf <- afDf[,c(1,11,2:10)]
totalETs <- length(unique(afDf[,1]))  #36627   #total expressed transcripts
totalEGs <- length(unique(afDf[,2]))  #36604   #total expressed genes
dupGenes <- unique(afDf[duplicated(afDf[,2]),2])
dataDup <- list()
for (i in 1:length(dupGenes)){
dataDup[[i]] <-afDf[which(afDf[,2]==dupGenes[i]),]
names(dataDup)[i] <- dupGenes[i]
i=i+1
}
#############################################################
#B. IDENTIFYING UP AND DOWN REGULATED DEGs
#############################################################
#############################################################
f.value <- afDf[,6]
id.up <- matrix(0, 1)
id.down <- matrix(0, 1)
alpha_up <- 0.975
alpha_dn <- 0.025

down.all <- 1;
up.all <- 1;
for (i in 1:totalETs) {
if (f.value[i] <= alpha_dn) { id.down[down.all] <- i; down.all <- down.all + 1; }
if (f.value[i] >= alpha_up) { id.up[up.all] <- i; up.all <- up.all + 1; }
}
data.down_all <- afDf[id.down[1:(down.all-1)], ];
ndown_all <- nrow(data.down_all);
data.up_all <- afDf[id.up[1:(up.all-1)], ];
nup_all <- nrow(data.up_all) 
fcrosDEGsStats <- data.frame(totalDEGs=ndown_all + nup_all,upregulatedDEGs=nup_all,downregulatedDEGs=ndown_all)
data.down_all$Regulation <- rep("DownRegulated", times=nrow(data.down_all))
data.up_all$Regulation <- rep("UpRegulated", times=nrow(data.up_all))
data.down_all$FoldChange <- as.numeric(2^data.down_all$FC2)*-1
data.up_all$FoldChange <- as.numeric(2^data.up_all$FC2)
options(scipen = 999) #turn off scientific notation

data.all <- rbind(data.up_all,data.down_all)
data.all <- data.all[,c(1:2,13,3:4,12,5:11)]
doubleRegGenes.tmp <- data.down_all[which(data.down_all[,2] %in% data.up_all[,2]),]
doubleRegGenes.tmp2 <- data.up_all[which(data.up_all[,2] %in% data.down_all[,2]),]
doubleRegGenes <- rbind(doubleRegGenes.tmp,doubleRegGenes.tmp2)
if (length(doubleRegGenes[,1]>0)){
data.all <- data.all[ !(data.all[,1] %in% doubleRegGenes[,1]), ]
double <-unique(doubleRegGenes[,2])
for (i in 1:length(double)){
doubleRegGenes[,2] <-paste0(doubleRegGenes[,2],"_",unique(doubleRegGenes[,1]))
i=i+1 
}
data.all <- rbind(data.all,doubleRegGenes)
}
res.deseq2Df <-data.frame(res.deseq2)
colnames(res.deseq2Df)<- paste0(colnames(res.deseq2Df),"_DESeq2")
colnames(res.deseq2Df)
res.deseq2Df$EnsembleGeneName <- rownames(res.deseq2Df) 
colnames(data.all)<- ifelse(colnames(data.all)=="EnsembleGeneName","EnsembleGeneName",ifelse(colnames(data.all)=="swissprot","swissprot",paste0(colnames(data.all),"_fcros")))
data.all.deseq.fcros <- left_join(data.all,res.deseq2Df,by="EnsembleGeneName")[,-20]

#############################################################
#creating output files
#############################################################
write.table(data.down_all,
 	file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_down, "DEA_fcros_downregulated_",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ")
write.table(data.up_all,
 	file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_up, "DEA_fcros_upregulated_",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ")
write.table(data.all, 
	file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_all, "DEA_fcros_",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ")

if(nConditions>2){
nameTested <- paste0(genoMutName,"_",genoControlName)
assign(paste0("data.up_all_",nameTested),data.up_all,.GlobalEnv)
assign(paste0("data.down_all_",nameTested),data.down_all,.GlobalEnv)
assign(paste0("res.deseq2_",nameTested),res.deseq2Df,.GlobalEnv)
assign(paste0("afDf_",nameTested),afDf,.GlobalEnv) #needed for volcano and correlation plot is all EGs 
assign(paste0("dataStats_",nameTested),dataStats,.GlobalEnv)
assign(paste0("totalETs_",nameTested),totalETs,.GlobalEnv)
assign(paste0("totalEGs_",nameTested),totalEGs,.GlobalEnv)
assign(paste0("fcrosDEGsStats_",nameTested),fcrosDEGsStats,.GlobalEnv)
assign(paste0("doubleRegGenes_",nameTested),doubleRegGenes,.GlobalEnv)
assign(paste0("data.all.deseq.fcros_",nameTested),data.all.deseq.fcros,.GlobalEnv)
save(afDf,data.all.deseq.fcros,file=paste0(wd, f_results,f_Rdata,name,"_",nameTested,"_tmp_Input_for_GAGE.RData")) #data.all.deseq.fcros, 

}else{
assign("data.up_all",data.up_all,.GlobalEnv)
assign("data.down_all",data.down_all,.GlobalEnv)
assign("res.deseq2",res.deseq2Df,.GlobalEnv)
assign("afDf",afDf,.GlobalEnv) #needed for volcano and correlation plot is all EGs 
assign("dataStats",dataStats,.GlobalEnv)
assign("totalETs",totalETs,.GlobalEnv)
assign("totalEGs",totalEGs,.GlobalEnv)
assign("fcrosDEGsStats",fcrosDEGsStats,.GlobalEnv)
assign("doubleRegGenes",doubleRegGenes,.GlobalEnv)
assign("data.all.deseq.fcros",data.all.deseq.fcros,.GlobalEnv)
save(afDf,data.all.deseq.fcros,file=paste0(wd, f_results,f_Rdata,name,"_tmp_Input_for_GAGE.RData")) #data.all.deseq.fcros, 

}

}
};

f_fcros<- "fcros/"
genoControlName <- "Euploid"
genoMutName <- "T21"


fcros_DEA.function(dds,genoControlName,genoMutName, 0.25,res.deseq2,annotation,name,name,nConditions) #for 2 genotypes wt and mut
head(data.DEA.results[,c(1,2,5,7,9,12,14,16,19,21)]) 
fcrosDEGsStats

#
#fetal_fibroblasts
#
#data.DEA.results[,c(1,2,5,7,9,12,14,16,19,21)]  
#                           data DEA_Deseq2_0.05 Number_upregulated_Deseq2_0.05
# monozygotic_twins_fibroblasts             174                            112
#  Number_downregulated_Deseq2_0.05 DEA_Deseq2_0.1 Number_upregulated_Deseq2_0.1
#                               62            266                           174
#  Number_downregulated_Deseq2_0.1 DEA_Deseq2_0.5 Number_upregulated_Deseq2_0.5
#                              92           1638                           866
#  Number_downregulated_Deseq2_0.5
#                             772
#
# fcrosDEGsStats
#  totalDEGs upregulatedDEGs downregulatedDEGs
#      1401             717               684


###############################
# save needed objects for 
# GAGE PATHWAYS ANALYSIS
###############################

save(RawcountData,colData,exp.fc,resTable005,resTable01,
	resTable05,data.DEA.results, dds,res.deseq2,resLFC,
	resALL.deseq2,afDf,fcrosDEGsStats,data.all.deseq.fcros,
	data.up_all,data.down_all,file=paste0(wd, f_results,f_Rdata,name,"_all_Input_for_GAGE.RData")) 

#SAVE SESSION
save.image(file=paste0(wd, f_results, f_Rdata, name,"_session.RData"))
#load(file=paste0(wd, f_results, f_Rdata, name,"_session.RData"))


#################################################################################################################
######################################## CODE ###################################################################


################################################################################################################################
################################################################################################################################
############################### END. Finished on the  19/02/2019 ####################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################
################################################################################################################################
################################################################################################################################


sessionInfo()
