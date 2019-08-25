###########################################################################################
###########################################################################################
##  Script for microarray Hgu133 DEA analysis, pipeline using Deseq2/fcros               ##
###########################################################################################
###########################################################################################
##README
#Done by Mar Muniz. YHerault team. Last update 041017.
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

source("http://bioconductor.org/biocLite.R")
#install.packages("collapseRow", type = "source")
#install.packages("collapseRow", dependencies=TRUE, repos='http://bioconductor.org/biocLite.R')
#options(download.file.method = "wget")
library("reshape2")
library("collapseRow")
library("gtools")
library("GEOquery")
library("ggfortify")
library("tidyr")
library("dplyr")
library("oligo")
library("ff")
library("nlme")
library("calibrate")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("genefilter")
library("scales")
library("limma")
library("GenomicRanges")
library("dendextend")
library("amap")
library("PoiClaClu")
library("apeglm")
library("gdata")
library("fcros")
library("biomaRt")
library("hgu133a.db")
library("pd.hg.u133a")
library("hgu133acdf")
library("hugene10stv1cdf")
library("hgu133aprobe")
library("affy")
library("openxlsx")
#####################################################

setwd("/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/onWorking/") 
wd <- getwd()


#set output directories
f_counts<- "/input/raw/"
f_results <- "/results/"
f_norm <- "normalization/"
f_qc <- "QC/"
f_plots <- "plots/"
f_DEGs <- "DEGs/"
f_DEGs <- "DEGs/"
f_EGs <- "EGs/"
f_corr <- "correlation/"
f_dist <- "SamplesDistance/"
f_heatmaps <- "heatmaps/"
f_data <- "data/"
f_edgeR <- "edgeR/"
f_Rdata <- "RData/"
f_deseq2 <-"deseq2/"
f_limma <-"limma/"
f_DEGs05 <- "fdr05/"
f_DEGs005 <- "fdr005/"
f_DEGs0005 <- "fdr0005/"
f_DEGs01 <- "fdr01/"
f_fcros <-"fcros/"
f_raw <- "raw/"
f_rma <- "rma/"


##########################
#folder for each analysis
##############################
f_name<- "human_fetal_brain/"
folder <- f_name
name <- "human_fetal_brain"

f_name<- "human_adult_PFC/"
folder <- f_name
name <- "human_adult_PFC"
##########################


#creating the folders
#######################################################
#dir.create(file.path(getwd (), f_results), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_raw), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_raw,f_name), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_rma), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_rma,f_name), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_Rdata), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm, f_plots), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_limma,f_name,f_Rdata), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_fcros), showWarnings = F)

dir.create(file.path(getwd (), f_results, f_norm,f_limma), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_limma,f_name), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_limma,f_name,f_DEGs01), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_limma,f_name,f_DEGs05), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_limma,f_name,f_DEGs005), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_norm,f_limma,f_name,f_DEGs0005), showWarnings = F)

###################################################################
#To change manually depending of the mice line we are analysing   #
###################################################################


#setings of plots colours:
#########################################
#col1 <- colorRampPalette(rev(brewer.pal(9, "RdPu")))(70)  # We choose the colours for the heatmap.
col1 <- colorRampPalette(rev(brewer.pal(9, "RdPu")))(12)  # We choose the colours for the heatmap.
col2 <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#load session
#load(file=paste0(wd,f_results, f_Rdata,name,"_session.RData"))

#################################################################################################################
######################################## CODE ###################################################################
#folllowing the guidelines to label the objects metadata and expression data from microarrays as described:
#https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
#inputfeatureData <- '/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/Mao_hgu133array/input/info/A-AFFY3-3adf.xlsx'
#inputphenoData<- '/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/Mao_hgu133array/input/info/E-GEOD-1397sdrf.xlsx'
#phenoData <- xlsx::read.xlsx(inputphenoData,1, header=T,stringsAsFactors=F)

#from processed data --> deprecated as we found the raw data
###################################################################
#inputphenoData<- '/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/Mao_hgu133array/input/info/E-GEOD-1397.sdrf.txt'
##SDRF <- read.delim(inputphenoData) #phenoData
#rownames(SDRF) <- SDRF$Array.Data.File
#SDRF <- AnnotatedDataFrame(SDRF)
#featureData <- xlsx::read.xlsx(inputfeatureData,1, header=T,stringsAsFactors=F)  #pData
#featureData <- unique(featureData[,c(1,7)])
#dim(featureData) #annotated probes with ENSEMBLID 22283 with the chip
#library("xlsx")
#exprs <- data.matrix(read.table(file=paste0(wd,f_counts,"E-GEOD-1397-processed-data-1626201579.txt" ), header=T, sep="\t",row.names=1,as.is=TRUE,quote=''))
#exprs <- exprs[,grep("\\.\\.1", colnames(exprs),perl=T)]


#################################################################################################################
######################################## CODE ###################################################################

inputTargets <- '/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/Mao_hgu133array/input/info/TargetsFile.xlsx'
colData <- xlsx::read.xlsx(inputTargets,1, header=T,stringsAsFactors=F) 
rownames(colData) <-colData[,'sample']
f_input <- "Mao_hgu133array/"

#is an old affymetrix array so dont import using oligo read.cellfiles it does not provide the annotation 
# of the probes as featureNames(affyRaw) dont work
celpath <- paste0("/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/",f_input,"input/raw/" )


cels = list.files(celpath, pattern = "cel",ignore.case =TRUE)
sapply(paste("data", cels, sep="/"), gunzip)
affyRaw <- affy::ReadAffy(celfile.path=celpath,filenames=cels,cdfname="hgu133acdf",verbose=TRUE) 
pData(affyRaw)
#exprs(affyRaw[c(1:3),])

data_info <- as.data.frame(pData(affyRaw))
data_info$fileName <- rownames(data_info)
data_full_info <- left_join(data_info,colData, by="fileName") # careful, the order of the rows have to be the index order!
pData(affyRaw) <-data_full_info

sampleNames(affyRaw) <- colnames(exprs(affyRaw[,]))
sampleNames(affyRaw) <- data_full_info[order(match(data_full_info$fileName,sampleNames(affyRaw))),5]
affyRaw <- affyRaw[,order(sampleNames(affyRaw))] #to order the files inside affyRaw in the order of their label name
Samples <- sampleNames(affyRaw) #storing the names in a new variable for later use
#pData(affyRaw)

#Retrieve the probes annotation as it was not imported in the .cel
feat = affyRaw@featureData
featureNames(affyRaw)
length(featureNames(affyRaw)) #nb of probes 22283
length(affy::probeNames(affyRaw)) #nb of probes in the array 247965 # average nb probes gene=11.12799

max(exprs(affyRaw)); #32617.5 fetal; adult: 46287 #number of rows and features. 
min(exprs(affyRaw)); #28fetal; adult: 29.3

#exprs(affyRaw);  #raw intensities the same as:
#oligo::intensity(affyRaw); #same as exprs(affyRaw)

  
###
################################################ NORMALIZATION ###########################################################
##########################################################################################################################
#rma method, proceeds with background subtraction, normalization and summarization using median-polish.
#using the processed data given by the authors in GEO ---> not to use, i will use the raw data
#exprs <- data.matrix(read.table(file=paste0("/home/munizmom/Documents/YH_Lab/DSmodels_RNASeq/projects/",f_input,"input/dataProcessed/E-GEOD-5390-processed-data-1756239242.txt" ), header=T, sep="\t",row.names=1,as.is=TRUE,quote=''))
#length(unique(rownames(exprs(exprs))))

#Normalization by RMA
###############################
affy_norm_t <- affy::rma(affyRaw, target = "core"); 
#length(unique(rownames(exprs(affy_norm_t))))
#dim(affy_norm_t)
write.csv(exprs(affyRaw), file=paste0(wd,f_results, f_raw,f_name, name, "_data_not_Normalized_probes.txt"))
write.exprs(affy_norm, file=paste0(wd, f_results, f_norm,f_rma,f_name, name, "_data_Normalized_probes.txt"))
write.exprs(affy_norm_t, file=paste0(wd,f_results, f_norm,f_rma,f_name, name, "_data_Normalized_transcript.txt"))
save(affyRaw,affy_norm,affy_norm_t, file=paste0(wd,f_results, f_Rdata, name, "_data_bforeNorm_Normalized_and_transcript.RData"))

f_rmaMatrix<- "/GEO_rma_matrix/"
dir.create(file.path(getwd (), f_rmaMatrix), showWarnings = FALSE) 
write.exprs(affy_norm_t, file=paste0(wd, f_rmaMatrix,name,"_subSeries_matrix.txt"))

f_metadata<- "/GEO_samples_metafile/"
dir.create(file.path(getwd (), f_metadata), showWarnings = FALSE) 
write.table(data_full_info,  file=paste0(wd, f_metadata,name,"_samples_metafile.txt"), col.names=T,row.names=F, sep="\t")


eset <- exprs(affy_norm_t)
summary(eset)
eset_df <- as.data.frame(eset)
colnames(eset_df) <- colnames(eset)
eset_df$PROBEID <- rownames(eset)


###
############################
#       ANNOTATION         #
############################
#pag9: man https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
#Symbols =data.frame(probes=probes,Gene_symbol=unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))) #this way of annotating gives less nb geneNames not to use! https://www.biostars.org/p/67294/
#Entrez_IDs = data.frame(probes=probes,ENTREZID=unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))) #this way of annotating gives less nb geneNames not to use!

probes=rownames(eset)
#columns(hgu133a.db) #to know the possible info to be annotated

annotation <- AnnotationDbi::select(hgu133a.db,
                                  keys = (probes),
                                  columns = c("SYMBOL", "ENTREZID","GENENAME","ENSEMBLTRANS"),
                                  keytype = "PROBEID")
annotation <- subset(annotation, !is.na(SYMBOL))  #13768
annotationUniq <- annotation[!duplicated(annotation$SYMBOL),] #13768

# filtering1: filtering out the probes that hybridate to different genes as they are not specific
annoGrouped <- group_by(annotation, PROBEID)
annoSummarized <- as.data.frame(dplyr::summarize(annoGrouped, no_of_matches = n_distinct(SYMBOL))) #none probe is mapping to diff genes, so is ok no need to filter by that
annoFiltered <- filter(annoSummarized, no_of_matches > 1) #only 1
probe_stats <- annoFiltered  
dim(affy_norm_t[-which(featureNames(affy_norm_t) %in% probe_stats$PROBEID),])  #21037
dim(affy_norm_t[which(featureNames(affy_norm_t) %in% probe_stats$PROBEID),])#1246

affy_norm_tF <- affy_norm_t[-which(featureNames(affy_norm_t) %in% probe_stats$PROBEID),]
########################################




###
############################
#           DEA            #
############################
#######Differential expression analysis.#######

###############################################################################################
############################ Methods to form  the design matrix.   ############################
###############################################################################################
# Method1: create a design matrix which includes a coefficient for the mutant vs wild type difference, 
##better to use METHOD1 WHEN is comparison between 2 groups, instead for the fetal data the method2

# method2 : create a design matrix which includes separate coefficients for wild type
#          and mutant mice andthen extract the difference as a contrast
#

#method2 : create a design matrix which includes separate coefficients for wild type
#          and mutant mice and then extract the difference as a contrast

###############################################################################################


#######################
#######################
#######################
# fetal: method2 as there is a multigroup comparison
#######################
#######################
#######################
colData$condition <- gsub("_", ".", colData$condition)  #reemplacing _ by . If not the contrast matrix does not work
colData$condition <- factor(colData$condition, levels =unique(colData$condition))
design <- model.matrix(~0+colData$condition + colData$bioSample)
colnames(design) <-gsub("colData\\$","",gsub("colData\\$condition","",colnames(design)))

#6. Contrast matrix.
contrast.matrix <- makeContrasts(Astrocyte=T21.Astrocyte-Euploid.Astrocyte, 
	Cerebellum=T21.Cerebellum-Euploid.Cerebellum,
	Cerebrum=T21.Cerebrum-Euploid.Cerebrum,
	Heart=T21.Heart-Euploid.Heart,
	levels=design)

fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
list_0 <- as.data.frame(topTable(fit2, number=Inf, adjust="BH"))
list_0$PROBEID <- rownames(list_0)

#
#LIST OF DIFFERENTIAL ESPRESSED GENES and Probes (DEGs):
#
##################################################################################

#NOTE: #the next three lines have the same results! I will use the topTableF because i prefer to not have results in scientific notation. 
#list_just_lines <- as.data.frame(topTableF(just_lines, number=Inf, adjust.method="BH")) 
## first column is log2FC for ech sample
#list_just_lines <- as.data.frame(topTable(just_lines, number=Inf, coef=NULL, adjust.method="BH"))
#list_just_lines <- as.data.frame(topTableF(just_lines, number=Inf))
#Also using the topTable function in all the ebayes data or just in a subset with the contrast of mic elines or whatever columns we want
# to subset does not change the results of TOPTABLE, makes sense, as the diff contrast are based in contrast already defined by mice line so the aveExpr, 
#F, padj and p val does not change after subseting if we chose all the mice lines,
#so i dont need to work with a subset of the ebayes list_0.



#FOR ALL DATA, lets do the annotation of all the genes in the chip
##########################################################################################
##########################################################################################
#Combine gene annotations with raw data
list_0 <- left_join(list_0,annotationUniq,by="PROBEID");
dim(list_0); #23151
list_0 <- subset(list_0, !is.na(SYMBOL));  #13768
dim(list_0);  #13768

list_0$GeneName_transcript <- paste0(list_0$SYMBOL, "_", list_0$ENSEMBLTRANS) 
list_0$GeneName_transcript <- gsub("_NA","", list_0$GeneName_transcript ) 

#stats for all EGs  eset/ebayes object annotated with Gene_symbol
########################################################################
min(ebayes$p.value) #3.399192e-09  #fetal
max(ebayes$p.value) #0.9998653
min(list_1b) #0.0003029768
min(list_0$adj.P.Val) #0.0006056694
min(ebayes$p.value) #3.399192e-09


#per model
##########################################
pdf(file=paste0(wd, f_results, f_norm,f_rma,f_name,name, "_histogram_DEGs_vs_EGs.pdf"))
op <- par(mar=c(7,10,10,7), mgp=c(63,1,0))  # to create the comparison between cell types 
for (i in 1:ncol(fit2$p.value)) 
{hist(fit2$p.value[,i], breaks=100, col="lightseagreen", xlab="p.value", ylab="number of genes", main=paste0("Proportion of DEGs vs EGs based on pval in: ", colnames(fit2$p.value)[i]))} #proportion of  DEGs vs non DEGs 
dev.off()


#############################
#remove dup gene names, do summary per gene
####In case of several probes to the same gene transcript, do not leverage take the one with stronger var between conditions as recommended
#This is ok https://www.biostars.org/p/51756/. https://www.researchgate.net/post/Which_probe_set_should_I_consider_for_each_gene_in_affymetrix_microarray
#We should not average them, just pick the one with more variance between the samples or leave all of them
repeated_limma <- unique(list_0[duplicated(list_0$GeneName_transcript),'GeneName_transcript']) #none. ook!
#############################
list_0 <-left_join(list_0, eset_df, by="PROBEID");

#fetal
list_0$Mean_expression_Euploid_Astrocyte=rowMeans(list_0[,grep("Euploid_Astrocyte",colnames(list_0))], na.rm=F)
list_0$Mean_expression_T21_Astrocyte=rowMeans(list_0[,grep("T21_Astrocyte",colnames(list_0))], na.rm=F)

list_0$Mean_expression_Euploid_Cerebellum=rowMeans(list_0[,grep("Euploid_Cerebellum",colnames(list_0))], na.rm=F)
list_0$Mean_expression_T21_Cerebellum=rowMeans(list_0[,grep("T21_Cerebellum",colnames(list_0))], na.rm=F)

list_0$Mean_expression_Euploid_Cerebrum=rowMeans(list_0[,grep("Euploid_Cerebrum",colnames(list_0))], na.rm=F)
list_0$Mean_expression_T21_Cerebrum=rowMeans(list_0[,grep("T21_Cerebrum",colnames(list_0))], na.rm=F)

list_0$Mean_expression_Euploid_Heart=rowMeans(list_0[,grep("Euploid_Heart",colnames(list_0))], na.rm=F)
list_0$Mean_expression_T21_Heart=rowMeans(list_0[,grep("T21_Heart",colnames(list_0))], na.rm=F)
list_0 <- list_0[,c(9:10,14,11:13,15:47,1:8)]
list_0$FC_Astrocyte <- 2^list_0$Astrocyte  #simplification of what limma does, limma_logFC=mean(log2(G1) -log2(Gr2))
list_0$FC_Cerebellum <- 2^list_0$Cerebellum  #simplification of what limma does, limma_logFC=mean(log2(G1) -log2(Gr2))
list_0$FC_Cerebrum<- 2^list_0$Cerebrum  #simplification of what limma does, limma_logFC=mean(log2(G1) -log2(Gr2))
list_0$FC_Heart <- 2^list_0$Heart  #simplification of what limma does, limma_logFC=mean(log2(G1) -log2(Gr2))
list_0 <- list_0[,c(1:39,48:51,40:47)]


# From the annotated EGS, with 
# ones have and fdr<0.5 or 
# fdr <0.05,or <0.1, or <0.005?
###############################
degs.fdr05<- subset(list_0, list_0$adj.P.Val < 0.5);
degs.fdr005<- subset(list_0, list_0$adj.P.Val < 0.05);
degs.fdr01<- subset(list_0, list_0$adj.P.Val < 0.1);
degs.fdr0005<- subset(list_0, list_0$adj.P.Val < 0.005);
#fetal:
list_0$Regulation_limma_Astrocyte <- ifelse(list_0$Astrocyte<0, "DownRegulated Astrocyte", "UpRegulated Astrocyte")
list_0$Regulation_limma_Cerebellum <- ifelse(list_0$Cerebellum<0, "DownRegulated Cerebellum", "UpRegulated Cerebellum")
list_0$Regulation_limma_Cerebrum <- ifelse(list_0$Cerebrum<0, "DownRegulated Cerebrum", "UpRegulated Cerebrum")
list_0$Regulation_limma_Heart <- ifelse(list_0$Heart<0, "DownRegulated Heart", "UpRegulated Heart")
list_0$Regulation_limma <- paste0(list_0$Regulation_limma_Cerebrum,", ",list_0$Regulation_limma_Cerebellum,", ",list_0$Regulation_limma_Astrocyte,", ",list_0$Regulation_limma_Heart)
list_0 <- list_0[,c(1:52,57)]

#############################
#annotate chr using ENSEMBL
anotate_chr.function <- function(list_0,annotation){
  ensembl =useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "hsapiens_gene_ensembl", mirror = "useast")  
  attributes = listAttributes(ensembl)
  inputNames <- list_0$SYMBOL
  
  #interesting attributes
  attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])))
  attrib <-attributes[attrib,1] # "ensembl_gene_id"    "external_gene_name" "entrezgene" 
  annotation = getBM(attributes=attrib,filters="external_gene_name",values=inputNames, mart=ensembl)
  colnames(annotation) <- c("SYMBOL","chr")
  chrSel <- c(seq(1:23), "X","Y","MT")
  annotation <- annotation[order(match(annotation$chr, chrSel )),]
  annotation <-  annotation[!duplicated(annotation$SYMBOL),] #all the duplicated genes are in scaffolds
  list_0 <- left_join(list_0,annotation,by="SYMBOL")
  as.data.frame(unique(list_0[,c('GeneName_transcript','chr')]) %>% group_by(chr)  %>% count())  #392 wo chr annotation in adult/fetal samples

  #reannotation absed in entrezid for the genes that failed the chr asignation base don gene name
  reAnnotate <- unique(list_0[is.na(list_0$chr),'ENTREZID'])
  toReAnnotate <- unique(list_0[is.na(list_0$chr),])
  toReAnnotate$chr <- NULL

  attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='entrezgene_id'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])))
  attrib <-attributes[attrib,1] # "ensembl_gene_id"    "external_gene_name" "entrezgene" 
  annotation2 = getBM(attributes=attrib,filters="entrezgene_id",values=reAnnotate, mart=ensembl)
  colnames(annotation2) <- c("ENTREZID","chr")
  annotation2 <-  annotation2[!duplicated(annotation2$ENTREZID),] #all the duplicated genes are in scaffolds
  annotation2$ENTREZID <- as.character(annotation2$ENTREZID)

  toReAnnotate <- left_join(toReAnnotate,annotation2,by="ENTREZID")
  list_0 <- list_0[!is.na(list_0$chr),] 
  list_0 <- rbind(list_0, toReAnnotate)

  #remove genes wo chr asignation by ensembl
  list_0 <- list_0[!is.na(list_0$chr),]#removed 202 genes cause change in technology

  #stats
  #length(unique(list_0$SYMBOL))
  #length(unique(list_0$GeneName_transcript)) 

  #for both analyses

  list_0 <- list_0[order(list_0$adj.P.Val,decreasing = TRUE),]
  write.table(list_0, file=paste0(wd,f_results, f_norm,f_limma,f_name, "list_all_EGs_F.txt"), col.names=T,row.names=T, sep="\t",quote=FALSE );
  save(list_0, file=paste0(wd,f_results, f_norm,f_limma,f_name, "list_all_EGs_F.RData"));
  save(list_0, file=paste0(wd,f_results, f_norm,f_limma,f_name,name,"_EGs_DEGs_full_annotation_limma.RData"));
  assign("list_0",list_0,.GlobalEnv);
};

anotate_chr.function(list_0,annotation)

degs.fdr05<- subset(list_0, list_0$adj.P.Val < 0.5);
degs.fdr005<- subset(list_0, list_0$adj.P.Val < 0.05);
degs.fdr01<- subset(list_0, list_0$adj.P.Val < 0.1);
degs.fdr0005<- subset(list_0, list_0$adj.P.Val < 0.005);

Ordered.degs.fdr0005 <- degs.fdr0005[order(degs.fdr0005$chr,decreasing = TRUE),]
Ordered.degs.fdr005 <- degs.fdr005[order(degs.fdr005$chr,decreasing = TRUE),]
Ordered.degs.fdr05 <- degs.fdr05[order(degs.fdr05$chr,decreasing = TRUE),]
Ordered.degs.fdr01 <- degs.fdr01[order(degs.fdr01$chr,decreasing = TRUE),]


#############
## fetal
########
# dim(degs.fdr0005 )
#14 54
# dim(degs.fdr005) 
#38 54
# dim(degs.fdr01)
#58 54
# dim(degs.fdr05)
#2624   54
#

allUp <- "UpRegulated Cerebrum, UpRegulated Cerebellum, UpRegulated Astrocyte, UpRegulated Heart"
NervSystemUp <- "UpRegulated Cerebrum, UpRegulated Cerebellum, UpRegulated Astrocyte,"
allDown <- "DownRegulated Cerebrum, DownRegulated Cerebellum, DownRegulated Astrocyte, DownRegulated Heart"
NervSystemDown <- "DownRegulated Cerebrum, DownRegulated Cerebellum, DownRegulated Astrocyte,"
##############################################
#for fetal samples:4 tissues so 4 contrast
##############################################

data.DEA.results<- data.frame(data=paste0("DEGs_limma_", name))
chr21.degs.fdr0005 <- degs.fdr0005[which(degs.fdr0005$chr==21),]
chr21.degs.fdr05 <- degs.fdr05[which(degs.fdr05$chr==21),]
chr21.degs.fdr005 <- degs.fdr005[which(degs.fdr005$chr==21),]
chr21.degs.fdr01 <- degs.fdr01[which(degs.fdr01$chr==21),]

data.DEA.results$`DEA_0.005` <- as.numeric(length(unique(degs.fdr0005[,2])))
data.DEA.results$`probes_0.005` <- paste(as.character(degs.fdr0005[,1]),collapse=" ")
data.DEA.results$`GeneSymbol_0.005` <-paste(as.character(degs.fdr0005[,2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_all_models` <- length(degs.fdr0005[which(degs.fdr0005$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.005_all_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_all_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.005_all_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_nervousSystem_models` <- length(degs.fdr0005[grep(NervSystemUp,degs.fdr0005$Regulation_limma),2])
data.DEA.results$`Upregulated_0.005_nervousSystem_models` <-paste(as.character(degs.fdr0005[grep(NervSystemUp,degs.fdr0005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_nervousSystem_models_chr21` <- length(chr21.degs.fdr0005[grep(NervSystemUp,chr21.degs.fdr0005$Regulation_limma),2])
data.DEA.results$`Upregulated_0.005_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr0005[grep(NervSystemUp,chr21.degs.fdr0005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_Cerebrum_models` <- length(degs.fdr0005[which(degs.fdr0005$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.005_Cerebrum_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_Cerebrum_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.005_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_Cerebellum_models` <- length(degs.fdr0005[which(degs.fdr0005$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.005_Cerebellum_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_Cerebellum_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.005_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_Astrocyte_models` <- length(degs.fdr0005[which(degs.fdr0005$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.005_Astrocyte_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Astrocyte>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.005_Astrocyte_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.005_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Astrocyte>0),2]),collapse=" ")

data.DEA.results$`Number_downregulated_0.005_all_models` <- length(degs.fdr0005[which(degs.fdr0005$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.005_all_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_all_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.005_all_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_nervousSystem_models` <- length(degs.fdr0005[grep(NervSystemDown, degs.fdr0005$Regulation_limma),2])
data.DEA.results$`Downregulated_0.005_nervousSystem_models` <-paste(as.character(degs.fdr0005[grep(NervSystemDown,degs.fdr0005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_nervousSystem_models_chr21` <- length(chr21.degs.fdr0005[grep(NervSystemDown,chr21.degs.fdr0005$Regulation_limma),2])
data.DEA.results$`Downregulated_0.005_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr0005[grep(NervSystemDown,chr21.degs.fdr0005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_Cerebrum_models` <- length(degs.fdr0005[which(degs.fdr0005$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.005_Cerebrum_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_Cerebrum_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.005_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_Cerebellum_models` <- length(degs.fdr0005[which(degs.fdr0005$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.005_Cerebellum_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_Cerebellum_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.005_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_Astrocyte_models` <- length(degs.fdr0005[which(degs.fdr0005$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.005_Astrocyte_models` <-paste(as.character(degs.fdr0005[which(degs.fdr0005$Astrocyte<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.005_Astrocyte_models_chr21` <- length(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.005_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr0005[which(chr21.degs.fdr0005$Astrocyte<0),2]),collapse=" ")

data.DEA.results$`DEA_0.05` <- as.numeric(length(unique(degs.fdr005[,2])))
data.DEA.results$`probes_0.05` <- paste(as.character(degs.fdr005[,1]),collapse=" ")
data.DEA.results$`GeneSymbol_0.05` <-paste(as.character(degs.fdr005[,2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_all_models` <- length(degs.fdr005[which(degs.fdr005$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.05_all_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_all_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.05_all_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_nervousSystem_models` <- length(degs.fdr005[grep(NervSystemUp,degs.fdr005$Regulation_limma),2])
data.DEA.results$`Upregulated_0.05_nervousSystem_models` <-paste(as.character(degs.fdr005[grep(NervSystemUp,degs.fdr005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_nervousSystem_models_chr21` <- length(chr21.degs.fdr005[grep(NervSystemUp,chr21.degs.fdr005$Regulation_limma),2])
data.DEA.results$`Upregulated_0.05_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr005[grep(NervSystemUp,chr21.degs.fdr005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_Cerebrum_models` <- length(degs.fdr005[which(degs.fdr005$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.05_Cerebrum_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_Cerebrum_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.05_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_Cerebellum_models` <- length(degs.fdr005[which(degs.fdr005$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.05_Cerebellum_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_Cerebellum_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.05_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_Astrocyte_models` <- length(degs.fdr005[which(degs.fdr005$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.05_Astrocyte_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Astrocyte>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.05_Astrocyte_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.05_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Astrocyte>0),2]),collapse=" ")

data.DEA.results$`Number_downregulated_0.05_all_models` <- length(degs.fdr005[which(degs.fdr005$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.05_all_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_all_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.05_all_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_nervousSystem_models` <- length(degs.fdr005[grep(NervSystemDown, degs.fdr005$Regulation_limma),2])
data.DEA.results$`Downregulated_0.05_nervousSystem_models` <-paste(as.character(degs.fdr005[grep(NervSystemDown,degs.fdr005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_nervousSystem_models_chr21` <- length(chr21.degs.fdr005[grep(NervSystemDown,chr21.degs.fdr005$Regulation_limma),2])
data.DEA.results$`Downregulated_0.05_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr005[grep(NervSystemDown,chr21.degs.fdr005$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_Cerebrum_models` <- length(degs.fdr005[which(degs.fdr005$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.05_Cerebrum_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_Cerebrum_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.05_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_Cerebellum_models` <- length(degs.fdr005[which(degs.fdr005$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.05_Cerebellum_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_Cerebellum_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.05_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_Astrocyte_models` <- length(degs.fdr005[which(degs.fdr005$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.05_Astrocyte_models` <-paste(as.character(degs.fdr005[which(degs.fdr005$Astrocyte<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.05_Astrocyte_models_chr21` <- length(chr21.degs.fdr005[which(chr21.degs.fdr005$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.05_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr005[which(chr21.degs.fdr005$Astrocyte<0),2]),collapse=" ")

data.DEA.results$`DEA_0.1` <- as.numeric(length(unique(degs.fdr01[,2])))
data.DEA.results$`probes_0.1` <- paste(as.character(degs.fdr01[,1]),collapse=" ")
data.DEA.results$`GeneSymbol_0.1` <-paste(as.character(degs.fdr01[,2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_all_models` <- length(degs.fdr01[which(degs.fdr01$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.1_all_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_all_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.1_all_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_nervousSystem_models` <- length(degs.fdr01[grep(NervSystemUp, degs.fdr01$Regulation_limma),2])
data.DEA.results$`Upregulated_0.1_nervousSystem_models` <-paste(as.character(degs.fdr01[grep(NervSystemUp, degs.fdr01$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_nervousSystem_models_chr21` <- length(chr21.degs.fdr01[grep(NervSystemUp,chr21.degs.fdr01$Regulation_limma),2])
data.DEA.results$`Upregulated_0.1_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr01[grep(NervSystemUp,chr21.degs.fdr01$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_Cerebrum_models` <- length(degs.fdr01[which(degs.fdr01$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.1_Cerebrum_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_Cerebrum_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.1_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_Cerebellum_models` <- length(degs.fdr01[which(degs.fdr01$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.1_Cerebellum_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_Cerebellum_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.1_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_Astrocyte_models` <- length(degs.fdr01[which(degs.fdr01$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.1_Astrocyte_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Astrocyte>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.1_Astrocyte_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.1_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Astrocyte>0),2]),collapse=" ")

data.DEA.results$`Number_downregulated_0.1_all_models` <- length(degs.fdr01[which(degs.fdr01$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.1_all_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_all_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.1_all_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_nervousSystem_models` <- length(degs.fdr01[grep(NervSystemDown, degs.fdr01$Regulation_limma),2])
data.DEA.results$`Downregulated_0.1_nervousSystem_models` <-paste(as.character(degs.fdr01[grep(NervSystemDown,degs.fdr01$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_nervousSystem_models_chr21` <- length(chr21.degs.fdr01[grep(NervSystemDown,chr21.degs.fdr01$Regulation_limma),2])
data.DEA.results$`Downregulated_0.1_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr01[grep(NervSystemDown,chr21.degs.fdr01$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_Cerebrum_models` <- length(degs.fdr01[which(degs.fdr01$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.1_Cerebrum_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_Cerebrum_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.1_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_Cerebellum_models` <- length(degs.fdr01[which(degs.fdr01$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.1_Cerebellum_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_Cerebellum_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.1_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_Astrocyte_models` <- length(degs.fdr01[which(degs.fdr01$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.1_Astrocyte_models` <-paste(as.character(degs.fdr01[which(degs.fdr01$Astrocyte<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.1_Astrocyte_models_chr21` <- length(chr21.degs.fdr01[which(chr21.degs.fdr01$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.1_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr01[which(chr21.degs.fdr01$Astrocyte<0),2]),collapse=" ")


data.DEA.results$`DEA_0.5` <- as.numeric(length(unique(degs.fdr05[,2])))
data.DEA.results$`probes_0.5` <- paste(as.character(degs.fdr05[,1]),collapse=" ")
data.DEA.results$`GeneSymbol_0.5` <-paste(as.character(degs.fdr05[,2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_all_models` <- length(degs.fdr05[which(degs.fdr05$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.5_all_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_all_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Regulation_limma==allUp),2])
data.DEA.results$`Upregulated_0.5_all_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Regulation_limma==allUp),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_nervousSystem_models` <- length(degs.fdr05[grep(NervSystemUp,degs.fdr05$Regulation_limma),2])
data.DEA.results$`Upregulated_0.5_nervousSystem_models` <-paste(as.character(degs.fdr05[grep(NervSystemUp,degs.fdr05$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_nervousSystem_models_chr21` <- length(chr21.degs.fdr05[grep(NervSystemUp,chr21.degs.fdr05$Regulation_limma),2])
data.DEA.results$`Upregulated_0.5_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr05[grep(NervSystemUp,chr21.degs.fdr05$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_Cerebrum_models` <- length(degs.fdr05[which(degs.fdr05$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.5_Cerebrum_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_Cerebrum_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebrum>0),2])
data.DEA.results$`Upregulated_0.5_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebrum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_Cerebellum_models` <- length(degs.fdr05[which(degs.fdr05$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.5_Cerebellum_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_Cerebellum_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebellum>0),2])
data.DEA.results$`Upregulated_0.5_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebellum>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_Astrocyte_models` <- length(degs.fdr05[which(degs.fdr05$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.5_Astrocyte_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Astrocyte>0),2]),collapse=" ")
data.DEA.results$`Number_upregulated_0.5_Astrocyte_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Astrocyte>0),2])
data.DEA.results$`Upregulated_0.5_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Astrocyte>0),2]),collapse=" ")

data.DEA.results$`Number_downregulated_0.5_all_models` <- length(degs.fdr05[which(degs.fdr05$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.5_all_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_all_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Regulation_limma==allDown),2])
data.DEA.results$`Downregulated_0.5_all_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Regulation_limma==allDown),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_nervousSystem_models` <- length(degs.fdr05[grep(NervSystemDown, degs.fdr05$Regulation_limma),2])
data.DEA.results$`Downregulated_0.5_nervousSystem_models` <-paste(as.character(degs.fdr05[grep(NervSystemDown,degs.fdr05$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_nervousSystem_models_chr21` <- length(chr21.degs.fdr05[grep(NervSystemDown,chr21.degs.fdr05$Regulation_limma),2])
data.DEA.results$`Downregulated_0.5_nervousSystem_models_chr21` <-paste(as.character(chr21.degs.fdr05[grep(NervSystemDown,chr21.degs.fdr05$Regulation_limma),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_Cerebrum_models` <- length(degs.fdr05[which(degs.fdr05$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.5_Cerebrum_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_Cerebrum_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebrum<0),2])
data.DEA.results$`Downregulated_0.5_Cerebrum_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebrum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_Cerebellum_models` <- length(degs.fdr05[which(degs.fdr05$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.5_Cerebellum_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_Cerebellum_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebellum<0),2])
data.DEA.results$`Downregulated_0.5_Cerebellum_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Cerebellum<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_Astrocyte_models` <- length(degs.fdr05[which(degs.fdr05$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.5_Astrocyte_models` <-paste(as.character(degs.fdr05[which(degs.fdr05$Astrocyte<0),2]),collapse=" ")
data.DEA.results$`Number_downregulated_0.5_Astrocyte_models_chr21` <- length(chr21.degs.fdr05[which(chr21.degs.fdr05$Astrocyte<0),2])
data.DEA.results$`Downregulated_0.5_Astrocyte_models_chr21` <-paste(as.character(chr21.degs.fdr05[which(chr21.degs.fdr05$Astrocyte<0),2]),collapse=" ")


write.table(all_info_degs05_limma, file=paste0(wd, f_results, f_norm,f_limma,f_name,f_DEGs05, "list_limma_DEGs05.txt"), col.names=T,row.names=F, sep="\t",quote=FALSE )
write.table(all_info_degs005_limma, file=paste0(wd, f_results, f_norm,f_limma,f_name,f_DEGs005, "list_limma_DEGS005.txt"), col.names=T,row.names=F, sep="\t",quote=FALSE )
write.table(all_info_degs01_limma, file=paste0(wd, f_results, f_norm,f_limma,f_name,f_DEGs01, "list_limma_DEGS01.txt"), col.names=T,row.names=F, sep="\t",quote=FALSE )
write.table(all_info_degs0005_limma, file=paste0(wd, f_results, f_norm,f_limma,f_name,f_DEGs0005, "list_limma_DEGS0005.txt"), col.names=T,row.names=F, sep="\t",quote=FALSE )


EGs_limma <- unique(data.frame(swissprot=list_0$SYMBOL, ensembleGeneName=list_0$ENSEMBLTRANS,GeneName_transcript=list_0$GeneName_transcript ))
EGs_Name_limma <- unique(data.frame(swissprot=list_0$SYMBOL))
DEGs_Name_limma005 <- unique(data.frame(swissprot=all_info_degs005_limma$SYMBOL))
DEGs_Name_limma0005 <- unique(data.frame(swissprot=all_info_degs0005_limma$SYMBOL))
save(DEGs_Name_limma0005, file=paste0(wd,f_results, f_norm,f_limma,f_name,f_DEGs0005, name,"_DEGs0005_limma_geneNames.RData"))
save(DEGs_Name_limma005, file=paste0(wd,f_results, f_norm,f_limma,f_name,f_DEGs005, name,"_DEGs005_limma_geneNames.RData"))
save(EGs_Name_limma, file=paste0(wd,f_results, f_norm,f_limma,f_name, name,"_EGs_limma_geneNamesInfo.RData"))
save(EGs_limma, file=paste0(wd,f_results, f_norm,f_limma,f_name,name,"_EGs_limma_transcript_annotation.RData"))
save(Ordered.degs.fdr0005,Ordered.degs.fdr005,Ordered.degs.fdr01,Ordered.degs.fdr05,list_0,data.DEA.results, file=paste0(wd,f_results, f_norm,f_limma,f_name,name,"_DEA_limma_all_final_objects.RData"))

sumupStatsLimma <- data.DEA.results[,c(grep("data",colnames(data.DEA.results)),grep("DEA",colnames(data.DEA.results)),grep("Number",colnames(data.DEA.results)))]  

save.image(file=paste0(wd,f_results, f_Rdata,name,"_session.RData"))



################################################
########################################
##fcros DEGs analysis
########################################
################################################
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

fcros_DEA.function<- function(eset_df,genoControlName,genoMutName, trim.opt,list_0,name,annotation,annotation2,annotationUniq){
	# This function perform the DEA using fcros package 
	# Doc: https://cran.r-project.org/web/packages/fcros/index.html
	# 	   https://www.ncbi.nlm.nih.gov/pubmed/24423217
	#  The input counts for fcros, should be the normalised counts, after adjusting by library size, possibly gene length,
	#  batch correction if needed. 
	#1. Using the normalised counts by Deseq2
	#2. We identify the DEGs by fcros, the up and regulated ones, and merge into one df. 
	#3. We also annotate these fcrosDEGs with the DESeq2 stats(FC, basemean, pval, fdr...)
	#4. This function also produces one of the input files needed to run gage data.all.deseq.fcros
	data <-  eset_df  #input data the normalised counts by deseq2 in log2
	ncolData <- ncol(data)		
	cont <- colnames(data[,grep(genoControlName,colnames(data))])
	print(cont)
	test <- colnames(data[,grep(genoMutName,colnames(data))])
	print(test)
	log2.opt <- 0 #data in the matrix "data" are expressed in a log2 scale
	data <- left_join(data,annotationUniq,by="PROBEID") 
	index <- 1:nrow(data)
	data$index <- index
	data <- data[, c("PROBEID" ,cont,test,"SYMBOL","index","ENSEMBLTRANS")]
	af <- fcros(data, cont, test, log2.opt, trim.opt) #The first column of the matrix "xdata" should contain the gene IDs or their names.
	af$fdr <- p.adjust(af$p.value, method="BH")
	dataStats <- data.frame(minPval=min(af$p.value),maxPval=max(af$p.value),minFDR=min(af$fdr),maxFDR=max(af$fdr))
	##############################################################
	cuack <-af
	cuack$bounds <- NULL
	cuack$params <- NULL
	cuack$params_t <- NULL
	afDf <- do.call(cbind.data.frame, cuack)
	colnames(afDf)[1] <- "PROBEID"
	afDf$bounds <- paste0(unlist(af[[7]][1]),":",unlist(af[[7]][2]))
	afDf$params <- paste0(unlist(af[[8]][1]),":",unlist(af[[8]][2]),":",unlist(af[[8]][3]))
	afDf$params_t <- paste0(unlist(af[[9]][1]),":",unlist(af[[9]][2]),":",unlist(af[[9]][3]))
	afDf <- left_join(afDf,annotationUniq,by="PROBEID") 
	afDf <- unique(afDf[,c(1,11:14,2:10)])
	afDf<- afDf[!is.na(afDf[,2]),]
	length(unique(afDf[,2]))  #13768   #total expressed genes annotated with gene Name
	
	#annotating with chr
	afDf.tmp <- left_join(afDf,annotation,by="SYMBOL")
	#reannotation absed in entrezid for the genes that failed the chr asignation base don gene name
	reAnnotate <- unique(afDf.tmp[is.na(afDf.tmp$chr),'ENTREZID'])
	toReAnnotate <- unique(afDf.tmp[is.na(afDf.tmp$chr),])
	toReAnnotate$chr <- NULL
	toReAnnotate <- left_join(toReAnnotate,annotation2,by="ENTREZID")
	afDf.tmp <- afDf.tmp[!is.na(afDf.tmp$chr),] 
	afDf.tmp <- rbind(afDf.tmp, toReAnnotate)
	#remove genes wo chr asignation by ensembl
	afDf <- afDf.tmp[!is.na(afDf.tmp$chr),] #13566   #total expressed genes annotated with gene Name and CHR 
	totalEGs <- length(unique(afDf[,2]))  #13566   #total expressed genes annotated with gene Name

	dupGenes <- unique(afDf[duplicated(afDf[,2]),2])
#	dataDup <- list()
#	for (i in 1:length(dupGenes)){
#		dataDup[[i]] <-afDf[which(afDf[,2]==dupGenes[i]),]
#		names(dataDup)[i] <- dupGenes[i]
#		i=i+1
#	};

	#############################################################
	#B. IDENTIFYING UP AND DOWN REGULATED DEGs
	#############################################################
	#############################################################
	#ading the alpha values manually, as it cannot be computed by the fcrostopN functionwith so many DEGs i our list.
	# So as this was the threshold defined by Dembele D. to be used equal to get a fdr <0,05

	f.value <- afDf[,'f.value']
	id.up <- matrix(0, 1)
	id.down <- matrix(0, 1)
	alpha_up <- 0.975
	alpha_dn <- 0.025

	down.all <- 1;
	up.all <- 1;
	for (i in 1:totalEGs) {
		if (f.value[i] <= alpha_dn) { id.down[down.all] <- i; down.all <- down.all + 1; }
		if (f.value[i] >= alpha_up) { id.up[up.all] <- i; up.all <- up.all + 1; }
	};
	data.down_all <- afDf[id.down[1:(down.all-1)], ];
	ndown_all <- nrow(data.down_all);
	data.up_all <- afDf[id.up[1:(up.all-1)], ];
	nup_all <- nrow(data.up_all) 
	nup_all_chr21 <- nrow(unique(data.up_all[which(data.up_all$chr==21),]))
	ndown_all_chr21 <- nrow(unique(data.down_all[which(data.down_all$chr==21),]))
	
	fcrosDEGsStats <- data.frame(totalDEGs=ndown_all + nup_all,upregulatedDEGs=nup_all,downregulatedDEGs=ndown_all,
		totalDEGs_chr21=nup_all_chr21 + ndown_all_chr21,upregulatedDEGs_chr21=nup_all_chr21,downregulatedDEGs_chr21=ndown_all_chr21)
	data.down_all$Regulation <- rep("DownRegulated", times=nrow(data.down_all))
	data.up_all$Regulation <- rep("UpRegulated", times=nrow(data.up_all))
	data.down_all$FoldChange <- as.numeric(2^data.down_all$FC2)*-1
	data.up_all$FoldChange <- as.numeric(2^data.up_all$FC2)
	options(scipen = 999) #turn off scientific notation
	
	data.all <- rbind(data.up_all,data.down_all)
	data.all <- data.all[,c(1:2,16,3:5,15,17,6:14)]
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
	};
	


	#############################################################
	#creating output files
	#############################################################

	if(name=="human_fetal"){
		res.limma <-data.frame(list_0)
		colnames(res.limma)[c(7:52)]<- paste0(colnames(res.limma)[c(7:52)],"_limma")
		colnames(res.limma)
		colnames(data.all)[c(3,8:17)] <- paste0(colnames(data.all)[c(3,8:17)],"_fcros")
		data.all.limma.fcros <- left_join(res.limma,data.all,by=c("ENSEMBLTRANS","PROBEID","ENTREZID","GENENAME","SYMBOL","chr"))		
		colnames(data.all.limma.fcros)[which(colnames(data.all.limma.fcros)=="ENSEMBLTRANS")] <- "EnsembleGeneName"
		colnames(data.all.limma.fcros)[which(colnames(data.all.limma.fcros)=="SYMBOL")] <- "swissprot"
		data.all.limma.fcros <- data.all.limma.fcros[gtools::mixedorder(as.character(data.all.limma.fcros$chr)),]
		nameTested <- gsub("T21_","",genoMutName)
		
		write.table(data.down_all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_down, "DEA_fcros_downregulated_",name,"_",nameTested,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		write.table(data.up_all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_up, "DEA_fcros_upregulated_",name,"_",nameTested,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		write.table(data.all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_all, "DEA_fcros_",name,"_",nameTested,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		#write.table(data.all.deseq.fcros, file=paste0(wd, f_results, f_norm, "DEA_fcros_deseq2_comparison",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		assign(paste0("data.up_all_",nameTested),data.up_all,.GlobalEnv)
		assign(paste0("data.down_all_",nameTested),data.down_all,.GlobalEnv)
		assign(paste0("res.limma_",nameTested),res.limma,.GlobalEnv)
		assign(paste0("afDf_",nameTested),afDf,.GlobalEnv) #needed for volcano and correlation plot is all EGs 
		assign(paste0("dataStats_",nameTested),dataStats,.GlobalEnv)
		assign(paste0("totalEGs_",nameTested),totalEGs,.GlobalEnv)
		assign(paste0("fcrosDEGsStats_",nameTested),fcrosDEGsStats,.GlobalEnv)
		#assign(paste0("doubleRegGenes_",nameTested),doubleRegGenes,.GlobalEnv)
		assign(paste0("data.all.limma.fcros_",nameTested),data.all.limma.fcros,.GlobalEnv)
		save(afDf,data.all.limma.fcros,file=paste0(wd, f_results,f_Rdata,name,"_",nameTested,"_tmp_Input_for_GAGE.RData")) #data.all.deseq.fcros, 

	}	else{
		res.limma <-data.frame(list_0)
		colnames(res.limma)[c(7:31)]<- paste0(colnames(res.limma)[c(7:31)],"_limma")
		colnames(res.limma)
		colnames(data.all)[c(3,8:17)] <- paste0(colnames(data.all)[c(3,8:17)],"_fcros")
		data.all.limma.fcros <- left_join(res.limma,data.all,by=c("ENSEMBLTRANS","PROBEID","ENTREZID","GENENAME","SYMBOL","chr"))
		
		colnames(data.all.limma.fcros)[which(colnames(data.all.limma.fcros)=="ENSEMBLTRANS")] <- "EnsembleGeneName"
		colnames(data.all.limma.fcros)[which(colnames(data.all.limma.fcros)=="SYMBOL")] <- "swissprot"
		data.all.limma.fcros <- data.all.limma.fcros[gtools::mixedorder(as.character(data.all.limma.fcros$chr)),]


		write.table(data.down_all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_down, "DEA_fcros_downregulated_",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		write.table(data.up_all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_up, "DEA_fcros_upregulated_",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		write.table(data.all, file=paste0(wd, f_results, f_norm,f_fcros, f_name,f_all, "DEA_fcros_",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		#write.table(data.all.deseq.fcros, file=paste0(wd, f_results, f_norm, "DEA_fcros_deseq2_comparison",name,".txt"), sep = "\t", col.names=TRUE, row.names=F,na = " ",)
		assign("data.up_all",data.up_all,.GlobalEnv)
		assign("data.down_all",data.down_all,.GlobalEnv)
		assign("res.limma",res.limma,.GlobalEnv)
		assign("afDf",afDf,.GlobalEnv) #needed for volcano and correlation plot is all EGs 
		assign("dataStats",dataStats,.GlobalEnv)
		assign("totalEGs",totalEGs,.GlobalEnv)
		assign("fcrosDEGsStats",fcrosDEGsStats,.GlobalEnv)
		#assign("doubleRegGenes",doubleRegGenes,.GlobalEnv)
		assign("data.all.limma.fcros",data.all.limma.fcros,.GlobalEnv)
		save(afDf,data.all.limma.fcros,file=paste0(wd, f_results,f_Rdata,name,"_tmp_Input_for_GAGE.RData")) #data.all.deseq.fcros, 	
		}	
};


f_fcros<- "fcros/"
trim.opt <- 0.25


###########
#fetal
###########
genoControlName <- "Euploid_Astrocyte"
genoMutName <- "T21_Astrocyte"
fcros_DEA.function(eset_df,genoControlName,genoMutName, trim.opt,list_0,name,annotation,annotation2,annotationUniq)
fcrosDEGsStats_Astrocyte
head(data.all.limma.fcros_Astrocyte)
#totalDEGs upregulatedDEGs downregulatedDEGs 	totalDEGs_chr21	  upregulatedDEGs_chr21 downregulatedDEGs_chr21
#       305             186               119			28   				27                      1 (ITGB2)


genoControlName <- "Euploid_Cerebrum"
genoMutName <- "T21_Cerebrum"
fcros_DEA.function(eset_df,genoControlName,genoMutName, trim.opt,list_0,name,annotation,annotation2,annotationUniq)
fcrosDEGsStats_Cerebrum
#  totalDEGs upregulatedDEGs downregulatedDEGs 	totalDEGs_chr21	  upregulatedDEGs_chr21 downregulatedDEGs_chr21
#       426             339                87			44   				44                      0

genoControlName <- "Euploid_Cerebellum"
genoMutName <- "T21_Cerebellum"
fcros_DEA.function(eset_df,genoControlName,genoMutName, trim.opt,list_0,name,annotation,annotation2,annotationUniq)
fcrosDEGsStats_Cerebellum
# totalDEGs upregulatedDEGs downregulatedDEGs 	totalDEGs_chr21	  upregulatedDEGs_chr21 downregulatedDEGs_chr21
#       558             230               328			36   				36                      0

genoControlName <- "Euploid_Heart"
genoMutName <- "T21_Heart"
fcros_DEA.function(eset_df,genoControlName,genoMutName, trim.opt,list_0,name,annotation,annotation2,annotationUniq)
fcrosDEGsStats_Heart
# totalDEGs upregulatedDEGs downregulatedDEGs 	totalDEGs_chr21	  upregulatedDEGs_chr21 downregulatedDEGs_chr21
#       279             179               100			27   				27                      0


colnames(fcrosDEGsStats_Astrocyte) <- paste0("Astrocyte: ",colnames(fcrosDEGsStats_Astrocyte))
colnames(fcrosDEGsStats_Cerebrum) <- paste0("Cerebrum: ",colnames(fcrosDEGsStats_Cerebrum))
colnames(fcrosDEGsStats_Cerebellum) <- paste0("Cerebellum: ",colnames(fcrosDEGsStats_Cerebellum))
colnames(fcrosDEGsStats_Heart) <- paste0("Heart: ",colnames(fcrosDEGsStats_Heart))
sumupStatsfcros  <- cbind(fcrosDEGsStats_Astrocyte,fcrosDEGsStats_Cerebrum,fcrosDEGsStats_Cerebellum,fcrosDEGsStats_Heart)



#sumup limma/fcros to later  use for the import EGs function
#adult
#############
data.EGs.limma.fcros <- data.all.limma.fcros[,c(1,33,6,22,23,31,2,37,32,34,35,41,24,29)]
DEGs.limma.fcros <- unique(data.EGs.limma.fcros[,c(7,10)])
DEGsAll <- list(limma_allEGs=list_0,limma_DEGs01=Ordered.degs.fdr01,limma_DEGs05=Ordered.degs.fdr05, limma_DEGs005=Ordered.degs.fdr005,
	limma_DEGs0005=Ordered.degs.fdr0005,LimmaStats=data.DEA.results,sumup_limma=sumupStatsLimma,EGsFcros=data.EGs.limma.fcros,
	DEGsFcros=DEGs.limma.fcros,fcrosStats=sumupStatsfcros)


#fetal
#############
#for EGs, the same in all the contrasts
data.EGs.limma.fcros_fetal <- unique(data.all.limma.fcros_Cerebrum[,c(1,54,6,36,37,52,2,58)])
#for DEGs: depends on the contrast
data.DEGs.limma.fcros_Astrocyte <- unique(data.all.limma.fcros_Astrocyte[,c(1,54,6,32,33,52,2,58,53,55,56,62,44,51)])
data.DEGs.limma.fcros_Astrocyte <- data.DEGs.limma.fcros_Astrocyte[!is.na(data.DEGs.limma.fcros_Astrocyte$Regulation_fcros),]
data.DEGs.limma.fcros_Cerebrum <- unique(data.all.limma.fcros_Cerebrum[,c(1,54,6,36,37,52,2,58,53,55,56,62,46,51)])
data.DEGs.limma.fcros_Cerebrum <- data.DEGs.limma.fcros_Cerebrum[!is.na(data.DEGs.limma.fcros_Cerebrum$Regulation_fcros),]
data.DEGs.limma.fcros_Cerebellum <- unique(data.all.limma.fcros_Cerebellum[,c(1,54,6,34,35,52,2,58,53,55,56,62,45,51)])
data.DEGs.limma.fcros_Cerebellum <- data.DEGs.limma.fcros_Cerebellum[!is.na(data.DEGs.limma.fcros_Cerebellum$Regulation_fcros),]
data.DEGs.limma.fcros_Heart <- unique(data.all.limma.fcros_Heart[,c(1,54,6,38,39,52,2,58,53,55,56,62,47,51)])
data.DEGs.limma.fcros_Heart <- data.DEGs.limma.fcros_Heart[!is.na(data.DEGs.limma.fcros_Heart$Regulation_fcros),]


DEGs.limma.fcros_Heart <- unique(data.DEGs.limma.fcros_Heart[,c(7,10)])
DEGs.limma.fcros_Cerebellum <- unique(data.DEGs.limma.fcros_Cerebellum[,c(7,10)])
DEGs.limma.fcros_Cerebrum <- unique(data.DEGs.limma.fcros_Cerebrum[,c(7,10)])
DEGs.limma.fcros_Astrocyte <- unique(data.DEGs.limma.fcros_Astrocyte[,c(7,10)])


DEGsAll <- list(limma_allEGs=list_0,limma_DEGs01=Ordered.degs.fdr01,limma_DEGs05=Ordered.degs.fdr05, limma_DEGs005=Ordered.degs.fdr005,
	limma_DEGs0005=Ordered.degs.fdr0005,LimmaStats=data.DEA.results,sumup_limma=sumupStatsLimma,EGsFcros=data.EGs.limma.fcros_fetal,
	DEGsFcros_Heart=DEGs.limma.fcros_Heart, DEGsFcros_Cerebellum=DEGs.limma.fcros_Cerebellum, DEGsFcros_Cerebrum=DEGs.limma.fcros_Cerebrum,
	 DEGsFcros_Astrocyte=DEGs.limma.fcros_Astrocyte,fcrosStats=sumupStatsfcros)




#for all 
###########
sumup <- list(sumup_limma=sumupStatsLimma,sumup_fcros=sumupStatsfcros)
write.xlsx(sumup, file = paste0(wd,f_results, f_norm,name, "_DEA_results_sumupStats.xlsx"))

DEGsLimma <- list(allEGs=list_0,DEGs01=Ordered.degs.fdr01,DEGs05=Ordered.degs.fdr05, DEGs005=Ordered.degs.fdr005,DEGs0005=Ordered.degs.fdr0005,stats=data.DEA.results)
library(openxlsx)
write.xlsx(DEGsLimma, file = paste0(wd,f_results, f_norm,f_limma,f_name,name, "_DEA_results.xlsx"))
write.xlsx(DEGsAll, file = paste0(wd,f_results, f_norm,name, "_ALL_DEA_results.xlsx"))

#fetal
###########
save(affyRaw,colData,Ordered.degs.fdr0005,Ordered.degs.fdr005,
	Ordered.degs.fdr01,Ordered.degs.fdr05,all_info_degs005_limma,
	all_info_degs05_limma,all_info_degs01_limma,all_info_degs0005_limma,
	afDf,sumupStatsfcros,eset,eset_df,data.all.limma.fcros_Astrocyte,data.all.limma.fcros_Cerebrum,
	data.all.limma.fcros_Cerebellum,data.all.limma.fcros_Heart,res.limma_Astrocyte,
	res.limma_Cerebrum,res.limma_Cerebellum,res.limma_Heart,
	data.up_all_Astrocyte,data.down_all_Astrocyte,data.up_all_Cerebrum,data.down_all_Cerebrum,
	data.up_all_Cerebellum,data.down_all_Cerebellum,data.up_all_Heart,data.down_all_Heart,data.DEA.results, 
	data.EGs.limma.fcros_fetal,DEGs.limma.fcros_Heart, DEGs.limma.fcros_Cerebellum, DEGs.limma.fcros_Cerebrum, 
	DEGs.limma.fcros_Astrocyte,
	file=paste0(wd, f_results,f_Rdata,name,"_all_Input_for_GAGE.RData")) 
save.image(file=paste0(wd, f_results, f_Rdata, name,"_session.RData"))



#For th DEGs venn diagram

#adult
write.table(DEGs.limma.fcros, file=paste0(wd,f_Result,f_gage,f_dfa,f_RData,name,"_DEGs_regSense.txt"), sep = "\t", col.names=TRUE, row.names=F)
write.table(EGs.limma.fcros, file=paste0(wd,f_inputDEGs,f_name,"_DEGs/list_all_EGs_annotated_plus.txt"), sep = "\t", col.names=TRUE, row.names=F)

#fetal
write.table(DEGs.limma.fcros_Heart, file=paste0(wd,f_Result,f_gage,f_dfa,f_RData,name,"_Heart_DEGs_regSense.txt"), sep = "\t", col.names=TRUE, row.names=F)
write.table(DEGs.limma.fcros_Cerebellum, file=paste0(wd,f_Result,f_gage,f_dfa,f_RData,name,"_Cerebellum_DEGs_regSense.txt"), sep = "\t", col.names=TRUE, row.names=F)
write.table(DEGs.limma.fcros_Cerebrum, file=paste0(wd,f_Result,f_gage,f_dfa,f_RData,name,"_Cerebrum_DEGs_regSense.txt"), sep = "\t", col.names=TRUE, row.names=F)
write.table(DEGs.limma.fcros_Astrocyte, file=paste0(wd,f_Result,f_gage,f_dfa,f_RData,name,"_Astrocyte_DEGs_regSense.txt"), sep = "\t", col.names=TRUE, row.names=F)

write.table(EGs.limma.fcros_fetal, file=paste0(wd,f_inputDEGs,f_name,"_DEGs/list_all_EGs_annotated_plus.txt"), sep = "\t", col.names=TRUE, row.names=F)


################################################################################################################################
################################################################################################################################
############################### END. Finished on the  19/02/2019 ####################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################
################################################################################################################################
################################################################################################################################

sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.5 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] affy_1.56.0            hgu133aprobe_2.18.0    hugene10stv1cdf_2.18.0
 [4] hgu133acdf_2.18.0      pd.hg.u133a_3.12.0     DBI_1.0.0             
 [7] RSQLite_2.1.1          hgu133a.db_3.2.3       org.Hs.eg.db_3.5.0    
[10] AnnotationDbi_1.40.0   biomaRt_2.34.2         fcros_1.5.6           
[13] gdata_2.18.0           apeglm_1.0.3           PoiClaClu_1.0.2.1     
[16] amap_0.8-16            dendextend_1.8.0       GenomicRanges_1.30.3  
[19] GenomeInfoDb_1.14.0    limma_3.34.9           scales_1.0.0          
[22] genefilter_1.60.0      RColorBrewer_1.1-2     gplots_3.0.1.1        
[25] calibrate_1.7.2        MASS_7.3-51.4          nlme_3.1-140          
[28] ff_2.2-14              bit_1.1-14             oligo_1.42.0          
[31] Biostrings_2.46.0      XVector_0.18.0         IRanges_2.12.0        
[34] S4Vectors_0.16.0       oligoClasses_1.40.0    dplyr_0.7.6           
[37] tidyr_0.8.1            ggfortify_0.4.7        ggplot2_3.1.0         
[40] GEOquery_2.46.15       Biobase_2.38.0         BiocGenerics_0.24.0   
[43] gtools_3.8.1           reshape2_1.4.3         openxlsx_4.1.0.1      
[46] BiocInstaller_1.28.0  

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1           class_7.3-15              
 [3] modeltools_0.2-22          mclust_5.4.2              
 [5] affyio_1.48.0              flexmix_2.3-15            
 [7] bit64_0.9-7                mvtnorm_1.0-8             
 [9] xml2_1.2.0                 codetools_0.2-16          
[11] splines_3.4.4              robustbase_0.93-2         
[13] annotate_1.56.2            cluster_2.1.0             
[15] kernlab_0.9-27             httr_1.4.0                
[17] readr_1.3.1                compiler_3.4.4            
[19] assertthat_0.2.1           Matrix_1.2-17             
[21] lazyeval_0.2.2             prettyunits_1.0.2         
[23] tools_3.4.4                bindrcpp_0.2.2            
[25] coda_0.19-2                gtable_0.3.0              
[27] glue_1.3.1                 GenomeInfoDbData_1.0.0    
[29] affxparser_1.50.0          Rcpp_1.0.1                
[31] bbmle_1.0.20               trimcluster_0.1-2.1       
[33] preprocessCore_1.40.0      iterators_1.0.10          
[35] fpc_2.1-11.1               stringr_1.4.0             
[37] XML_3.98-1.20              DEoptimR_1.0-8            
[39] zlibbioc_1.24.0            hms_0.4.2                 
[41] SummarizedExperiment_1.8.1 memoise_1.1.0             
[43] gridExtra_2.3              emdbook_1.3.11            
[45] stringi_1.2.4              foreach_1.4.4             
[47] caTools_1.17.1.1           zip_2.0.2                 
[49] rlang_0.2.2                pkgconfig_2.0.2           
[51] prabclus_2.2-7             matrixStats_0.54.0        
[53] bitops_1.0-6               lattice_0.20-38           
[55] purrr_0.2.5                bindr_0.1.1               
[57] tidyselect_0.2.5           plyr_1.8.4                
[59] magrittr_1.5               R6_2.4.0                  
[61] DelayedArray_0.4.1         pillar_1.3.0              
[63] whisker_0.3-2              withr_2.1.2               
[65] survival_2.44-1.1          RCurl_1.95-4.11           
[67] nnet_7.3-12                tibble_1.4.2              
[69] crayon_1.3.4               KernSmooth_2.23-15        
[71] progress_1.2.2             viridis_0.5.1             
[73] grid_3.4.4                 blob_1.1.1                
[75] digest_0.6.17              diptest_0.75-7            
[77] xtable_1.8-4               numDeriv_2016.8-1.1       
[79] munsell_0.5.0              viridisLite_0.3.0

