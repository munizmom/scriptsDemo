###########################################################################################
###########################################################################################
##  Script to represent the FC along the chr duplicated regions for the transcriptomics  ## 
##  data  pipeline     																	 ##
##   																				     ##
###########################################################################################
###########################################################################################
#
# Done by Mar Muniz. YHerault team. Last update 010219.
##########################################################################################
# NOTES:
###########################################################################################
# Be sure that you are working with an up to date version of R and bioconductor environment,
# The hidden parts of the code are there for QC purposes or to adapt the script to multiple
# analysis
############################################################################################

############################################################################################
##################################### Packages needed ######################################
############################################################################################
#update R https://www.r-bloggers.com/updating-r-on-ubuntu/

#############################################
############# Packages needed ###############
source("http://bioconductor.org/biocLite.R");
#update.packages(ask = FALSE)
library("PoiClaClu");
library("DESeq2");
library("apeglm");
library("genefilter");
library("GenomicRanges");
library("tidyr");
library("dplyr");
library("gtools");
library("gdata");
library("biomaRt");
library("plotly");
library("rgl");
library("plotrix");
library("grid");
library("RColorBrewer");
library("gtable");
library("DESeq2");
require("vsn");
require("hexbin");
library("scatterplot3d");
library("ggplot2");
library("DESeq2");
library("gridGraphics");
library("cowplot");
library("ggpubr");
library("gridBase");
library("gridExtra");
library("gtable");
library("gplots");
library("grid");
library("xlsx");
library("gage");
library("biomaRt"); 
library("data.table");
library("mygene");
library("stringr");
library("PoiClaClu");
library("VennDiagram");
library("lettercase");
library("grDevices");
library("openxlsx");
library("psych");

#############################################
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()
#wd:
#####################################################

setwd("/Users/marmmoreno/Desktop/onWorking/"); 
wd <- getwd();

#inputs
f_input <- "/RData/RNASeq/" 
f_inputEGs <- "/RData/fcrosEGs/" 
f_dese2Input <- "/RData/deseq2EGs/"

#output directories
f_results <- "/results/"
f_Rdata <- "RData/"
f_FC<- "FC/"
f_hsa21<- "Hsa21/"
datePlotFC<- "220819"
f_Figure <- paste0("Figure_plot_", datePlotFC)
dir.create(file.path(getwd (),f_results, f_FC), showWarnings = FALSE)
dir.create(file.path(getwd (),f_results, f_FC,f_hsa21), showWarnings = FALSE)
dir.create(file.path(getwd (),f_results, f_FC,f_Figure), showWarnings = FALSE)
dir.create(file.path(getwd (),f_results, f_FC,f_Rdata), showWarnings = FALSE)

#################################################################################################################
######################################## CODE ###################################################################
#load sessions
#load(file=paste0(wd, f_results, f_FC,f_Rdata,"rats_FCplot_session.RData"))
#save.image(file=paste0(wd, f_results, f_FC,f_Rdata,"rats_FCplot_session.RData"))


######################################################
###########################
#FC plots
###########################
######################################################
#######################################################################
# Doing the same figure but with all the genes of the HSA21 regions
#######################################################################

nameEGsRNASeqToImport <-c("rat_Rno20delDup_Del_Wt","rat_Rno20delDup_Dup_Wt","rat_Rno20delDup_DelDup_Wt","rat_DupCBS","rat_Dyrk1a_cko","rat_Rno11Rno20_Rno11Dup_Wt","rat_Rno11Rno20_Rno20Dup_Wt","rat_Rno11Rno20_DoubleDup_Wt", "rat_Rno20Dup")
nameNew <-  c("rat_Rno20_Del_Wt","rat_Rno20_Dup_Wt","rat_Rno20_DelDup_Wt","rat_DupCBS","rat_Dyrk1a_cko","rat_Rno11Dup_Wt","rat_Rno20Dup_Wt","rat_Rno11Dup_Rno20Dup", "rat_Rno20Dup")

import_fcrosFC_EGs.function <- function(nameEGsRNASeqToImport,nameNew,f_inputEGs){
      ###############################################################
      # NOTE: this function will help to import all the EGs data coming 
      # from different experiment (RNASeq/microarrays) if all the models
      # per technology are stored in the same directories (one for rnaseq, 
      # another for the arrays)
      # output: 2 list files will be created, 
      #  1. EGsListReg: list  containing the EGs name plus a
      #     second column with the regulation
      #  2. EGsList: EGs input ready to be feed to venn function
      ###############################################################
	EGsList <-list()
 
    for (i in 1:length(nameEGsRNASeqToImport)){
            
        if (nameEGsRNASeqToImport[i]=="TgDyrk1a"){
              nameNew <- "TgDyrk1a_E15.5"
        } else if (nameEGsRNASeqToImport[i]=="G8Dlx"){
              nameNew <- "G8Dlx_E15.5"
        } else if (nameEGsRNASeqToImport[i]=="G8Het"){
              nameNew <- "Dyrk1a_het_E15.5_seRNASeq"
        } else if (nameEGsRNASeqToImport[i]=="G7Dlx"){
              nameNew <- "G7Dlx_E15.5"
        } else{
              nameNew <- nameEGsRNASeqToImport[i]
    	}
	    EGs <- local({load(paste0(wd,f_inputEGs,nameEGsRNASeqToImport[i],"_fcros_EGs_results.RData")) 
	    	stopifnot(length(ls())==1) 
	        environment()[[ls()]]
	        });
	    if (i==1){
	        EGsFC <- unique(EGs[,c(1,2,4)]);
	        colnames(EGsFC)[3] <- paste0("FC_",nameNew);
	        EGs$model <- paste0(nameNew,"_seRNASeq");
	        EGsList[[i]]<- EGs;
	        names(EGsList)[i] <- paste0(nameNew,"_seRNASeq");
	    } else {
	        EGsFC.tmp <- unique(EGs[,c(1,2,4)]);
	        colnames(EGsFC.tmp)[3] <- paste0("FC_",nameNew);
	        EGsFC <- full_join(EGsFC,EGsFC.tmp,by=c("EnsembleGeneName", "swissprot"));
	        EGs$model <- paste0(nameNew,"_seRNASeq");
	        EGsList[[i]]<- EGs;
	        names(EGsList)[i] <- paste0(nameNew,"_seRNASeq");
	    }
	    i=i+1
	}
    assign("EGsList_afDf_objects",EGsList,.GlobalEnv)
    assign("EGsFC",EGsFC,.GlobalEnv)

};
import_fcrosFC_EGs.function(nameEGsRNASeqToImport,nameNew,f_inputEGs)


import_deseq2FC_EGs.function <- function(nameEGsRNASeqToImport,nameNew,f_dese2Input,organism){
      ###############################################################
      # NOTE: this function will help to import all the EGs data coming 
      # from different experiment (RNASeq/microarrays) if all the models
      # per technology are stored in the same directories (one for rnaseq, 
      # another for the arrays)
      # output: 2 list files will be created, 
      #  1. EGsListReg: list  containing the EGs name plus a
      #     second column with the regulation
      #  2. EGsList: EGs input ready to be feed to venn function
      ###############################################################
	EGsList <-list()
	ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = organism,mirror="uswest"); 
	attributes = listAttributes(ensembl);
	attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])));
	attrib <-attributes[attrib,1];  

    for (i in 1:length(nameEGsRNASeqToImport)){
            
        if (nameEGsRNASeqToImport[i]=="TgDyrk1a"){
              nameNew <- "TgDyrk1a_E15.5"
        } else if (nameEGsRNASeqToImport[i]=="G8Dlx"){
              nameNew <- "G8Dlx_E15.5"
        } else if (nameEGsRNASeqToImport[i]=="G8Het"){
              nameNew <- "Dyrk1a_het_E15.5_seRNASeq"
        } else if (nameEGsRNASeqToImport[i]=="G7Dlx"){
              nameNew <- "G7Dlx_E15.5"
        } else{
              nameNew <- nameEGsRNASeqToImport[i]
    	}
	    EGs <- local({load(paste0(wd,f_dese2Input,nameEGsRNASeqToImport[i],"_deseq2_EGs_results.RData")) 
	    	stopifnot(length(ls())==1) 
	        environment()[[ls()]]
	        });
	    EGs <- as.data.frame(EGs);
	    EGs$EnsembleGeneName <- rownames(EGs);
	    EGs$FC <- 2^EGs$log2FoldChange;
	    colnames(EGs)[8] <- paste0("FC_",nameNew);
		egsInput <-unique(EGs$EnsembleGeneName);
		egsAnnot = getBM(attributes=attrib,filters="ensembl_gene_id",values=egsInput, mart=ensembl);
		colnames(egsAnnot) <- c("EnsembleGeneName","Gene_symbol"); 
		egsAnnot$EnsembleGeneName <- as.character(egsAnnot$EnsembleGeneName)
		egsAnnotF<- left_join(EGs,egsAnnot,by="EnsembleGeneName");

	    if (i==1){

	        EGsFC <- unique(egsAnnotF[,c(7,9,8,2)]);
	        colnames(EGsFC)[4] <- paste0("log2FC_",nameNew);
	        EGs$model <- paste0(nameNew,"_seRNASeq");
	        EGsList[[i]]<- EGs;
	        names(EGsList)[i] <- paste0(nameNew[i],"_seRNASeq");
	    } else {
	    	print(i);
	        EGsFC.tmp <- unique(egsAnnotF[,c(7,9,8,2)]);
	        colnames(EGsFC.tmp)[4] <- paste0("log2FC_",nameNew);
	        EGsFC <- full_join(EGsFC,EGsFC.tmp,by=c("EnsembleGeneName", "Gene_symbol"));
	        EGs$model <- paste0(nameNew,"_seRNASeq");
	        EGsList[[i]]<- EGs;
	        names(EGsList)[i] <- paste0(nameNew,"_seRNASeq");
	    }
	    i=i+1
	}
    assign("EGsList_res.Deseq_objects",EGsList,.GlobalEnv)
    colnames(EGsFC)[2] <- "swissprot"
    assign("EGsFC.deseq2",EGsFC,.GlobalEnv)

};
organism<-"rnorvegicus_gene_ensembl"
import_deseq2FC_EGs.function(nameEGsRNASeqToImport,nameNew,f_dese2Input,organism)

inputFiles_for_FC_plot.function<- function(EGsFC,duplicatedRegion,organism,modelName,datePlotFC,typeFC){
	# Add the duplicated region in the form of chr:beginingPos:EndPos, or a vector with several regions like:
	# mouseHsa21selFC<-c("10:76207174:78464975", "16:75548585:97967873", "17:30953888:32065118");

	f_biocatInput <- "/Users/marmmoreno/Documents/YH_Lab/collabs_outside/DS_models_book/results/RData/"	
	load(file=paste0(f_biocatInput, "/biotypeCategory_pertenenceData.RData") ) #biotypeCat

	nameDupArea <- paste0("chr",gsub(":.*","",duplicatedRegion))
	ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = organism,mirror="useast"); 
	attributes = listAttributes(ensembl);
	
	# A. Hsa21 syntenic genes: Creating input files #######   
	#######################################################
	#
	# rats:
	#########
	# rno11Region <- "11:Lipi:Zbtb21" 11:13960230:38457380
	#### lipi: 11:13960246:13999862
	#### Zbtb21: 11:38442742:38457373
	#
	#mouse:
	#########
	# mouseHsa21selFC<-c("10:76207174:78464975", "16:75548585:97967873", "17:30953888:32065118");
	# IMPORTANT NOTE:
	# I want only the annotations of genes that are inside the synteny regions of human chr21 in mouse/rat. For mouse: 
	# these regions are defined in ENSEMBL webpage , but modifying the info of the mm10 to include the info Li et al 
	# 2010, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2893810/
	# The end of mm10 synteny region is supposed to be after Pdxk after 78464975
	# regions<-c("10:76207174:78307744", "16:75548585:97967873", "17:30953888:32065118") #ensembl defined

	#attributes interesting: 
	attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='start_position'),][1,])), 
		as.numeric(rownames(attributes[which(attributes[,1]=='end_position'),][1,])), 
		as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,]))
		);
	attrib <-attributes[attrib,1];  
	dupRegion = unique(getBM(attributes=attrib,filters="chromosomal_region",values=duplicatedRegion, mart=ensembl));
	colnames(dupRegion) <- c("EnsembleGeneName","Gene_symbol","chr","start","end","gene_biotype"); 
 	dupRegion$chr<- as.character(dupRegion$chr);
 	dupRegion <- left_join(dupRegion,biotypeCat,by="gene_biotype");
 	dupRegion_full <-dupRegion

 	statsDupRegion <- as.data.frame(dupRegion[,c('Gene_symbol','chr','gene_biotype','biotype_category')] %>% group_by(chr,biotype_category)  %>% count());
	stats2<- data.frame(syntenic_genes=length(unique(dupRegion$Gene_symbol)));
	stats2$syntenic_isoforms <-length(unique(dupRegion$EnsembleGeneName));	
	listRes0 <- list("dupRegion" = dupRegion, "stats_DupRegion" = statsDupRegion,stats2=stats2);

	#if the region  have several chrs, to know how many genes per chr:
	chrHSA21 <- unique(dupRegion$chr);
	ensembl_chr <- list();
	for (i in chrHSA21){
	ensembl_chr[[i]] <- dupRegion[which(dupRegion$chr == i), ];
	};  
	print(paste0("Number of HSA21 syntenic genes defined by Ensembl: ",sapply(ensembl_chr, nrow)));

	#saving these trisomic genes in a table
 	assign(paste0(nameDupArea,".","dupRegion_full","_",typeFC),dupRegion,.GlobalEnv);
 	save(dupRegion, file=paste0(wd,f_results, f_FC,f_hsa21, "Ensembl_Annotation_HSA21_",datePlotFC,"_",typeFC,".RData"));
	write.xlsx(listRes0, file = paste0(wd,f_results, f_FC,f_hsa21, "Ensembl_Annotation_HSA21_syntenic_region_",modelName,"_",nameDupArea,"_",datePlotFC,"_",typeFC,".xlsx"));
	#write.xlsx(ensembl_chr, file = paste0(wd,f_results, f_FC,f_hsa21, "Ensembl_Annotation_perChr_HSA21_",datePlotFC,".xlsx")); #if more than one chr region is used 
	
	dupRegion<- unique(dupRegion[,c(2,3,4)]);
 	assign(paste0(nameDupArea,".","dupRegion_short","_",typeFC),dupRegion,.GlobalEnv);

	# B. annotating all EGs coming from  
	#    fcros analysis with chr location  
	#    and start/end: 			#######   
	#######################################
	egsInput <-unique(EGsFC$EnsembleGeneName);

	egsAnnot = getBM(attributes=attrib,filters="ensembl_gene_id",values=egsInput, mart=ensembl);
	colnames(egsAnnot) <- c("EnsembleGeneName","Gene_symbol","chr","start","end","gene_biotype"); 
 	egsAnnot$chr<- as.character(egsAnnot$chr);
	print(paste0("Number of genes EGs: ",length(unique(egsAnnot$Gene_symbol)),"; Number of isoforms EGs: ",length(unique(egsAnnot$EnsembleGeneName))));
	egsAnnotF<- left_join(EGsFC,egsAnnot,by="EnsembleGeneName");
	egsAnnotF<- left_join(egsAnnotF,biotypeCat,by="gene_biotype");
	#ordering: by start first, and by numeric order of the chr second
	egsAnnotF <- egsAnnotF[order(egsAnnotF$start,decreasing=FALSE),];
	egsAnnotF <- egsAnnotF[mixedorder(egsAnnotF$chr,decreasing=FALSE),];
	rownames(egsAnnotF) <- 1:nrow(egsAnnotF); #reorder the rownames, to use as gene init/end coordinates for the plot

	#per chr
	chrEGs <- unique(egsAnnot$chr);
	ensembl_chr_EGs <- list();
	for (i in chrEGs){
	ensembl_chr_EGs[[i]] <- egsAnnot[which(egsAnnot$chr == i), ];
	};  

	#stats by gene not by isoform
	statsEGs <- as.data.frame(egsAnnotF[,c('Gene_symbol','chr','gene_biotype','biotype_category')] %>% group_by(chr)  %>% count());
	statsEGs <- statsEGs[mixedorder(statsEGs$chr,decreasing=FALSE),];
	statsChrEgs <- as.data.frame(egsAnnotF[,c('Gene_symbol','chr','gene_biotype','biotype_category')] %>% group_by(chr,biotype_category)  %>% count());
	statsChrEgs <- statsChrEgs[mixedorder(statsChrEgs$chr,decreasing=FALSE),];

	statsChrEgsBioType <- as.data.frame(egsAnnotF[,c('Gene_symbol','chr','gene_biotype','biotype_category')] %>% group_by(chr,biotype_category,gene_biotype)  %>% count());
	statsChrEgsBioType <- statsChrEgsBioType[mixedorder(statsChrEgsBioType$chr,decreasing=FALSE),];

	#saving these trisomic genes in a table
	save(egsAnnotF, file=paste0(wd,f_results, f_FC, "Ensembl_Annotation_EGs_",datePlotFC,".RData"));
 	listRes <- list("ensembl_chr" = ensembl_chr, "dupRegion" = dupRegion);
	write.xlsx(listRes, file = paste0(wd,f_results, f_FC,f_hsa21, "Ensembl_Annotation_perChr_HSA21_",datePlotFC,".xlsx"));
	write.xlsx(egsAnnotF, file = paste0(wd,f_results, f_FC,f_hsa21, "Ensembl_Annotation_fcros_EGs_",datePlotFC,".xlsx"));
	statsEGsList<- list(statsEGs=statsEGs,statsChrEgs=statsChrEgs,statsChrEgsBioType=statsChrEgsBioType);
	write.xlsx(statsEGsList, file = paste0(wd,f_results, f_FC,f_hsa21, "Ensembl_Annotation_fcros_EGs_STATS_",datePlotFC,".xlsx"));
 	assign(paste0("egsAnnot_full","_",typeFC),egsAnnotF,.GlobalEnv);
 	assign("biotypeCat",biotypeCat,.GlobalEnv);
 	assign(paste0("ensembl_chr_EGs","_",typeFC),ensembl_chr_EGs,.GlobalEnv);
	assign(paste0("statsEGsList","_",typeFC),statsEGsList,.GlobalEnv);
 	assign(paste0(nameDupArea,".","DupRegion_listRes0","_",typeFC),listRes0,.GlobalEnv);
	
	# C.  EGs: get 3 
	#    border genes to the HSA21 
	#    syntenic region in the specie
	#######################################

	#positions Yann wants: 3 genes before and after the pointbreaks
	posInit <-  as.numeric(rownames(egsAnnotF[which(egsAnnotF$Gene_symbol==dupRegion_full[1,2]),]))-3;    
	posEnd <-  as.numeric(rownames(egsAnnotF[which(egsAnnotF$Gene_symbol==dupRegion_full[nrow(dupRegion_full),2]),]))+3;      

	Hsa21PlotData <- egsAnnotF[c(posInit:posEnd),]; #data to plot
	rownames(Hsa21PlotData) <- 1:nrow(Hsa21PlotData); #reorder the rownames, to use as gene init/end coordinates for the plot
	
	# To know how many genes we have FC annotated in at least one model 
	# and how many we did not have any expression info:
	#########################################################################
	print(paste0("The number of genes we have FC info in at least 1 model is : ",dim(Hsa21PlotData)[1] ));#66
	print(paste0("The number of  genes inside the duplicated region : ",length(unique(dupRegion_full[,2])))); #382

	# defining star and end of hsa21 genes
	#########################################################################
	dupInit	 <- as.numeric(rownames(Hsa21PlotData[which(Hsa21PlotData$Gene_symbol==dupRegion_full[1,2]),]));
	DupEnd <- as.numeric(rownames(Hsa21PlotData[which(Hsa21PlotData$Gene_symbol==dupRegion_full[nrow(dupRegion_full),2]),]));

	#Tables for the paper ####
	#############################
	#############################
	FC_objectsList<- list(ratModelHsa21genes=dupRegion_full, hsa21_genes_plus_borders=Hsa21PlotData);
	save(dupRegion_full,Hsa21PlotData,egsAnnotF, file=paste0(wd,f_results, f_FC,f_hsa21, "FC_objects_generationInputFiles_",modelName,"_",nameDupArea,"_",datePlotFC,"_",typeFC,".Rdata"));
	write.xlsx(FC_objectsList, file=paste0(wd,f_results, f_FC,f_hsa21, "TablesPaper_chrs_genes_FC_Hsa21Plot",modelName,"_",nameDupArea,"_",datePlotFC,"_",typeFC,".xlsx"), col.names=T,row.names=F, sep="\t",quote=FALSE );
 	assign(paste0(nameDupArea,".","dupInit","_",typeFC),dupInit,.GlobalEnv);
 	assign(paste0(nameDupArea,".","DupEnd","_",typeFC),DupEnd,.GlobalEnv);
 	assign(paste0(nameDupArea,".","Hsa21PlotData","_",typeFC),Hsa21PlotData,.GlobalEnv);
};

#######################################################

#NOTE:
#for the FC calculated using DESeq2, 
#the input file has the same format 
#than for fcros, so we just use the same function
# Just specify the typeFC to DESeq2 instead of fcros
#######################################################

typeFC<- "fcros"
duplicatedRegion <- "11:13960230:38457380";
modelName <- "Rno11"
inputFiles_for_FC_plot.function(EGsFC,duplicatedRegion,organism,modelName,datePlotFC,typeFC)

#
typeFC<- "DESeq2"
duplicatedRegion <- "11:13960230:38457380";
modelName <- "Rno11"
inputFiles_for_FC_plot.function(EGsFC.deseq2,duplicatedRegion,organism,modelName,datePlotFC,typeFC)

# "Number of HSA21 syntenic genes defined by Ensembl: 211"
# "Number of genes EGs: 23184; Number of isoforms EGs: 24835"       
# "The number of genes we have FC info in at least 1 model is : 160"
# "The number of  genes inside the duplicated region : 194"

################################################################################33
################################################################################33
#
#B) ploting the genes in HSA21 region
#
################################################################################33
################################################################################33

df1<- chr11.Hsa21PlotData_fcros
chr11.dupInit_pos <- as.numeric(chr11.dupInit_fcros) - 0.5 #before this gene
chr11.DupEnd_pos <- as.numeric(chr11.DupEnd_fcros) + 0.5 #after this gene
colnames(df1) <- gsub(" Wt","/Wt",gsub("_"," ",gsub("FC_rat_","",colnames(df1))))
positionLegend <-"topright"
nameSavePlot <- "Hsa21_chr11_fcros"
namePlot <-"Rno6 chr11"
cexNb <- 12
cexGeneNames=0.41
cexNb1<- 0.56

plot_FC_chr11.function <- function(df1,chr11.dupInit_pos,chr11.DupEnd_pos,positionLegend,namePlot, nameSavePlot,cexNb,cexNb1,cexGeneNames){
	
	pdf(file=paste0(wd,f_results, f_FC,f_Figure,"/",nameSavePlot,".pdf"), width=11.69, height=8.27)
	op <- par(mfrow = c(7,1),
	          oma = c(5.4,1.3,1.4,0) + 0.1,
	          mar = c(0,0.7,0.7,0.7) + 0.1,
	          cex=1.2)

	#1 Rno20delDup Del/Wt
	plot(factor(df1[,2], levels=df1[,2]), df1[,3], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,3], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	text(chr11.dupInit_pos+0.5, 3.8,labels=paste0("                            ",gsub("delDup","",gsub(" .*","",colnames(df1)[3]))), col = "lightpink4",cex=0.9, font=2)
	legend(positionLegend, legend=gsub(".* ","",colnames(df1)[3]), bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03)


	#2 Rno20delDup Dup/Wt
	plot(factor(df1[,2], levels=df1[,2]), df1[,4], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,4], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.5, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	text(chr11.dupInit_pos+0.5, 3.8, labels=paste0("                           ",gsub("delDup","",gsub(" .*","",colnames(df1)[4]))), col = "darkolivegreen4",cex=0.9, font=2)
	legend(positionLegend, legend=gsub(".* ","",colnames(df1)[4]), bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03)

	#3 Rno20delDup DelDup/Wt
	plot(factor(df1[,2], levels=df1[,2]), df1[,5], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,5], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.5, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	text(chr11.dupInit_pos+0.5, 3.8,labels=paste0("                          ",gsub("delDup","",gsub(" .*","",colnames(df1)[5]))), col = "cornflowerblue",cex=0.9, font=2)
	legend(positionLegend, legend=gsub(".* ","",colnames(df1)[5]), bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03)
	mtext("FC", side=2, line=1.4, cex=0.8)


	#4 Rno11Rno20 Rno20Dup/Wt
	plot(factor(df1[,2], levels=df1[,2]), df1[,9], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	panel.first = rect(chr11.dupInit_pos, -1e6, chr11.DupEnd_pos, 1e6, col=rgb(1,0.713,0.756, alpha=0.37), border=NA, lwd=2) #Rno20DeldUP
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,9], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.5, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	text(chr11.dupInit_pos+0.5, 3.8, labels=paste0("                           ",gsub("delDup","",gsub(" .*","",colnames(df1)[9]))), col = "olivedrab4",cex=0.9, font=2)
	legend(positionLegend, legend=gsub(".* ","",colnames(df1)[9]), bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03)


	#5 Rno11Rno20 DoubleDup/Wt
	plot(factor(df1[,2], levels=df1[,2]), df1[,10], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	panel.first = rect(chr11.dupInit_pos, -1e6, chr11.DupEnd_pos, 1e6, col=rgb(1,0.713,0.756, alpha=0.37), border=NA, lwd=2) #Rno20DeldUP
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,10], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.5, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	text(chr11.dupInit_pos+0.5, 3.8, labels=paste0("                           ",gsub(" .*","",colnames(df1)[10])), col = "springgreen4",cex=0.9, font=2)
	legend(positionLegend, legend=gsub(".* ","",colnames(df1)[10]), bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03)


	#6 Rno11Rno20 Rno11Dup/Wt
	plot(factor(df1[,2], levels=df1[,2]), df1[,8], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	panel.first = rect(chr11.dupInit_pos, -1e6, chr11.DupEnd_pos, 1e6, col=rgb(1,0.713,0.756, alpha=0.37), border=NA, lwd=2) #Rno20DeldUP
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,8], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.5, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	text(chr11.dupInit_pos+0.5, 3.8, labels=paste0("                           ",gsub("/Wt","",gsub("Rno11Rno20 ","",colnames(df1)[8]))), col = "snow4",cex=0.9, font=2)
	legend(positionLegend, legend=gsub(".* ","",colnames(df1)[8]), bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03);

	#7 Dyrk1a cko
	plot(factor(df1[,2], levels=df1[,2]), df1[,7], col="black", pch=19,  frame=TRUE, cex=cexNb, las=2, xaxt = "n",cex.lab=1,cex.axis=0.5, ylab="FC", xlab="",  ylim=c(0.1,5))
	#panel.first = rect(dyrkcko_pos -0.5, -1e6, dyrkcko_pos+0.5, 1e6, col=rgb(0.913, 0.588, 0.478, alpha=0.47), border=NA) #Dp1Rhr
	par(new=TRUE)
	points(factor(df1[,2], levels=df1[,2]), df1[,7], col="black", pch=19,  cex=cexNb1, las=2,xaxt = "n", cex.lab=1,cex.axis=0.5, ylab="", xlab="",  ylim=c(0.1,5))
	abline(h=1, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=1.4, col= "deeppink",lty = 6,lwd=0.9)
	abline(h=2.5, col= "deeppink",lty = 6,lwd=1.3)
	legend(positionLegend, legend=" Dyrk1a +/-", bty='n', col=rgb(1,1,0.4, alpha=0.002), text.col="black",pch=26, cex=0.75,lwd=6)
	axis(side=1, at=1:length(df1[,2]), tick=TRUE, labels=rep("", times=length(df1[,2])),tck=-0.03)
	op <- par(cex=cexGeneNames,mgp=c(3,0.35,0))
	axis(side=1, at=1:length(df1[,2]), cex=1.5,las=2,tick=TRUE, labels=df1[,2],tck=-0.06)
	title(namePlot, line = -0.5, outer=TRUE, cex.main=3.7)
	par(op)
	dev.off()

}



################################################################################################################################
################################################################################################################################
############################### END.  11/02/2019 ####################################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################
################################################################################################################################
################################################################################################################################


