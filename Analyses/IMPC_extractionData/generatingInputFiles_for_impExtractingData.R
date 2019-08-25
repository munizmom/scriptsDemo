

###################################################################
###################################################################
###################################################################
#
# Script to generate the input files needed to extract 
# genotype-phenotye info from the IMPC webpage. 
#   Two secuential input files are needed.
###################################################################
###################################################################
###################################################################
#
# part1. A dataframe or txt file containing the genes of interest
# we want to extract the phenotype info from the IMPC data. 
# Language: R
# Example: getting the name of the genes included in the mouse 
# Hsa21 synteny region
#
###################################################################

############# Packages needed ###############
source("http://bioconductor.org/biocLite.R");
library("tidyr");
library("dplyr");
library("biomaRt");
library("xlsx");
#############################################

################## CODE #####################


GeneSymbol_inside_GenomicRegion.function <- function(genomicRegions,organismNameEnsemblDataset,nameData,f_output, date){

	#README: this function will generate two files:
	###################################################################
	a) Ensembl_Annotation_ : from this file open in excel & copy as 
	txt only the gene Names :) That is the real input file to be feed 
	to the IMPC web.
	b) Ensembl_Annotation_perChr
	###################################################################
	
	ensembl_mo = useMart("ensembl", dataset = organismNameEnsemblDataset); 
	attributes_m = listAttributes(ensembl_mo);

	#interesting attributes start_position
	attrib <-c(as.numeric(rownames(attributes_m[which(attributes_m[,1]=='ensembl_gene_id'),][1,])),
		as.numeric(rownames(attributes_m[which(attributes_m[,1]=='external_gene_name'),][1,])),
		as.numeric(rownames(attributes_m[which(attributes_m[,1]=='chromosome_name'),][1,])),
		as.numeric(rownames(attributes_m[which(attributes_m[,1]=='start_position'),][1,])), 
		as.numeric(rownames(attributes_m[which(attributes_m[,1]=='end_position'),][1,])), 
		as.numeric(rownames(attributes_m[which(attributes_m[,1]=='gene_biotype'),][1,]))
			);

	attrib <-attributes_m[attrib,1]; # 
	GenesInfo = getBM(attributes=attrib,filters="chromosomal_region",values=genomicRegions, mart=ensembl_mo);
	colnames(GenesInfo) <- c("EnsembleGeneName","Gene_symbol","chr","start","end","gene_biotype"); 
 	GenesInfo$chr<- as.character(GenesInfo$chr);
	print(paste0("Number of ",nameData,": ",length(unique(GenesInfo$Gene_symbol))));

	chrGenesInfo <- unique(GenesInfo$chr);
	ensembl_chr <- list();

	for (i in chrGenesInfo){
	ensembl_chr[[i]] <- GenesInfo[which(GenesInfo$chr == i), ]
	};  

	print(paste0("Number of ",nameData," : ",sapply(ensembl_chr, nrow)));
	
	#chr16 -->263 -->271 -->384-->386
	#chr17 -->29
	#chr10 --> 84 -->98
	#saving these trisomic genes in a table

 	assign("GenesInfo_full",GenesInfo,.GlobalEnv);
 	save(GenesInfo, file=paste0(f_output, "Ensembl_Annotation_",nameData,date,".RData"))
 	#from this file, open in excel & copy as txt only the gene Names :) Thatis the real input file to be feed to the IMPC web.

	write.xlsx(ensembl_chr, file = paste0(f_output,  "Ensembl_Annotation_perChr_",nameData,date,".xlsx"))
	write.xlsx(GenesInfo, file=paste0(f_output,  "Ensembl_Annotation_",nameData,date,".xlsx"))
}	

#running the function for the example: mouse Hsa21 synteny region
organismNameEnsemblDataset <- "mmusculus_gene_ensembl"
genomicRegions<-c("10:76207174:78464975", "16:75548585:97967873", "17:30953888:32065118"); #mod

f_output<- ('/Users/marmmoreno/Desktop/onWorking/FC_plot/Hsa21/')
GeneSymbol_inside_GenomicRegion.function(genomicRegions,"mmusculus_gene_ensembl","hsa21IMPC1",f_output, "190319")

#####################################################################################################
############################### END. Finished on the  250319 ########################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################

###################################################################
#
# Part2. Downloading the genotype-phenotype data from the IMPC 
# using the API:   
# Language: bash
#
###################################################################
wd='/Users/marmmoreno/Documents/DS_models_book/ongoing/results/IMPC/dataRetrieval/'
#example:
wget 'http://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:Akt2&wt=json' -O  Akt2_ex.txt

#For the mouse synteny Hsa21 region including genes: #####
##########################################################
#Input list of hsa21 genes: /Users/marmmoreno/Desktop/onWorking/FC/Hsa21/Hsa21Genes_190319.txt
#Input list with the urls to download:  awk '{OFS="\t"} {print "http://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:" $1 "&wt=json";}' > IMPCurlList.txt

#code to run:
wd='/Users/marmmoreno/Documents/DS_models_book/ongoing/results/IMPC/dataRetrieval/'
cd $wd
awk  '{OFS="\t"} {print "http://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:" $1 "&wt=json";}' /Users/marmmoreno/Desktop/onWorking/FC/Hsa21/Hsa21Genes_190319.txt | xargs -L1  wget -O - >>hsa21IMPC.txt

gsed -e '/numFound":0/{n;N;d}' hsa21IMPC.txt > hsa21IMPC_formatted.tmp;
sed -s hsa21IMPC_formatted.tmp <<< $'g/numFound":0/-6,.d\n,p' > hsa21IMPC_formatted.txt;

#####################################################################################################
############################### END. Finished on the  250319 ########################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################



