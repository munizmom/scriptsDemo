###########################################################################################
###########################################################################################
##  Script for the DS models book: lessons learnt from the Genomic data 	             ##
##   																				     ##
###########################################################################################
###########################################################################################
##README
# Author: Mar Muniz. YHerault team. Last update 300319
# This scrip have a serie of functions that will allow to re produce/redo the genome wide
# analysis of  prot coding genes and other regulatory features in the different species
# as described in the book chapter5.
########################################################################
# The reference to the book: 
#  Muniz Moreno M.M., Brault V., Birling M.C., Pavlovic G. 
#  Progress in Brain Research Series: Down Syndrome animal models. 
#  Chapter 5. “Modelling Down Syndrome in animals from the early stage 
#			   to the 4.0 models and next”
########################################################################
###########################################################################################
# NOTES:
###########################################################################################
# Be sure that you are working with an up to date version of R and bioconductor environment,
# update R --> https://www.r-bloggers.com/updating-r-on-ubuntu/
# The hidden parts of the code are there for QC purposes or to take a look at the data
############################################################################################

############################################################################################
##################################### Packages needed ######################################
############################################################################################

#update R https://www.r-bloggers.com/updating-r-on-ubuntu/

############################################################################################
############# Packages needed ###############
#in case of errors:
######################################
#unlink("/home/me/src/Rlibs/00LOCK-Rcpp", recursive = TRUE)
#install.packages("lubridate", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#
#   "gplots",
#   repos = c("http://rstudio.org/_packages",
#   "http://cran.rstudio.com", dependencies = TRUE)
#)
######################################
source("http://bioconductor.org/biocLite.R");

library("biomaRt");
library("tidyr");
library("dplyr");
library("RColorBrewer");
library("gtable");
library("ggplot2");
library("gplots");
library("openxlsx");
library("gridExtra");
library("cowplot");
library("grid");
library("ggpubr");
#library("data.table")
#library("mygene")
#library("stringr")
##################################################################################
#Before begining:
set.seed(22) # need to set it to get always the same random results and plots
#sessionInfo()
wd <- '/Users/marmmoreno/Documents/YH_Lab/collabs_outside/DS_models_book/'
setwd(wd)
wd <- getwd();

#set up output directories
##########################

f_results <- "/results/";
f_Rdata <- "RData/";
f_tables <- "tables/";
f_figures <-"figures/";
f_homologs <- "homologs/"
f_stats <-"stats/"
f_trisomic <- "trisomicModels/"
load(file=paste0(wd, f_results, f_Rdata,"session.RData"))

# Creating the folders
##################################################################################
dir.create(file.path(getwd (), f_results), showWarnings = F);
dir.create(file.path(getwd (), f_results, f_Rdata), showWarnings = F);
dir.create(file.path(getwd (), f_results, f_tables), showWarnings = F);
dir.create(file.path(getwd (), f_results, f_figures), showWarnings = F);
dir.create(file.path(getwd (), f_results,f_homologs), showWarnings = F);
dir.create(file.path(getwd (), f_results,f_homologs, f_Rdata), showWarnings = F);
dir.create(file.path(getwd (), f_results,f_homologs, f_tables), showWarnings = F);
dir.create(file.path(getwd (), f_results,f_homologs, f_figures), showWarnings = F);
dir.create(file.path(getwd (), f_results,f_trisomic), showWarnings = F);

##################################################################################
#DEFINING BIOTYPES CATEGORIES in a df
#ProtCot: IGC gene, IGD gene, IG gene, IGJ gene, IGLV gene, IGM gene, IGV gene, IGZ gene, nonsense mediated decay, nontranslating CDS, non stop decay, polymorphic pseudogene, TRC gene, TRD gene, TRJ gene.
#Pseudogene: disrupted domain, IGC pseudogene, IGJ pseudogene, IG pseudogene, IGV pseudogene, processed pseudogene, transcribed processed pseudogene, transcribed unitary pseudogene, transcribed unprocessed pseudogene, translated processed pseudogene, TRJ pseudogene, unprocessed pseudogene
#Long noncoding: 3prime overlapping ncrna, ambiguous orf, antisense, antisense RNA, lincRNA, ncrna host, processed transcript, sense intronic, sense overlapping
#Short noncoding: miRNA, miRNA_pseudogene, miscRNA, miscRNA pseudogene, Mt rRNA, Mt tRNA, rRNA, scRNA, snlRNA, snoRNA, snRNA, tRNA, tRNA_pseudogene


#to know the numb and name of chrchecking both the biomaRt info and the ensembl caryotipes ex. 
# http://www.ensembl.org/Danio_rerio/Location/Genome
speciesBiomaRt <- c("hsapiens_gene_ensembl","mmusculus_gene_ensembl","rnorvegicus_gene_ensembl",
	"dmelanogaster_gene_ensembl","ptroglodytes_gene_ensembl","ggorilla_gene_ensembl",
	"pabelii_gene_ensembl","drerio_gene_ensembl","celegans_gene_ensembl");

speciesNames <- c("Human","Mouse","Rat","Drosophila","Chimpanzee","Gorilla",
	"Orangutan","Zebrafish","C.elegands");

chrsSpecieList <- list()
genesDf_all <- list()
for (i in 1:length(speciesBiomaRt)){
	ensembl = useMart("ensembl", dataset = speciesBiomaRt[i]) 
	attributes = listAttributes(ensembl)
	attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='start_position'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,])))
	attrib <-attributes[attrib,1] # "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "gene_biotype"  
	genesDf_all = getBM(attributes=attrib,mart=ensembl); 
	#print(unique(genesDf_all$chromosome_name));
	chrsSpecieList[[i]]<- unique(genesDf_all$chromosome_name);
	names(chrsSpecieList)[i] <- speciesBiomaRt[i];
	genesDf_allCuack[[i]] <- genesDf_all
	names(genesDf_allCuack)[i] <- speciesBiomaRt[i];
	i=i+1
};	

chrHuman <- c(rep(1:22, times=1),"X","Y","MT");
chrMouse <- c(rep(1:19, times=1),"X","Y","MT");
chrRattu <- c(rep(1:19, times=1),"X","Y","MT");
chrDroso <- c("2R","2L","3R","3L","4","X","Y","mitochondrion_genome");
chrChimp <- c("1","2A","2B",rep(3:22, times=1),"X","Y","MT");
chrGoril <- c("1","2A","2B",rep(3:22, times=1),"X","MT");
chrOrang <- c("1","2a","2b",rep(3:22, times=1),"X","MT");
chrZebrf <- c(rep(1:25, times=1),"MT");
chrCeleg <- c("I", "IV", "V", "X", "II", "III", "MtDNA");

##################################################################################
#A. GENOME WIDE ANALYSIS
##################################################################################
speciesBiomaRt <- c("hsapiens_gene_ensembl","mmusculus_gene_ensembl","rnorvegicus_gene_ensembl",
	"dmelanogaster_gene_ensembl","ptroglodytes_gene_ensembl","ggorilla_gene_ensembl",
	"pabelii_gene_ensembl","drerio_gene_ensembl","celegans_gene_ensembl");

speciesNames <- c("Human","Mouse","Rat","Drosophila","Chimpanzee","Gorilla",
	"Orangutan","Zebrafish","C.elegands");

chrList <- list(chrHuman=chrHuman,chrMouse=chrMouse,chrRattu=chrRattu,chrDroso=chrDroso,chrChimp=chrChimp,chrGoril=chrGoril,chrOrang=chrOrang,chrZebrf=chrZebrf,chrCeleg=chrCeleg);
ensemblLengthList_270219<- list(chrHuman=3609003417,chrMouse=3486944526,chrRattu=3042335753,chrDroso=142573024,chrChimp=3385800935,chrGoril=2917385452,chrOrang=3109347532,chrZebrf=1674207132,chrCeleg=103022290);

# #celegans
# Genebuild released	Oct 2014
# Genebuild last updated/patched	Jun 2017
# WBcel235, INSDC Assembly GCA_000002985.3, Dec 2012
# #orangutan
# Genebuild released	Mar 2008
# Genebuild last updated/patched	Aug 2012
# PPYG2, Sep 2007
# #mus
# Genebuild released	Jul 2012
# Genebuild last updated/patched	Sep 2018
# GRCm38.p6 (Genome Reference Consortium Mouse Reference 38), INSDC Assembly GCA_000001635.8, Jan 2012
# #human
# Genebuild released	Jul 2014
# Genebuild last updated/patched	Jul 2018
# GRCh38.p12 (Genome Reference Consortium Human Build 38), INSDC Assembly GCA_000001405.27, Dec 2013
# #gorilla
# Genebuild released	Dec 2017
# Genebuild last updated/patched	Jan 2018
# gorGor4, INSDC Assembly GCA_000151905.3, Dec 2014
#chimp
# assembly Pan_tro_3.0, INSDC Assembly GCA_000001515.5, May 2016
# Genebuild released	Dec 2017
# Genebuild last updated/patched	Jan 2018
# #Droso
# Genebuild by	FlyBase
# BDGP6.22, INSDC Assembly GCA_000001215.4
# #zebraf
# GRCz11 (Genome Reference Consortium Zebrafish Build 11), INSDC Assembly GCA_000002035.4, May 2017
# Genebuild released	Mar 2018
# Genebuild last updated/patched	Apr 2018
# rat
# Rnor_6.0, INSDC Assembly GCA_000001895.4, Jul 2014
# Genebuild released	Jun 2015
# Genebuild last updated/patched	Jan 2017


genomeWide_geneAnnot.function <- function(speciesBiomaRt,speciesNames,chrList,ensemblLengthList){
	##################################################################################
	# Function to annotate all the genes of the genomes of interest with prot coding/non coding annot
	# and gives the stats genome wide. Input files are:
	###1. The ensembl biomaRt species name
	###2. The names we will use for annotation.In the same order of species than the input file nb 1 
	###3. A list were each element  is the vector of the chromosomes each specie has.
	###   In the same order of species than the other 2 input files 
	##################################################################################
	specieResList <- list();
	for (i in 1:length(speciesBiomaRt)){
		
		f_specie <- paste0(gsub(" ","_",speciesNames[i]),"/");
		name <- speciesNames[i];
		chrSel <- chrList[[i]];
		specieAnalysis <- speciesBiomaRt[i];

		#1. Creating the folders
		######################
		dir.create(file.path(getwd (), f_results,f_specie, f_figures), showWarnings = F);
		dir.create(file.path(getwd (), f_results,f_specie, f_tables), showWarnings = F);
		dir.create(file.path(getwd (), f_results,f_specie, f_Rdata), showWarnings = F);
		
		#2. Quering ENSEMBL 
		######################
		#ensembl <- useEnsembl(biomart = "ensembl")
		#listDatasets(mart = ensembl)
		
		#ensembl = useMart("ensembl", dataset = specieAnalysis) 
		ensembl = useMart("ensembl", dataset = speciesBiomaRt[i]); 

		attributes = listAttributes(ensembl);
		attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),
			as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='start_position'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='end_position'),][1,])),
			as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,])));
		attrib <-attributes[attrib,1]; # "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "gene_biotype"  
		genesDf_all = getBM(attributes=attrib,mart=ensembl); #64914 nb of genes including isoforms and genes in scaffolds
		genesDf_isof <- genesDf_all[which(genesDf_all$chromosome_name %in% chrList[[i]]), ]; #58676 nb of genes including isoforms
		genesDf_uniq <- genesDf_isof[!duplicated(genesDf_isof$external_gene_name), ]; #57152 #nb of uniquely genes, excluding isoforms
		#genesDf[which(genesDf$external_gene_name =="IGF2"), ]  
		genesDf_uniq$model <- speciesNames[i];
		specieResList[[i]] <- genesDf_isof;
		names(specieResList)[i] <- speciesNames[i];
	
		if(i==1){
			AllGenesDf <- genesDf_isof;
		} else {
			AllGenesDf <- rbind(AllGenesDf,genesDf_isof);
			};
		i+i+1
	};
	assign("AllGenesDf",AllGenesDf,.GlobalEnv);

	#2 Creating the biotype categories
	# Classification following ENSEMBL and GENCODE: using the info contained in 
	# https://www.ensembl.org/Help/Faq?id=468, 
	# http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
	# https://www.ensembl.org/info/genome/genebuild/ncrna.html
	# https://www.gencodegenes.org/pages/biotypes.html

	biotypes <- unique(AllGenesDf$gene_biotype);
	biotype.coding <-data.frame(gene_biotype=unique(biotypes[c(grep("_gene",biotypes),grep("protein_coding",biotypes), grep("polymorphic_pseudogene",biotypes),grep("nonsense_mediated_decay",biotypes),
	grep("nontranslating_CDS",biotypes),grep("non_stop_decay",biotypes))]));

	biotype.pseudo <-unique(biotypes[c(grep("pseudogene",biotypes), grep("disrupted_domain",biotypes))]);
	biotype.pseudo <-data.frame(gene_biotype=biotype.pseudo[-c(grep("miRNA_pseudogene", biotype.pseudo),grep("polymorphic_pseudogene", biotype.pseudo),
		grep("miscRNA_pseudogene", biotype.pseudo),grep("tRNA_pseudogene", biotype.pseudo))]);

	biotype.lncRNA <- data.frame(gene_biotype=unique(biotypes[c(grep("3prime_overlapping_ncRNA",biotypes), grep("ambiguous_orf",biotypes),grep("antisense",biotypes),
		grep("lincRNA",biotypes),grep("ncRNA_host",biotypes),grep("^lincRNA$",biotypes),grep("processed_transcript",biotypes),
		grep("sense_intronic",biotypes),grep("sense_overlapping",biotypes),grep("retained_intron",biotypes),grep("lncRNA",biotypes),grep("^non_coding$",biotypes))]));

	biotype.shortNC <- data.frame(gene_biotype=unique(biotypes[c(grep("^miRNA$",biotypes), grep("miscRNA",biotypes),grep("Mt_rRNA",biotypes),
		grep("Mt_tRNA",biotypes),grep("^rRNA$",biotypes),grep("scRNA",biotypes),grep("snlRNA",biotypes),
		grep("snoRNA",biotypes),grep("snRNA",biotypes),grep("^tRNA$",biotypes),grep("tRNA_pseudogene",biotypes),
		grep("vaultRNA",biotypes),grep("siRNA",biotypes),grep("piRNA",biotypes),grep("^ncRNA$",biotypes),
		grep("scaRNA",biotypes),grep("ribozyme",biotypes),grep("^misc_RNA$",biotypes),grep("sRNA",biotypes),grep("pre_miRNA",biotypes) )]));

	biotype.TEC <- data.frame(gene_biotype=unique(biotypes[grep("TEC",biotypes)]));
	rest <- biotypes[-(biotypes %in% c(biotype.coding[,1], biotype.pseudo[,1],biotype.lncRNA[,1],biotype.shortNC[,1],biotype.TEC[,1]))]; 
	biotype.coding$biotype_category <-  "Protein coding"
	biotype.pseudo$biotype_category <-  "Pseudogenes"
	biotype.lncRNA$biotype_category <-  "Long noncoding"
	biotype.shortNC$biotype_category <-  "Short noncoding"
	biotype.TEC$biotype_category <-  "TEC"
	biotypeCat <- unique(rbind(biotype.coding, biotype.pseudo,biotype.lncRNA,biotype.shortNC,biotype.TEC));
	assign("rest",rest,.GlobalEnv);
	assign("biotypeCat",biotypeCat,.GlobalEnv);
	assign("specieResList",specieResList,.GlobalEnv);

	#3. Annotating all the species data with the biotype category

	chrmBioCatList <- list();
	chrmBioCatTypesList<- list(); 
	for (i in 1:length(names(specieResList))){
		df <- specieResList[[i]];
		df <- unique(left_join(df,biotypeCat,by="gene_biotype"));
		specieResList[[i]] <- df;
		nameDf <- names(specieResList)[i]


		# Calcualting stats:
		## A. GenomeWide stats
		######################

		if (i==1){
			
			#stats per biotype 
			bioCatTypes <-as.data.frame(df %>% group_by(biotype_category,gene_biotype) %>% count());
			bioCatTypes$perc <- round((bioCatTypes[,3]/sum(bioCatTypes[,3]))*100,digits=2);			#calculating the %
			colnames(bioCatTypes)[c(3,4)] <- c(nameDf, paste0("% ",nameDf));
			bioCatTypes<- bioCatTypes[order(match(bioCatTypes$gene_biotype,unique(biotypeCat[,1]))),];

			#stats per chr
			chr <-as.data.frame(df %>% group_by(chromosome_name) %>% count());
			lengthChr<- as.data.frame(df %>% group_by(chromosome_name) %>%  summarise(lengthLastgenomicElement = max(end_position)));
			chr<- left_join(lengthChr,chr,by="chromosome_name");
			chr$perc <- round((chr[,3]/sum(chr[,3]))*100,digits=2);			#calculating the %
			colnames(chr)[c(2,3,4)] <- c(paste0("position last genomic element ",nameDf),nameDf, paste0("% ",nameDf));
			chr<- chr[order(match(chr$chromosome_name, chrList[[i]])),];	
			
			#stats per biotype category
			bioCat <-as.data.frame(df %>% group_by(biotype_category) %>% count());
			bioCat$perc <- round((bioCat[,2]/sum(bioCat[,2]))*100,digits=2);			#calculating the %
			bioCat$length <- sum(chr[,2]);		
			bioCat$`ENSEMBL Length`	<- ensemblLengthList[[i]];
			colnames(bioCat)[c(2,3,4,5)] <- c(nameDf, paste0("% ",nameDf),paste0("position last genomic element ",nameDf),paste0("ENSEMBL length ",nameDf));
            bioCat<- bioCat[order(match(bioCat$biotype_category,unique(biotypeCat[,2]))),];
		

		} else {
			#stats per biotype 
			bioCatTypes.tmp <-as.data.frame(df %>% group_by(biotype_category,gene_biotype) %>% count());
			bioCatTypes.tmp$perc <- round((bioCatTypes.tmp[,3]/sum(bioCatTypes.tmp[,3]))*100,digits=2);			#calculating the %
			colnames(bioCatTypes.tmp)[c(3,4)] <- c(nameDf, paste0("% ",nameDf ));
			bioCatTypes <- unique(full_join(bioCatTypes,bioCatTypes.tmp, by=c( "biotype_category", "gene_biotype")));

			#stats per chr
			chr.tmp <-as.data.frame(df %>% group_by(chromosome_name) %>% count());
			lengthChr<- as.data.frame(df %>% group_by(chromosome_name) %>%  summarise(lengthLastgenomicElement = max(end_position)))
			chr.tmp<- left_join(lengthChr,chr.tmp,by="chromosome_name")
			chr.tmp$perc <- round((chr.tmp[,3]/sum(chr.tmp[,3]))*100,digits=2)			#calculating the %
			colnames(chr.tmp)[c(2,3,4)] <- c(paste0("position last genomic element ",nameDf),nameDf, paste0("% ",nameDf));
			chr <- unique(full_join(chr,chr.tmp, by="chromosome_name"));
			
			#stats per biotype
			bioCat.tmp <-as.data.frame(df %>% group_by(biotype_category) %>% count());
			bioCat.tmp$perc <- round((bioCat.tmp[,2]/sum(bioCat.tmp[,2]))*100,digits=2);			#calculating the %
			bioCat.tmp$length <- sum(chr.tmp[,2]);		
			bioCat.tmp$`ENSEMBL Length`	<- ensemblLengthList[[i]];
			colnames(bioCat.tmp)[c(2,3,4,5)] <- c(nameDf, paste0("% ",nameDf),paste0("position last genomic element ",nameDf),paste0("ENSEMBL length ",nameDf));
			bioCat <- unique(full_join(bioCat,bioCat.tmp, by="biotype_category"));
		}
		
		## B. Per chr stats:
		######################
		#stats per chr per biotype category
		#print("cuack")
		chrmBioCat<- as.data.frame(df %>% group_by(chromosome_name,biotype_category) %>% count()); #stats per chr per biotype category
		chrmBioCat<- chrmBioCat[order(match(chrmBioCat[,1], chrList[[i]])),];		
		#print("cuack")

		chrmBioCatList[[i]] <- chrmBioCat;
		names(chrmBioCatList)[i] <- nameDf;
		#print("cuack")

		chrmBioCatTypes<- as.data.frame(df %>% group_by(chromosome_name,biotype_category,gene_biotype) %>% count()); 	#stats per chr per biotype
		chrmBioCatTypes<- chrmBioCatTypes[order(match(chrmBioCatTypes[,1], chrList[[i]])),];	
		totalElements <- as.data.frame(chrmBioCatTypes  %>% group_by(chromosome_name,biotype_category) %>%  summarise(totalElements = sum(n)));
		chrmBioCatTypes <- left_join(chrmBioCatTypes,totalElements, by=c("chromosome_name", "biotype_category"));
		chrmBioCatTypes$`% elements per biotype category`<- round((chrmBioCatTypes[,4]/chrmBioCatTypes[,5] )*100,digits=2);
		#print("cuack")

		chrmBioCatTypesList[[i]] <- chrmBioCatTypes;
		names(chrmBioCatTypesList)[i] <- nameDf;
		#print("cuack")

		i=i+1
	};

	#Saving the results
	genomewideStats <- list(Stats_perCategory=bioCat,Stats_perBiotype=bioCatTypes,Stats_perchr=chr);
	write.xlsx(chrmBioCatList, file = paste0(wd,f_results,  "per_model_chrCategory_Stats.xlsx"));
	write.xlsx(chrmBioCatTypesList, file = paste0(wd,f_results, "per_model_chrBiotype_Stats.xlsx"));
	write.xlsx(genomewideStats, file = paste0(wd,f_results, "wholeGenomeStats.xlsx"));



	assign("specieResList",specieResList,.GlobalEnv);
	assign("chrmBioCatList",chrmBioCatList,.GlobalEnv);
	assign("chrmBioCatTypesList",chrmBioCatTypesList,.GlobalEnv);
	assign("bioCat",bioCat,.GlobalEnv);
	assign("bioCatTypes",bioCatTypes,.GlobalEnv);
	assign("biotypeCat",biotypeCat,.GlobalEnv);
	assign("chr",chr,.GlobalEnv);

};

genomeWide_geneAnnot.function(speciesBiomaRt,speciesNames,chrList,ensemblLengthList_270219);
save.image(file=paste0(wd, f_results, f_Rdata,"session.RData"))
 save(biotypeCat,file=paste0(wd, f_results, f_Rdata,"biotypeCategory_pertenenceData.RData") )


##################################################################################
#B. Trisomic genes  analysis of the lab DS models
##################################################################################

chrMouse <- c(rep(1:19, times=1),"X","Y","MT");
chrRattu <- c(rep(1:19, times=1),"X","Y","MT");


speciesBiomaRt <- c("mmusculus_gene_ensembl","mmusculus_gene_ensembl","mmusculus_gene_ensembl",
	"mmusculus_gene_ensembl","mmusculus_gene_ensembl","mmusculus_gene_ensembl",
	"mmusculus_gene_ensembl");

modelsNames <- c("Ts65Dn","Dp1Yey","Dp3Yah","Dp5/Dp1","Dp5Yah","Dp1Rhr","Tg(Dyrk1a)");

reg_Ts65dn <-c("16:84717500:97962621", "17:3076578:9167836")
reg_Dp1yey <-c("16:75540500:97962631")
reg_Dp3yah <-c("16:75755160:84954430")
reg_Dp5yah <-c("16:85173717:92826149")
reg_Dp1rhr <-c("16:93607720:97947423")
reg_Dp5Dp1 <-c("16:85173717:92826149","16:93607720:97947423")
reg_tgDyrk1a <-c("16:94570010-94695520")
modelsTrisomicRegions <- list(reg_Ts65dn=reg_Ts65dn,reg_Dp1yey=reg_Dp1yey,reg_Dp3yah=reg_Dp3yah,reg_Dp5Dp1=reg_Dp5Dp1,reg_Dp5yah=reg_Dp5yah,reg_Dp1rhr=reg_Dp1rhr)

#Tgdyrk1a_Trisomic <-1

trisChr <- c("16","17")
pg1.to.bps<- 978000000 
ensemblLengthList_270219<- list(chrHuman=3609003417,chrMouse=3486944526,chrRattu=3042335753,chrDroso=142573024,chrChimp=3385800935,
	chrGoril=2917385452,chrOrang=3109347532,chrZebrf=1674207132,chrCeleg=103022290);
#chrList <- list(chrHuman=chrHuman,chrMouse=chrMouse,chrRattu=chrRattu,chrDroso=chrDroso,chrChimp=chrChimp,chrGoril=chrGoril,
	#chrOrang=chrOrang,chrZebrf=chrZebrf,chrCeleg=chrCeleg);

#2 Ts65Dn
##################
#Ts65Dn: Trisomic region in ch16. 
#previous gene mir155 finish 84714204
#   begin: included  Mrpl39: 84717576:84735742  (described in https://www.jax.org/research-and-faculty/tools/cytogenetic-and-down-syndrome-models-resource/protocols)
#   end: After Zbtb21 97947435-97962621
#   Plus the translocation breakpoint area in  chr17  (reported Duchon et al,. 2011)
#       start: pisd-ps2 17: 3076578-3084183 
#       end: 6530411M01Rik 17:9147719-9167836



#Annotation of all genes in the trisomic region of Ts65Dn
###########################################################
reg_Ts65dn <-c("16:84717500:97962621", "17:3076578:9167836")
trisChr <- c("16","17")
attributes = listAttributes(ensembl);
		attrib <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),
			as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='start_position'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='end_position'),][1,])),
			as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,])));
		attrib <-attributes[attrib,1]; # "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "gene_biotype"  
		genesDf_allTs65Dn = getBM(attributes=attrib,filters="chromosomal_region",values=reg_Ts65dn, mart=ensembl); 
		genesDf_isof.Ts65Dn <- genesDf_allTs65Dn[which(genesDf_allTs65Dn$chromosome_name %in% trisChr), ]; #58676 nb of genes including isoforms
		genesDf_uniq.Ts65Dn <- genesDf_isof.Ts65Dn[!duplicated(genesDf_isof.Ts65Dn$external_gene_name), ]; #57152 #nb of uniquely genes, excluding isoforms
		genesDf.Ts65Dn <- left_join(genesDf_uniq.Ts65Dn,biotypeCat,by="gene_biotype")
		aRes<- as.data.frame(genesDf.Ts65Dn[,c(2,3,7)] %>%group_by(chromosome_name,biotype_category) %>% count())
		bRes<- as.data.frame(genesDf.Ts65Dn[,c(2,3,7)] %>%group_by(biotype_category) %>% count())
		cRes<- as.data.frame(genesDf.Ts65Dn[,c(2,3,6)] %>%group_by(chromosome_name,gene_biotype) %>% count())
		Ts65Dn_homRes<- list(genesDf.Ts65Dn=genesDf.Ts65Dn,statsByChr=aRes, statsBybiotypeCat=bRes, statsByGeneBiotype=cRes)
	write.xlsx(Ts65Dn_homRes, file = paste0(wd,f_results, f_trisomic,"Ts65Dn_chrBiotype_Stats_genes.xlsx"));

#are these genes considered essential genes? as published by dickinson2016_essentialGenes

input_essentialG<- "/results/IMPC/dickinson2016_essentialGenes/"

lethalG <-read.xlsx(paste0(wd,input_essentialG, "NIHMS809908-supplement-supp_table1.xlsx"),2)
subviableG <-read.xlsx(paste0(wd,input_essentialG, "NIHMS809908-supplement-supp_table1.xlsx"),3)
essentialG <- rbind(lethalG,subviableG)

dfmouse<-resultsHomologs_Mouse$genesDf_specie
essentialMouse<- dfmouse[which(dfmouse$associated_gene_name %in% essentialG$gene_symbol),]
essentialMouse$associated_gene_name

##################################################
##################################################
# Figures:
##################################################
##################################################

#mod from LeilaHC  
### https://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot/26749522#26749522

donuts_plot <- function(
                        panel = runif(3), # counts
                        pctr = c(.5,.2,.9), # percentage in count
                        legend.label='',
                        cols = c('chartreuse', 'chocolate','deepskyblue'), # colors
                        outradius = 1, # outter radius
                        radius = .7,   # 1-width of the donus 
                        add = F,
                        innerradius = .5, # innerradius, if innerradius==innerradius then no suggest line
                        legend = F,
                        pilabels=F,
                        legend_offset=.25, # non-negative number, legend right position control
                        borderlit=c(F,F,T,T)
                        ){
    par(new=add,mar=c(14,14,14,16))
    if(sum(legend.label=='')>=1) legend.label=paste("Series",1:length(pctr))
    if(pilabels){
        pie(panel, col=cols,border = borderlit[1],labels = legend.label,radius = outradius)
    }
    panel = panel/sum(panel)

    pctr2= panel*(1 - pctr)
    pctr3 = c(pctr,pctr)
    pctr_indx=2*(1:length(pctr))
    pctr3[pctr_indx]=pctr2
    pctr3[-pctr_indx]=panel*pctr
    cols_fill = c(cols,cols)
    cols_fill[pctr_indx]='white'
    cols_fill[-pctr_indx]=cols
    par(new=TRUE,mar=c(14,14,14,16))
    pie(pctr3, col=cols_fill,border = borderlit[2],labels = '',radius = outradius)
    par(new=TRUE,mar=c(14,14,14,16))
    pie(panel, col='white',border = borderlit[3],labels = '',radius = radius)
    par(new=TRUE,mar=c(14,14,14,16))
    pie(1, col='white',border = borderlit[4],labels = '',radius = innerradius)
    if(legend){
        par(mar=c(6,13,32,6), xpd=TRUE)
        legend("topright",inset=c(-legend_offset,0),pt.cex=1.2,legend=legend.label, cex=0.8,pch=rep(15,'.',length(pctr)), 
               col=cols,bty='n')
    }
    par(new=FALSE)
}
## col- > subcor(change hue/alpha)
subcolors <- function(.dta,main,mainCol){
    tmp_dta = cbind(.dta,1,'col')
    tmp1 = unique(.dta[[main]])
    for (i in 1:length(tmp1)){
        tmp_dta$"col"[.dta[[main]] == tmp1[i]] = mainCol[i]
    }
    u <- unlist(by(tmp_dta$"1",tmp_dta[[main]],cumsum))
    n <- dim(.dta)[1]
    subcol=rep(rgb(0,0,0),n);
    for(i in 1:n){
        t1 = col2rgb(tmp_dta$col[i])/256
        subcol[i]=rgb(t1[1],t1[2],t1[3],1/(1+u[i]))
    }
    return(subcol);
}


#donutplot for %categories
catDf<- bioCat[,c(1,3,7,11,15,19,23,27,31,35)]
catDf[is.na(catDf)] <- 0
catDf$biotype_category <- factor(catDf$biotype_category, levels = unique(catDf$biotype_category))

#as there is some superposition, changing the order or categories to plot
catDf1<-catDf[c(1:4),c(1,5)]
catDf1<-catDf1[c(1,2,4,3),]
catDf2<-catDf[,c(1,6)]
catDf2<-catDf2[c(1,2,4,3,5),]

colDonutCat<- c("#00C5CD","#FFFF99","#993399","#CC0066","#FF9933")
colDonutCat1<- c("#00C5CD","#FFFF99","#CC0066","#993399")
colDonutCat2<- c("#00C5CD","#FFFF99","#CC0066","#993399","#FF9933")


pdf(paste0(wd, f_results,"_categories_donut.pdf"),width = 8.27, height=11.69)  #8.27 × 11.69)
#1 "% Human"            
donuts_plot(catDf[,2],rep(1,5),catDf[,1],
        cols=colDonutCat,pilabels=F,legend=T,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  
plot.new()
#2"% Human" 
donuts_plot(catDf[,2],rep(1,5),paste0(catDf[,2]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  
plot.new()
# 3"% Mouse"  
donuts_plot(catDf[,3],rep(1,5),paste0(catDf[,3]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  


plot.new()
#"4% Rat" 
donuts_plot(catDf[,4],rep(1,5),paste0(catDf[,4]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  


plot.new()
#"5% Gorilla"
donuts_plot(catDf[,7],rep(1,5),paste0(catDf[,7]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  


plot.new()
#"6% Orangutan" 
donuts_plot(catDf[,8],rep(1,5),paste0(catDf[,8]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  


plot.new()
#"7% Zebrafish"
donuts_plot(catDf[,9],rep(1,5),paste0(catDf[,9]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  



plot.new()
#8"% C.elegands"
donuts_plot(catDf[,10],rep(1,5),paste0(catDf[,10]," %"),
        cols=colDonutCat,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  


plot.new()
#9"% Drosophila" 
donuts_plot(catDf1[,2],rep(1,4),paste0(catDf1[,2]," %"),
        cols=colDonutCat1,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,4) )  


plot.new()
#10"% Chimpanzee" 
donuts_plot(catDf2[,2],rep(1,5),paste0(catDf2[,2]," %"),
        cols=colDonutCat2,pilabels=T,legend=F,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,5) )  


dev.off()
library(openxlsx)
orderedSpecies <- c("Human","Gorilla","Chimpanzee","Rat","Mouse","Orangutan","Zebrafish","Drosophila","C.elegands")
labels <- c("Human","Gorilla","Chimpanzee","Rat","Mouse","Orangutan","Zebrafish","Drosophila","c.elegands")
size <-read.xlsx("/Users/marmmoreno/Documents/YH_Lab/collabs_outside/DS_models_book/results/tables/TablesBook.xlsx",2)
size <- size[order(match(size$Specie,orderedSpecies )),]
size$labels <-  c("Human","Gorilla","Chimpanzee","Rat","Mouse","Orangutan","Zebrafish","FruitFly","Worm")
size$sizeCol <- "#ffaaf0"	

colnames(sizes)[2] <- "GenomeSize"
pathsTopAll$q.valMod <- 1 -as.numeric(pathsTopAll$Specie)
pathsTopAll$q.valMod <- factor(pathsTopAll$q.valMod, levels = unique(pathsTopAll$q.valMod))
pathsTopAll$q.val <- factor(pathsTopAll$q.val, levels = unique(pathsTopAll$q.val))
pathsTopAll$GeneSetName <- factor(pathsTopAll$GeneSetName, levels = pathsTopAll$GeneSetName)

#library(ggplot2)
regCol <- c("#C51B7D","#762A83")  #colour #1:down #2:up
legend_title <- "Genome size (Bps)"

plotBar <- ggplot(size, aes(Specie,GenomeSize, fill = "#CC0066")) +
		geom_bar(position="dodge", stat="identity") + 
		scale_fill_manual(" " , labels="" ,values="#CC0066") + 
        coord_flip() + labs(y =paste0(size$labels)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plotBar
	
labels = scales::scientific,
breaks <- c("1000000","10000000","100000000","1000000000","10000000000")

ggsave(paste0(wd,f_results, f_gage,f_name,f_dfa, name,"_top",numberTopPaths,"_",GeneSet,".pdf", sep=""),  plotTopPaths, width = width, height = height,scale=scale, units = "in")



##################################################################################
#C. Human chr21 genes: orthologs genes in the different species 
##################################################################################


chr_homologs.function <- function(specieReferenceBiomaRt,specieReferenceName,chrNameSpecieReference,speciesQueryBiomaRt,speciesQueryNames,biotypeCat,dateHom){
	##################################################################################
	# Function to annotate all the genes of the chr of the reference specie with
	# its homologs in the species of interest defined as query
	# Its produce the stats. Input files are:
	#
	##################################################################################	

	#2. Quering ENSEMBL 
	######################
	#ensembl <- useEnsembl(biomart = "ensembl")
	#listDatasets(mart = ensembl)
	#for the reference specie:
	
	ensembl = useMart("ensembl", dataset = specieReferenceBiomaRt); 
	attributes = listAttributes(ensembl);
	attributesHom = attributes[which(attributes[,3]=="homologs"),];	
	rownames(attributesHom) <- 1:nrow(attributesHom);
	mouseInfo<-  attributesHom[1342:1357,c(1,2)]

	attribRef <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='start_position'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='end_position'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,]))  );
	attribRef <-attributes[attribRef,1];  
	genesDf_ref = getBM(attributes=attribRef,filters="chromosome_name",values=chrNameSpecieReference, mart=ensembl); 
	values<-unique(genesDf_ref[,1]); 
	assign(paste0("genesDf_ref", "_", specieReferenceBiomaRt),genesDf_ref,.GlobalEnv);

	#paralogs human
	###########################3
	f_paralogs <- "paralogs/"
	dir.create(file.path(getwd (), f_results,f_homologs, f_paralogs), showWarnings = F);

	attribPar <-c(as.numeric(rownames(attributesHom[which(attributesHom[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributesHom[which(attributesHom[,1]=='external_gene_name'),][1,])),
		as.numeric(rownames(attributesHom[grep("hsapiens",attributesHom[,1]),]))  );
	attribPar <-attributesHom[attribPar,1];  
	
	genesDf_paralogs = getBM(attributes=attribPar,filters="ensembl_gene_id",values=values, mart=ensembl); 
	write.xlsx(genesDf_paralogs, file = paste0(wd,f_results,f_homologs,f_paralogs, specieReferenceName, "_specieRef_paralogs_defined_by_ensembl_",dateHom, ".xlsx"));
	assign(paste0("genesDf_paralogs_",specieReferenceName),genesDf_paralogs,.GlobalEnv);

	perchrHomologsFList <- list();
	for (i in 1:length(speciesQueryBiomaRt)){	
		f_specie <- paste0(gsub(" ","_",speciesQueryNames[i]),"/");
		name <- speciesQueryNames[i];
		specieAnalysis <- speciesQueryBiomaRt[i];

		#1. Creating the folders
		######################
		dir.create(file.path(getwd (), f_results,f_homologs, f_specie), showWarnings = F);
		dir.create(file.path(getwd (), f_results,f_homologs, f_stats), showWarnings = F);
		dir.create(file.path(getwd (), f_results,f_homologs, f_specie, f_figures), showWarnings = F);
		dir.create(file.path(getwd (), f_results,f_homologs, f_specie, f_tables), showWarnings = F);
		#dir.create(file.path(getwd (), f_results,f_homologs, f_specie, f_Rdata), showWarnings = F);
		
		#2. Quering ENSEMBL 
		######################
		#ensembl <- useEnsembl(biomart = "ensembl")
		#listDatasets(mart = ensembl)
		#for the reference specie:

		#per specie orthologs
		##################################################
		#meaning of each column: ex

		attribDSpecie <-c(as.numeric(rownames(attributesHom[which(attributesHom[,1]=='ensembl_gene_id'),][1,])),
			as.numeric(rownames(attributesHom[grep(speciesQueryBiomaRt[i],attributesHom[,1]),]))  );
		attribSpecie <-attributesHom[attribDSpecie,1];  

		genesDf_specie = getBM(attributes=attribSpecie,filters="ensembl_gene_id",values=values, mart=ensembl);
		genesDf_specie <- left_join(genesDf_specie,genesDf_ref[,c(1,2,6)], by="ensembl_gene_id");
		genesDf_specie <- left_join(genesDf_specie,biotypeCat,by="gene_biotype");
		genesDf_specie<- genesDf_specie[,c(grep("ensembl_gene_id" ,colnames(genesDf_specie)),grep("external_gene_name" ,colnames(genesDf_specie)),grep("biotype_category" ,colnames(genesDf_specie)),grep("gene_biotype" ,colnames(genesDf_specie)),grep("homolo" ,colnames(genesDf_specie))  )];
		colnames(genesDf_specie) <- gsub(".*_homolog_","",colnames(genesDf_specie));
		genesDf_specie <-genesDf_specie[-which(genesDf_specie$ensembl_gene==""),];# not homolog found #all homologs		
		remoteHomologs <- genesDf_specie[which(genesDf_specie[,20]==0 & genesDf_specie[,14] >15),]; #remote, confidence =0, seq similarity > 15% i the query animal
		remoteHomologs$HomologyCategory <- "remote homologs";
		homologsHighConf<- genesDf_specie[which(genesDf_specie[,20]==1),]; # high confidence: conf score=1
		homologsHighConf$HomologyCategory <- "High confidence homologs";
		predGenes <- genesDf_specie[is.na(genesDf_specie[,20]),]; #creating this category because, the conf level of these genes is 0, but the % seq shared is high, 
		predGenes$HomologyCategory <- "predicted genes";
		genesDf_specieF<- rbind(homologsHighConf,remoteHomologs,predGenes);

		#and it may be cause most of these genes are predicted in the specie of reference and we dont have enough proof of them yet that is why the conf score is 0 in the end
		#colnames(genesDf_specieF) <- gsub(".*_homolog_","",colnames(genesDf_specieF));
		
		#per chr homologs with names in the query specie
		##################################################
		df<- genesDf_specieF;
		orthology_typeOrder <- c("ortholog_one2one","ortholog_one2many","ortholog_many2many");
		homology_typeOrder <- c("High confidence homologs","predicted genes","remote homologs");

		HomologyCategoryBioM <-as.data.frame(unique(df[,c(2,3,21)]) %>% group_by(HomologyCategory,biotype_category) %>% count());
		colnames(HomologyCategoryBioM)[3] <- paste0("Nb ", specieReferenceName, " homologs found in ",name," per biotype");
		HomologyCategoryBioM2 <-as.data.frame(df %>% group_by(HomologyCategory,biotype_category) %>% count());
		colnames(HomologyCategoryBioM2)[c(3)] <- paste0("Nb ",name," homologs to ",specieReferenceName," per biotype");
		HomologyCategoryBioM2 <- HomologyCategoryBioM2[order(match(HomologyCategoryBioM2[,2],unique(biotypeCat[,2]))),];
		HomologyCategoryBioM2 <- HomologyCategoryBioM2[order(match(HomologyCategoryBioM2[,1],homology_typeOrder)),];
		HomologyCategoryBioM <- left_join(HomologyCategoryBioM,HomologyCategoryBioM2,by=c("HomologyCategory","biotype_category"));
		perchrHomologsF <- as.data.frame(unique(df[,c(2,3,6,13,21,8)]) %>% group_by(HomologyCategory,biotype_category,orthology_type,chromosome) %>% count());
		#colnames(perchrHomologs)[c(4,5)] <- c(paste0(name," chromosome "), paste0("Nb ",name," homologs to ",specieReferenceName," per chr"));

		perchrHomologsNamesRef <-  as.data.frame(unique(df[,c(2,3,13,21,8)]) %>% group_by(HomologyCategory,biotype_category,orthology_type,chromosome) %>% count());
		colnames(perchrHomologsNamesRef)[5] <- paste0("Nb of ", specieReferenceName ," homologs found in ",name, " per chr");
		
		perchrHomologsNames <- as.data.frame(unique(df[,c(2,3,6,13,21,8)]) %>% group_by(HomologyCategory,biotype_category,orthology_type,chromosome));
		perchrHomologsNames$sumup <- paste0(perchrHomologsNames[,2],":",perchrHomologsNames[,4],":",perchrHomologsNames[,5],":",perchrHomologsNames[,6]);
		#and it may be cause most of these genes are predicted in the specie of reference and we dont have enough proof of them yet that is why the conf score is 0 in the end
		sumup <- unique(aggregate(associated_gene_name ~sumup , data = perchrHomologsNames, paste, collapse = ", "));
		sumup<- separate(sumup, sumup, c("biotype_category" ,"orthology_type", "HomologyCategory", "chromosome"  ),sep=":");
		#and it may be cause most of these genes are predicted in the specie of reference and we dont have enough proof of them yet that is why the conf score is 0 in the end
		perchrHomologsF <- left_join(perchrHomologsF,sumup, by=colnames(sumup[,c(4)]))
		#and it may be cause most of these genes are predicted in the specie of reference and we dont have enough proof of them yet that is why the conf score is 0 in the end
		sumup2 <- unique(aggregate(external_gene_name ~sumup , data = perchrHomologsNames, paste, collapse = ", "));
		sumup2<- separate(sumup2, sumup, c("biotype_category" ,"orthology_type", "HomologyCategory", "chromosome"  ),sep=":");
		perchrHomologsF <- left_join(perchrHomologsF,sumup2, by=colnames(sumup[,c(4)]));
		perchrHomologsF$totalHomologs<- length(unique(df[,c(2)])) #total nb of unique reference specie genes 
		colnames(perchrHomologsF)[c(5,8)] <- c(paste0("Nb of ", name," homologs to ",specieReferenceName ," per chr"), paste0("Total nb of ", specieReferenceName," homologs found in ",name));	
		perchrHomologsF<- left_join(perchrHomologsF,HomologyCategoryBioM,by=c("HomologyCategory","biotype_category"));
		perchrHomologsF$perc1 <-round((perchrHomologsF[,5]/perchrHomologsF[,8])*100,digits=2);
		perchrHomologsF$perc2 <-round((perchrHomologsF[,5]/perchrHomologsF[,10])*100,digits=2);
		colnames(perchrHomologsF)[c(11,12)] <- c(paste0("% ", name," homologs to per total ",specieReferenceName, " homologs"), paste0("%  ", name," homologs to ",specieReferenceName," per biotype "));


		perchrHomologsF <- perchrHomologsF[order(perchrHomologsF[,5],decreasing = TRUE ),];
		perchrHomologsF <- perchrHomologsF[order(match(perchrHomologsF[,3],orthology_typeOrder )),];
		perchrHomologsF <- perchrHomologsF[order(match(perchrHomologsF[,2],unique(biotypeCat[,2]) )),];
		perchrHomologsF <- perchrHomologsF[order(match(perchrHomologsF[,1],homology_typeOrder)),];
		perchrHomologsF <-perchrHomologsF[,c(1:5,8:12,6,7)]
		perchrHomologsF <- left_join(perchrHomologsF,perchrHomologsNamesRef,by=colnames(perchrHomologsNamesRef)[c(1:4)])
		perchrHomologsF$perc1 <-round((perchrHomologsF[,13]/perchrHomologsF[,6])*100,digits=2);
		perchrHomologsF$perc2 <-round((perchrHomologsF[,13]/perchrHomologsF[,7])*100,digits=2);
		colnames(perchrHomologsF)[c(14,15)] <- c(paste0("% ",specieReferenceName," homologs to per total ",name," homologs"),paste0("% ",specieReferenceName," homologs per total ",name," homologs per chr"))
		perchrHomologsF<- perchrHomologsF[,c(1:4,13,6,7,14:15,5,8,9:10,11,12)]
		resultsHomologs <- list(genesDf_specie=genesDf_specie,homologsHighConf=homologsHighConf,remoteHomologs=remoteHomologs,predGenes=predGenes,perchrHomologsF=perchrHomologsF);
		write.xlsx(resultsHomologs, file = paste0(wd,f_results,f_homologs, f_specie, f_tables, specieReferenceName, "_specieRef_homologs_with_query_specie_",speciesQueryNames[i],"_",dateHom,".xlsx"));
		perchrHomologsFList[[i]] <- perchrHomologsF
		names(perchrHomologsFList)[i] <- paste0(name, "_perchrHomologsF")
		#saving the variables for each animal
		##################################################

		assign(paste0("genesDf_specie","_",name),genesDf_specie,.GlobalEnv);
		assign(paste0("homologsHighConf","_",name),homologsHighConf,.GlobalEnv);
		assign(paste0("remoteHomologs","_",name),remoteHomologs,.GlobalEnv);
		assign(paste0("predGenes","_",name),predGenes,.GlobalEnv);
		assign(paste0("resultsHomologs", "_", name),resultsHomologs,.GlobalEnv);
		assign(paste0("perchrHomologsF", "_", name),perchrHomologsF,.GlobalEnv);



		# Calculating stats:
		##################################################

		## A. GenomeWide stats
		######################

		if (i==1){
			
			#per homology category
			##############################
			#length(unique(df[,c(2)])) #total nb of unique reference specie genes
			#length(unique(df[,6]))  #total nb of unique query specie homologs genes
			#1
			HomologyCategory <-as.data.frame(unique(df[,c(2,21)]) %>% group_by(HomologyCategory) %>% count()); #nb of unique human genes with at least one homolog. 
			HomologyCategory$total <- length(unique(df[,c(2)])) #total nb of unique reference specie genes
			
			HomologyCategory$perc <- round((HomologyCategory[,2]/HomologyCategory[,3])*100,digits=2);			#calculating the %
			colnames(HomologyCategory)[c(2:4)] <- c(paste0("Nb ",specieReferenceName, " homologs found in ",name), paste0("Total nb ", specieReferenceName, " homologs found in ",name), paste0("% ",specieReferenceName," homologs found in ",name));
			#2
			HomologyCategory2 <-as.data.frame(df %>% group_by(HomologyCategory) %>% count()); #nb of unique human genes with at least one homolog. 
			HomologyCategory2$total <- length(unique(df[,c(6)])) #total nb of unique query specie homologs genes
			HomologyCategory2$perc <- round((HomologyCategory2[,2]/HomologyCategory2[,3])*100,digits=2);			#calculating the %
			colnames(HomologyCategory2)[c(2:4)] <- c(paste0("Nb ",name," homologs to ",specieReferenceName ), paste0("Total nb ", name , " homologs found in ",specieReferenceName), paste0("% ",name," homologs to ",specieReferenceName));
			#3

			HomologyCategory <- left_join(HomologyCategory,HomologyCategory2,by="HomologyCategory");
			#HomologyCategory[nrow(HomologyCategory)+1,] <- c("total", round(colSums(HomologyCategory[,c(2:5)]), digits=0));

			#per bio category
			##############################
			
			#1
			HomologyCategoryBio <-as.data.frame(unique(df[,c(2,3,21)]) %>% group_by(HomologyCategory,biotype_category) %>% count());
			HomologyCategoryBio$total <- length(unique(df[,c(2)])) #total nb of unique reference specie genes
			HomologyCategoryBio$perc <- round((HomologyCategoryBio[,3]/HomologyCategoryBio[,4])*100,digits=2);			#calculating the %
			colnames(HomologyCategoryBio)[c(3:5)] <- c(paste0("Nb ", specieReferenceName, " homologs found in ",name," per biotype"), paste0("Total nb ", specieReferenceName, " homologs found in ",name), paste0("% ",specieReferenceName," homologs found in ",name," per biotype"));
			#2
			HomologyCategoryBio2 <-as.data.frame(df %>% group_by(HomologyCategory,biotype_category) %>% count());
			HomologyCategoryBio2$total <- length(unique(df[,c(6)])) #total nb of unique query specie homologs genes
			HomologyCategoryBio2$perc <- round((HomologyCategoryBio2[,3]/HomologyCategoryBio2[,4])*100,digits=2);			#calculating the %
			colnames(HomologyCategoryBio2)[c(3:5)] <- c(paste0("Nb ",name," homologs to ",specieReferenceName," per biotype" ), paste0("Total nb ", name , " homologs found in ",specieReferenceName), paste0("% ",name," homologs to ",specieReferenceName," per biotype"));
			#3
			HomologyCategoryBio <- left_join(HomologyCategoryBio,HomologyCategoryBio2,by=c("HomologyCategory","biotype_category"));
			HomologyCategoryBio <- HomologyCategoryBio[order(match(HomologyCategoryBio[,2],unique(biotypeCat[,2]))),];
			HomologyCategoryBio <- HomologyCategoryBio[order(match(HomologyCategoryBio[,1],homology_typeOrder)),];
			#HomologyCategoryBio[nrow(HomologyCategoryBio)+1,] <- c("total","", round(colSums(HomologyCategoryBio[,c(3:6)]), digits=0));

			#stats per  type of homology 
			##############################
			#1
			HomologyCategoryBio.queryHomologs <-as.data.frame(unique(df[,c(2,3,13,21)]) %>% group_by(HomologyCategory,biotype_category,orthology_type) %>% count());
			HomologyCategoryBio.queryHomologs <- left_join(HomologyCategoryBio.queryHomologs,HomologyCategory[,c(1:2,4)],by="HomologyCategory");
			HomologyCategoryBio.queryHomologs <- left_join(HomologyCategoryBio.queryHomologs,HomologyCategoryBio[,c(1:3,5)],by=c("HomologyCategory","biotype_category"));
			HomologyCategoryBio.queryHomologs$perc1 <- round((HomologyCategoryBio.queryHomologs[,4]/as.numeric(HomologyCategoryBio.queryHomologs[,5]))*100,digits=2);			#calculating the % respect to nb genes per homology category
			HomologyCategoryBio.queryHomologs$perc2 <- round((HomologyCategoryBio.queryHomologs[,4]/as.numeric(HomologyCategoryBio.queryHomologs[,7]))*100,digits=2);			#calculating the %respect to nb genes per homology category&  biotype category
			colnames(HomologyCategoryBio.queryHomologs)[c(4,9,10)] <- c(paste0("Nb ",specieReferenceName, " homologs found in ",name," per orthology type"), paste0("% ",specieReferenceName," homologs found in ",name," per Homology Category & orthology type"),paste0("% ",specieReferenceName," homologs found in ",name," per biotype Category & orthology type"));
			
			#2
			HomologyCategoryBio.queryHomologs2 <-as.data.frame(df %>% group_by(HomologyCategory,biotype_category,orthology_type) %>% count());
			HomologyCategoryBio.queryHomologs2 <- left_join(HomologyCategoryBio.queryHomologs2,HomologyCategory[,c(1:2,4)],by="HomologyCategory");
			HomologyCategoryBio.queryHomologs2 <- left_join(HomologyCategoryBio.queryHomologs2,HomologyCategoryBio[,c(1:3,5)],by=c("HomologyCategory","biotype_category"));
			HomologyCategoryBio.queryHomologs2$perc1 <- round((HomologyCategoryBio.queryHomologs2[,4]/as.numeric(HomologyCategoryBio.queryHomologs2[,6]))*100,digits=2);			#calculating the % respect to nb genes per homology category
			HomologyCategoryBio.queryHomologs2$perc2 <- round((HomologyCategoryBio.queryHomologs2[,4]/as.numeric(HomologyCategoryBio.queryHomologs2[,8]))*100,digits=2);			#calculating the %respect to nb genes per homology category&  biotype category
			#3
			#reorder
			HomologyCategoryBio.queryHomologs <- HomologyCategoryBio.queryHomologs[order(HomologyCategoryBio.queryHomologs[,9],decreasing = TRUE ),];
			HomologyCategoryBio.queryHomologs <- HomologyCategoryBio.queryHomologs[order(match(HomologyCategoryBio.queryHomologs[,2],unique(biotypeCat[,2]) )),];
			HomologyCategoryBio.queryHomologs <- HomologyCategoryBio.queryHomologs[order(match(HomologyCategoryBio.queryHomologs[,1],homology_typeOrder )),];
			#HomologyCategoryBio.queryHomologs[nrow(HomologyCategoryBio.queryHomologs)+1,] <- c("total","","", round(colSums(as.numeric(HomologyCategoryBio.queryHomologs[,c(4:13)])), digits=0))
			#stats per type of homology chr
			##############################
			perchrHomologs <- as.data.frame(unique(df[,c(2,3,6,13,21,8)]) %>% group_by(HomologyCategory,biotype_category,orthology_type,chromosome) %>% count());
			colnames(perchrHomologs)[c(4,5)] <- c("chromosome", paste0("Nb ",name," homologs to ",specieReferenceName," per chr"));
			orthology_typeOrder <- c("ortholog_one2one","ortholog_one2many","ortholog_many2many")
			perchrHomologs <- perchrHomologs[order(perchrHomologs[,5],decreasing = TRUE ),];
			perchrHomologs <- perchrHomologs[order(match(perchrHomologs[,3],orthology_typeOrder )),];
			perchrHomologs <- perchrHomologs[order(match(perchrHomologs[,2],unique(biotypeCat[,2]) )),];
			perchrHomologs <- perchrHomologs[order(match(perchrHomologs[,1],homology_typeOrder )),];

			perchrHomologsNoOrtCat <- as.data.frame(unique(df[,c(2,3,6,21,8)]) %>% group_by(HomologyCategory,biotype_category,chromosome) %>% count());
			colnames(perchrHomologsNoOrtCat)[c(3,4)] <- c("chromosome", paste0("Nb ",name," homologs to ",specieReferenceName," per chr"));
			perchrHomologsNoOrtCat$total <-length(unique(df[,c(2)])); #total nb of unique reference specie genes
			perchrHomologsNoOrtCat <- left_join(perchrHomologsNoOrtCat,HomologyCategoryBio[,c(1:2,5,3)], by=c("HomologyCategory","biotype_category"));
			colnames(perchrHomologsNoOrtCat)[5] <-  paste0(name, " total number of homologous genes found to ",specieReferenceName);
			perchrHomologsNoOrtCat<- perchrHomologsNoOrtCat[,c(1:3,5,4,6,7)];
			perchrHomologsNoOrtCat <- perchrHomologsNoOrtCat[order(perchrHomologsNoOrtCat[,5],decreasing = TRUE ),];
			perchrHomologsNoOrtCat <- perchrHomologsNoOrtCat[order(match(perchrHomologsNoOrtCat[,2],unique(biotypeCat[,2]) )),];				
			perchrHomologsNoOrtCat <- perchrHomologsNoOrtCat[order(match(perchrHomologsNoOrtCat[,1],homology_typeOrder )),];				
					
	
		} else {

			#per homology category
			##############################
			#length(unique(df[,c(2)])) #total nb of unique reference specie genes
			#length(unique(df[,6]))  #total nb of unique query specie homologs genes
			#1
			HomologyCategory.tmp <-as.data.frame(unique(df[,c(2,21)]) %>% group_by(HomologyCategory) %>% count()); #nb of unique human genes with at least one homolog. 
			HomologyCategory.tmp$total <- length(unique(df[,c(2)])) #total nb of unique reference specie genes
			HomologyCategory.tmp$perc <- round((HomologyCategory.tmp[,2]/HomologyCategory.tmp[,3])*100,digits=2);			#calculating the %
			colnames(HomologyCategory.tmp)[c(2:4)] <- c(paste0("Nb ",specieReferenceName, " homologs found in ",name), paste0("Total nb ", specieReferenceName, " homologs found in ",name), paste0("% ",specieReferenceName," homologs found in ",name));
			#2
			HomologyCategory2.tmp <-as.data.frame(df %>% group_by(HomologyCategory) %>% count()); #nb of unique human genes with at least one homolog. 
			HomologyCategory2.tmp$total <- length(unique(df[,c(6)])) #total nb of unique query specie homologs genes
			HomologyCategory2.tmp$perc <- round((HomologyCategory2.tmp[,2]/HomologyCategory2.tmp[,3])*100,digits=2);			#calculating the %
			colnames(HomologyCategory2.tmp)[c(2:4)] <- c(paste0("Nb ",name," homologs to ",specieReferenceName ),paste0("Total nb ", name , " homologs found in ",specieReferenceName),  paste0("% ",name," homologs to ",specieReferenceName));
			#3
			HomologyCategory.tmp <- left_join(HomologyCategory.tmp,HomologyCategory2.tmp,by="HomologyCategory");
			HomologyCategory <- unique(full_join(HomologyCategory,HomologyCategory.tmp, by="HomologyCategory"));
			#per bio category
			##############################
			
			#1
			HomologyCategoryBio.tmp <-as.data.frame(unique(df[,c(2,3,21)]) %>% group_by(HomologyCategory,biotype_category) %>% count());
			HomologyCategoryBio.tmp$total <- length(unique(df[,c(2)])) #total nb of unique reference specie genes
			HomologyCategoryBio.tmp$perc <- round((HomologyCategoryBio.tmp[,3]/HomologyCategoryBio.tmp[,4])*100,digits=2);			#calculating the %
			colnames(HomologyCategoryBio.tmp)[c(3:5)] <- c(paste0("Nb ", specieReferenceName, " homologs found in ",name," per biotype"), paste0("Total nb ", specieReferenceName, " homologs found in ",name), paste0("% ",specieReferenceName," homologs found in ",name," per biotype"));
			#2
			HomologyCategoryBio2.tmp <-as.data.frame(df %>% group_by(HomologyCategory,biotype_category) %>% count());
			HomologyCategoryBio2.tmp$total <- length(unique(df[,c(6)])) #total nb of unique query specie homologs genes
			HomologyCategoryBio2.tmp$perc <- round((HomologyCategoryBio2.tmp[,3]/HomologyCategoryBio2.tmp[,4])*100,digits=2);			#calculating the %
			colnames(HomologyCategoryBio2.tmp)[c(3:5)] <- c(paste0("Nb ",name," homologs to ",specieReferenceName," per biotype" ), paste0("Total nb ", name , " homologs found in ",specieReferenceName), paste0("% ",name," homologs to ",specieReferenceName," per biotype"));
			#3
			HomologyCategoryBio.tmp <- left_join(HomologyCategoryBio.tmp,HomologyCategoryBio2.tmp,by=c("HomologyCategory","biotype_category"));
			HomologyCategoryBio <- unique(full_join(HomologyCategoryBio,HomologyCategoryBio.tmp, by=c("HomologyCategory", "biotype_category")));			
			#HomologyCategoryBio.tmp[nrow(HomologyCategoryBio.tmp)+1,] <- c("total","", round(colSums(HomologyCategoryBio.tmp[,c(3:6)]), digits=0))

			#stats per  type of homology 
			##############################
			#1
			HomologyCategoryBio.queryHomologs.tmp <-as.data.frame(unique(df[,c(2,3,13,21)]) %>% group_by(HomologyCategory,biotype_category,orthology_type) %>% count());
			HomologyCategoryBio.queryHomologs.tmp <- left_join(HomologyCategoryBio.queryHomologs.tmp,HomologyCategory.tmp[,c(1:2,4)],by="HomologyCategory");
			HomologyCategoryBio.queryHomologs.tmp <- left_join(HomologyCategoryBio.queryHomologs.tmp,HomologyCategoryBio.tmp[,c(1:3,5)],by=c("HomologyCategory","biotype_category"));
			HomologyCategoryBio.queryHomologs.tmp$perc1 <- round((HomologyCategoryBio.queryHomologs.tmp[,4]/as.numeric(HomologyCategoryBio.queryHomologs.tmp[,5]))*100,digits=2);			#calculating the % respect to nb genes per homology category
			HomologyCategoryBio.queryHomologs.tmp$perc2 <- round((HomologyCategoryBio.queryHomologs.tmp[,4]/as.numeric(HomologyCategoryBio.queryHomologs.tmp[,7]))*100,digits=2);			#calculating the %respect to nb genes per homology category&  biotype category
			colnames(HomologyCategoryBio.queryHomologs.tmp)[c(4,9,10)] <- c(paste0("Nb ",specieReferenceName, " homologs found in ",name," per orthology type"), paste0("% ",specieReferenceName," homologs found in ",name," per Homology Category & orthology type"),paste0("% ",specieReferenceName," homologs found in ",name," per biotype Category & orthology type"));
			
			#2
			HomologyCategoryBio.queryHomologs2.tmp <-as.data.frame(df %>% group_by(HomologyCategory,biotype_category,orthology_type) %>% count());
			HomologyCategoryBio.queryHomologs2.tmp <- left_join(HomologyCategoryBio.queryHomologs2.tmp,HomologyCategory.tmp[,c(1:2,4)],by="HomologyCategory");
			HomologyCategoryBio.queryHomologs2.tmp <- left_join(HomologyCategoryBio.queryHomologs2.tmp,HomologyCategoryBio.tmp[,c(1:3,5)],by=c("HomologyCategory","biotype_category"));
			HomologyCategoryBio.queryHomologs2.tmp$perc1 <- round((HomologyCategoryBio.queryHomologs2.tmp[,4]/as.numeric(HomologyCategoryBio.queryHomologs2.tmp[,6]))*100,digits=2);			#calculating the % respect to nb genes per homology category
			HomologyCategoryBio.queryHomologs2.tmp$perc2 <- round((HomologyCategoryBio.queryHomologs2.tmp[,4]/as.numeric(HomologyCategoryBio.queryHomologs2.tmp[,8]))*100,digits=2);			#calculating the %respect to nb genes per homology category&  biotype category
			colnames(HomologyCategoryBio.queryHomologs2.tmp)[c(4,9,10)] <- c(paste0("Nb ",name," homologs to ",specieReferenceName," per orthology type"), paste0("% ",name," homologs to ",specieReferenceName," per Homology Category & orthology type"),paste0("% ",name," homologs to ",specieReferenceName," per biotype Category & orthology type"));
			#3
			#reorder
			HomologyCategoryBio.queryHomologs.tmp <- HomologyCategoryBio.queryHomologs.tmp[order(HomologyCategoryBio.queryHomologs.tmp[,9],decreasing = TRUE ),];
			HomologyCategoryBio.queryHomologs.tmp <- HomologyCategoryBio.queryHomologs.tmp[order(match(HomologyCategoryBio.queryHomologs.tmp[,2],unique(biotypeCat[,2]) )),];
			HomologyCategoryBio.queryHomologs.tmp <- HomologyCategoryBio.queryHomologs.tmp[order(match(HomologyCategoryBio.queryHomologs.tmp[,1],homology_typeOrder )),];
			HomologyCategoryBio.queryHomologs <- unique(full_join(HomologyCategoryBio.queryHomologs,HomologyCategoryBio.queryHomologs.tmp, by=c("HomologyCategory", "biotype_category","orthology_type")));
			#HomologyCategoryBio.queryHomologs.tmp[nrow(HomologyCategoryBio.queryHomologs.tmp)+1,] <- c("total","","", round(colSums(as.numeric(HomologyCategoryBio.queryHomologs.tmp[,c(4:13)])), digits=0))
			
			#stats per type of homology chr
			##############################
			perchrHomologs.tmp <- as.data.frame(unique(df[,c(2,3,6,13,21,8)]) %>% group_by(HomologyCategory,biotype_category,orthology_type,chromosome) %>% count());
			colnames(perchrHomologs.tmp)[c(4,5)] <- c("chromosome", paste0("Nb ",name," homologs to ",specieReferenceName," per chr"));
			orthology_typeOrder <- c("ortholog_one2one","ortholog_one2many","ortholog_many2many");
			perchrHomologs.tmp <- perchrHomologs.tmp[order(perchrHomologs.tmp[,5],decreasing = TRUE ),];
			perchrHomologs.tmp <- perchrHomologs.tmp[order(match(perchrHomologs.tmp[,3],orthology_typeOrder )),];
			perchrHomologs.tmp <- perchrHomologs.tmp[order(match(perchrHomologs.tmp[,2],unique(biotypeCat[,2]) )),];				
			perchrHomologs.tmp <- perchrHomologs.tmp[order(match(perchrHomologs.tmp[,1],homology_typeOrder )),];				
			perchrHomologs <- unique(full_join(perchrHomologs,perchrHomologs.tmp, by=c("HomologyCategory", "biotype_category", "orthology_type","chromosome")));


			perchrHomologsNoOrtCat.tmp <- as.data.frame(unique(df[,c(2,3,6,21,8)]) %>% group_by(HomologyCategory,biotype_category,chromosome) %>% count());
			colnames(perchrHomologsNoOrtCat.tmp)[c(3,4)] <- c("chromosome", paste0("Nb ",name," homologs to ",specieReferenceName," per chr"));
			perchrHomologsNoOrtCat.tmp$total <-length(unique(df[,c(2)])); #total nb of unique reference specie genes
			perchrHomologsNoOrtCat.tmp <- left_join(perchrHomologsNoOrtCat.tmp,HomologyCategoryBio.tmp[,c(1:2,5,3)], by=c("HomologyCategory","biotype_category"));
			colnames(perchrHomologsNoOrtCat.tmp)[5] <- paste0(name, " total number of homologous genes found to ",specieReferenceName);
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[,c(1:3,5,4,6,7)];
			#perchrHomologsNoOrtCat.tmp$perc <- round((perchrHomologsNoOrtCat.tmp[,5]/as.numeric(perchrHomologsNoOrtCat.tmp[,7]))*100,digits=2);
			#colnames(perchrHomologsNoOrtCat.tmp)[8] <- paste0(name, " total number of homologous genes found to ",specieReferenceName);				
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[order(perchrHomologsNoOrtCat.tmp[,5],decreasing = TRUE ),];
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[order(match(perchrHomologsNoOrtCat.tmp[,2],unique(biotypeCat[,2]) )),];				
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[order(match(perchrHomologsNoOrtCat.tmp[,1],homology_typeOrder )),];	
			perchrHomologsNoOrtCat <- unique(full_join(perchrHomologsNoOrtCat,perchrHomologsNoOrtCat.tmp, by=c("HomologyCategory", "biotype_category","chromosome")));		
		}
		 
		i=i+1
	};
	#Saving the results
	resListStats <- list(HomologyCat=HomologyCategory,HomologyCatBio=HomologyCategoryBio,HomCatBio.queryHom=HomologyCategoryBio.queryHomologs,perchrHom=perchrHomologsNoOrtCat, perchrHomologsOrthoType=perchrHomologs)
	write.xlsx(resListStats, file = paste0(wd,f_results,f_homologs, f_stats, specieReferenceName,"_homology_per_labModel_Stats","_",dateHom,".xlsx"));
	assign("Homologs_Stats",resListStats,.GlobalEnv);
	assign("perchrHomologsFList",perchrHomologsFList,.GlobalEnv);
	assign("homology_typeOrder",homology_typeOrder,.GlobalEnv);

};

#run the function
#########################
specieReferenceBiomaRt<- "hsapiens_gene_ensembl"
specieReferenceName <- "Human"
speciesQueryBiomaRt <- c("mmusculus","rnorvegicus","dmelanogaster","ptroglodytes","ggorilla",
	"pabelii","drerio","celegans");
speciesQueryNames <- c("Mouse","Rat","Drosophila","Chimpanzee","Gorilla","Orangutan","Zebrafish","C.elegands");
chr_homologs.function(specieReferenceBiomaRt,specieReferenceName,21,speciesQueryBiomaRt,speciesQueryNames,biotypeCat,"280319")

#data preapration for figures plots:
#####################################
dataFigureHomCat <- Homologs_Stats$HomologyCat
dataFigureHomCat <- dataFigureHomCat[,c(grep("HomologyCategory", colnames(dataFigureHomCat)),grep("Nb Human homologs ", colnames(dataFigureHomCat)),
	grep("% Human", colnames(dataFigureHomCat)))]

dataFigureHom <- Homologs_Stats$HomologyCatBio
dataFigureHom <- dataFigureHom[,c(grep("HomologyCategory", colnames(dataFigureHom)),grep("biotype_category", colnames(dataFigureHom)),grep("Nb Human homologs ", colnames(dataFigureHom)),
	grep("Total nb Human ", colnames(dataFigureHom)),grep("% Human", colnames(dataFigureHom)))]
allHom.tmp <- left_join(dataFigureHom,dataFigureHomCat,by="HomologyCategory")

#1. simple homolog and category sumup
allHom <- dataFigureHomCat[,c(grep("HomologyCategory", colnames(dataFigureHomCat)),grep("% Human homologs ", colnames(dataFigureHomCat)))]
allHom$colPlot <- ifelse(allHom[,1]==homology_typeOrder[1],"#b4b2d7",
		ifelse(allHom[,1]==homology_typeOrder[2],"#efa824","#2de115"
			))
allHom[,1]<-gsub("predicted genes","Predicted genes",gsub("remote homologs","Remote homologs",allHom[,1]))
allHom[,1] <- factor(allHom[,1], levels = unique(allHom[,1]))

#2.with chr info per query specie
plotOrthotypesList<- list()
for (i in 1:length(perchrHomologsFList)){
	df <- perchrHomologsFList[[i]]
	name <- gsub("_.*","",names(perchrHomologsFList)[i])
	df<- unique(df[,c(1:3,8,6,4)])
	df$colPlot <- ifelse(df[,1]==homology_typeOrder[1],"#b4b2d7",
		ifelse(df[,1]==homology_typeOrder[2],"#efa824","#2de115"
			))
	df[,1]<-gsub("High confidence homologs","HCH",gsub("predicted genes","PG",gsub("remote homologs","RH",df[,1])))
	df[,2]<-gsub("Protein coding","PC",gsub("Short noncoding","shortNC",df[,2]))
	df[,3]<-gsub("ortholog_","",df[,3])

	df$category <- paste0(df[,2],":",df[,3],":chr",df[,6],":" )
	df$category <- factor(df$category, levels = unique(df$category))

	plotOrthotypesList[[i]] <- df
	names(plotOrthotypesList)[i] <- name

};

write.xlsx(plotOrthotypesList, file = paste0(wd,f_results,f_homologs, f_stats, specieReferenceName,"_homology_per_labModel_sumupPlot","_",dateHom,".xlsx"));

#donutplot for %categories

catDf$category <- factor(catDf$biotype_category, levels = unique(catDf$biotype_category))

#as there is some superposition, changing the order or categories to plot

pdf(paste0(wd,f_results,f_homologs,f_figures,"_HomologyCategories_donut.pdf"),width = 8.27, height=11.69)  #8.27 × 11.69)

# 3"% Mouse"  
donuts_plot(allHom[,2],rep(1,3), allHom[,1],
        cols=allHom[,10],pilabels=T,legend=T,legend_offset=0.4,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  


plot.new()
donuts_plot(allHom[,2],rep(1,3), paste0(allHom[,2], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  


plot.new()
#"4% Rat" 
donuts_plot(allHom[,3],rep(1,3), paste0(allHom[,3], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  


plot.new()
#"5% Gorilla"
donuts_plot(allHom[,6],rep(1,3), paste0(allHom[,6], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  
 


plot.new()
#"6% Orangutan" 
donuts_plot(allHom[,7],rep(1,3), paste0(allHom[,7], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  



plot.new()
#"7% Zebrafish"
donuts_plot(allHom[,8],rep(1,3), paste0(allHom[,8], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  



plot.new()
#8"% C.elegands"
donuts_plot(allHom[,9],rep(1,3), paste0(allHom[,9], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  



plot.new()
#9"% Drosophila" 
donuts_plot(allHom[,4],rep(1,3), paste0(allHom[,4], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  


plot.new()
#10"% Chimpanzee" 
donuts_plot(allHom[,5],rep(1,3), paste0(allHom[,5], "%"),
        cols=allHom[,10],pilabels=T,legend=F,legend_offset=1.3,
        outradius = 1,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,3) )  


dev.off()

#################################################
#################################################
# Human without homologs in the species
# plot : genomic location of the found homologs
#################################################
#################################################

speciesHomologsFoundList<- list(Gorilla=genesDf_specie_Gorilla,Chimpanzee=genesDf_specie_Chimpanzee,
	Orangutan=genesDf_specie_Orangutan,Mouse=genesDf_specie_Mouse,Rat=genesDf_specie_Rat,
	Zebrafish=genesDf_specie_Zebrafish,Drosophila=genesDf_specie_Drosophila,C.elegands=genesDf_specie_C.elegands)
specieQueryNames <- c("Gorilla","Chimpanzee","Orangutan","Mouse","Rat","Zebrafish","Drosophila","C.elegands");

homologs_no_homologs_found.function <- function(speciesHomologsFoundList,genesDf_ref,specieReferenceBiomaRt,specieQueryNames,date){
	statsList <- list();
	noHomologsList <- list();
	genesDf_ref <- genesDf_ref[order(genesDf_ref$start_position),];
	genesDf_ref <- genesDf_ref[gtools::mixedorder(as.character(genesDf_ref$chromosome_name)),];
	genesDf_ref<- left_join(genesDf_ref,biotypeCat,by="gene_biotype")

	for (i in 1:length(speciesHomologsFoundList)){
		df_specie <- speciesHomologsFoundList[[i]];
		name_specie<- names(speciesHomologsFoundList)[i];

		#reordering main df
	
		genesHomologs <- genesDf_ref[which(genesDf_ref$external_gene_name %in% df_specie$external_gene_name),];
		genesHomologs <- genesHomologs[order(genesHomologs$start_position),];
		genesHomologs <- genesHomologs[gtools::mixedorder(as.character(genesHomologs$chromosome_name)),];
		genesHomologs$specie <- name_specie
		names(genesHomologs)[8]<- paste0("homologs_found in ",name_specie)

		genesNoHomologs <- genesDf_ref[-which(genesDf_ref$external_gene_name %in% df_specie$external_gene_name),];
		genesNoHomologs <- genesNoHomologs[order(genesNoHomologs$start_position),];
		genesNoHomologs <- genesNoHomologs[gtools::mixedorder(as.character(genesNoHomologs$chromosome_name)),];
		
		noHomologsList[[i]]<-genesNoHomologs;
		names(noHomologsList)[i] <- name_specie;
		
		genesHomologsStatsBiotype <- as.data.frame(unique(genesHomologs[,c(2,7)]) %>% group_by(biotype_category) %>% count());
		colnames(genesHomologsStatsBiotype)[2]<- "Nb total genes found homolog per biotype category";

		genesNoHomologsStatsBiotype <- as.data.frame(unique(genesNoHomologs[,c(2,7)]) %>% group_by(biotype_category) %>% count());
		colnames(genesNoHomologsStatsBiotype)[2]<- "Nb total genes without homolog per biotype category";
		
		stats <- full_join(genesHomologsStatsBiotype,genesNoHomologsStatsBiotype,by=colnames(genesHomologsStatsBiotype)[1]);
		stats$totalHumanGenesChr21<- length(unique(genesDf_ref[,2])); #815
		stats$`Nb total homologs found`<-  length(unique(genesDf_ref[which(genesDf_ref$external_gene_name %in% df_specie$external_gene_name),2])); #208
		stats$`Nb total genes without homolog`<-  length(unique(genesNoHomologs[,2])); #208
		stats<- stats[,c(1,4:6,2:3)]
		#results
		statsList[[i]] <- stats
		names(statsList)[i] <- name_specie
				
		
		if(i==1){
			HomologsFound.NotFound<- left_join(genesDf_ref,genesHomologs,by=colnames(genesDf_ref))

		} else {
			HomologsFound.NotFound.tmp<- left_join(genesDf_ref,genesHomologs,by=colnames(genesDf_ref))
			HomologsFound.NotFound<- left_join(HomologsFound.NotFound,HomologsFound.NotFound.tmp,by=colnames(genesDf_ref))

		}
		i=i+1
	};
	assign("noHomologsList",noHomologsList,.GlobalEnv);
	assign("statsList_homologs",statsList,.GlobalEnv);
	assign("HomologsFound.NotFound",HomologsFound.NotFound,.GlobalEnv);
	assign("genesDf_ref",genesDf_ref,.GlobalEnv);
	write.xlsx(statsList, file = paste0(wd,f_results,f_homologs,f_tables, specieReferenceName, "_specieRef_homologs_no_homologs_defined_by_ensembl_",date, ".xlsx"));
};
date<- "310319"
homologs_no_homologs_found.function(speciesHomologsFoundList,genesDf_ref,specieReferenceBiomaRt,specieQueryNames,date)



####################################
# Visualization of the homologs 
# respect to the genomic  human 
# chr21 location
####################################


####################################
# A. If we dont need to split the data
####################################
	GenPlot <-unique(HomologsFound.NotFound[,c(2,4,7:15)]);
	GenPlot<- GenPlot[!duplicated(GenPlot[,1]),];

	#adding colour code per biotype category as defined before:
	GenPlot$biotype_category <- factor(GenPlot$biotype_category, levels = catDf$biotype_category )
	GenPlot$colPlotCat <- ifelse(GenPlot[,3]==levels(GenPlot$biotype_category)[1],"turquoise3",
			ifelse(GenPlot[,3]==levels(GenPlot$biotype_category)[2],"khaki",
			ifelse(GenPlot[,3]==levels(GenPlot$biotype_category)[3],"purple",
			ifelse(GenPlot[,3]==levels(GenPlot$biotype_category)[4],"deeppink3",
			ifelse(GenPlot[,3]==levels(GenPlot$biotype_category)[5],"orange",
				NA)))))
	colDonutCat<- c("#00C5CD","#FFFF99","#993399","#CC0066","#FF9933")
	GenPlot$human <- "Human"
	#GenPlot$external_gene_name <- factor(GenPlot$external_gene_name, levels = unique(GenPlot$external_gene_name ))

	GenPlot$colPlotCat2<- ifelse(GenPlot[,3]==levels(catDf$biotype_category)[1],"#00C5CD",
		ifelse(GenPlot[,3]==levels(catDf$biotype_category)[2],"#cccc7a",
		ifelse(GenPlot[,3]==levels(catDf$biotype_category)[3],"#993399",
		ifelse(GenPlot[,3]==levels(catDf$biotype_category)[4],"#CC0066",
		ifelse(GenPlot[,3]==levels(catDf$biotype_category)[5],"#FF9933",
			NA)))));


	groups_Allsumup <- factor(GenPlot$colPlotCat2, levels = c("#00C5CD","#cccc7a","#993399","#CC0066","#FF9933")) #100
	GenPlot[,13] <-  as.numeric(gsub("Human","9",GenPlot[,13]))
	GenPlot[,4] <-  as.numeric(gsub("Gorilla","7",GenPlot[,4]))
	GenPlot[,5] <-  as.numeric(gsub("Chimpanzee","8",GenPlot[,5]))
	GenPlot[,6] <-  as.numeric(gsub("Orangutan","6",GenPlot[,6]))
	GenPlot[,7] <-  as.numeric(gsub("Mouse","4",GenPlot[,7]))
	GenPlot[,8] <-  as.numeric(gsub("Rat","5",GenPlot[,8]))
	GenPlot[,9] <-  as.numeric(gsub("Zebrafish","3",GenPlot[,9]))
	GenPlot[,10] <-  as.numeric(gsub("Drosophila","2",GenPlot[,10]))
	GenPlot[,11] <-  as.numeric(gsub("C.elegands","1",GenPlot[,11]))
	bgin_pos <- as.numeric(rownames(GenPlot)[1])
	end_pos <- as.numeric(length(rownames(GenPlot)))


#################################################
#  Plor R basic graphics, only BLACK COLOR!!  #########
#################################################

#################################################
#now the strategy to represent each in one plot 
# will be to give them an y value of 1 to n depending of the model
# data represented ordered by decrease number of homologs genes
#^
#9. human
#8. chimpanzee
#7. gorilla
#6. orangutan
#5 rat
#4 mouse
#3. zebrafish
#2.Fruitfly
#1. C.elegans
#
#y-axis
#
# x-axis Ordered gene names by genomic coordinates 
# in chr21 >
#################################################

library("extrafont")
font_import() # import all your fonts
fonts() #get a list of fonts
fonttable()


pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologsScatterplot_with_GeneNames2.pdf"), width=48.69, height=9.27)
labels <- c("Human","Chimpanzee","Gorilla","Orangutan","Rat","Mouse","Zebrafish","Drosophila","C.elegands")
op <- par(las=1,mar=c(12, 15, 4, 2))
#plot(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,13], col=c("turquoise3","khaki1","purple","deeppink","darkorange1")[GenPlot$biotype_category], type="p", pch=15,  frame=TRUE, cex=12, las=2, xaxt = "n",yaxt = "n",cex.lab=1.5,cex.axis=1.3, ylab="", xlab="",  ylim=c(1,9))
#plot(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,13], col=factor(GenPlot$colPlotCat2), type="p", pch=15,  frame=TRUE, cex=12, las=2, xaxt = "n",yaxt = "n",cex.lab=1.5,cex.axis=1.3, ylab="", xlab="",  ylim=c(1,9))
#plot(GenPlot[,2]), GenPlot[,13], col="red",   pch=15,  frame=TRUE, cex=12, las=2, xaxt = "n",yaxt = "n",cex.lab=1.5,cex.axis=1.3, ylab="", xlab="",  ylim=c(1,9))
plot(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,13], col=c("turquoise3","khaki1","purple","deeppink","darkorange1")[groups_Allsumup],  pch=15,  frame=TRUE, cex=12, las=2, xaxt = "n",yaxt = "n",cex.lab=1.5,cex.axis=1.3, ylab="", xlab="",  ylim=c(1,9))
title(ylab="Specie", line=6, cex.lab=1.5, family="Arial")
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,5], col="black", type="p", pch=15,  cex=0.54, las=2,xaxt = "n", yaxt = "n",cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,4], col="black",type="p",  pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,6], col="black",type="p",  pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,8], col="black",type="p",  pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,7], col="black", type="p", pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,9], col="black", type="p", pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,10], col="black",type="p",  pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
points(factor(GenPlot[,1], levels=GenPlot[,1]), GenPlot[,11], col="black",type="p",  pch=15,  cex=0.54, las=2,xaxt = "n",yaxt = "n", cex.lab=1,cex.axis=0.7, ylab="", xlab="",  ylim=c(1,9))
axis(side=2, at=seq(1, 9, by=1), labels = rev(labels))
abline(v=bgin_pos- 0.5,  col= "black",lty = 18)
abline(v=end_pos +0.5,  col= "black",lty = 18)
op <- par(cex=0.3,mgp=c(4,1.2,0))#,las=2
title(main="Homologs to human chr21", outer=TRUE,line = -2, cex.main=3.4,family="Arial") #
axis(side=1, at=1:length(GenPlot[,1]), tick=TRUE, las=2, labels=GenPlot[,1],tck=-0.03);
op <- par(cex=1.4);
dev.off()

#################################################
#  ggplot option #########
#################################################

#plot 1:  BY GENOMIC LCOATION  SHOWING GENE NAME 

	legend_title <- "Biotype Category"
	myTitlePlot <- "Homologs to human chr21"

	plotFC_bio <- ggplot(GenPlot, aes(x=GenPlot[,1], y=value,color =groups_Allsumup)) +
		geom_point(aes(y=GenPlot[,13]), size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,5]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,4]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,6]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,8]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,7]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,9]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,10]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,11]),size = 1.6,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain"),plot.margin=unit(c(2,1,1,4),"cm")) + guides(color = guide_legend(override.aes = list(size=5))) 
			


plotFC_bio <-plotFC_bio + scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9),labels = rev(labels)) + theme(axis.text.y = element_text(size = 12, face = "plain"))

#plot 2:  BY GENOMIC LCOATION NOT SHOWING GENE NAME INSTEAD SHOWING GENOMIC COORDINATES
options(scipen=10000)
	plotFC_bio_loc <- ggplot(GenPlot, aes(x=GenPlot[,2], y=value,color =groups_Allsumup)) +
		geom_point(aes(y=GenPlot[,13]), size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,5]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,4]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,6]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,8]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,7]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,9]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,10]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=GenPlot[,11]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_Allsumup)) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) 
			


plotFC_bio_loc <-plotFC_bio_loc + scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9),labels = rev(labels)) + theme(axis.text.y = element_text(size = 12, face = "plain"),axis.text.x = element_text(size = 12, face = "plain"))



pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologsScatterplot_with_GeneNamesAll.pdf"), width=60.69, height=9.27)
print(plotFC_bio)
plot.new()
print(plotFC_bio_loc)
dev.off()


pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologsScatterplot_with_GeneNames_byGenomicLocation.pdf"), width=18.69, height=9.27)
print(plotFC_bio_loc)
dev.off()


####################################
# B. If we need to split the data
####################################
formatting_plotData_Splitted.function<- function(data,catDf,nbParts,date){

	GenPlotF <-unique(data[,c(1,2,4:15)]);
	GenPlot <-unique(data[,c(2,4,7:15)]);
	GenPlot<- GenPlot[!duplicated(GenPlot[,1]),];

	genPLotList <- list();
	bgin_posList <- list();
	end_posList <- list(); 
	begin_pos_locList <- list() ;  
	end_pos_locList <- list();
	groups_sumupList<- list();
	lol=1
	for (i in 1:nbParts){		
		n=round(length(GenPlot[,1])/nbParts)*i;
		df <- GenPlot[lol:n,];
		rownames(df) <- 1:nrow(df);
		#adding colour code per biotype category as defined before:
		df$biotype_category <- factor(df$biotype_category, levels = catDf$biotype_category )
		df$colPlotCat<- ifelse(df[,3]==levels(catDf$biotype_category)[1],"turquoise3",
				ifelse(df[,3]==levels(catDf$biotype_category)[2],"khaki",
				ifelse(df[,3]==levels(catDf$biotype_category)[3],"purple",
				ifelse(df[,3]==levels(catDf$biotype_category)[4],"deeppink3",
				ifelse(df[,3]==levels(catDf$biotype_category)[5],"orange",
					NA)))));
		#colDonutCat<- c("#00C5CD","#FFFF99","#993399","#CC0066","#FF9933")
		df$human <- "Human";
		#df$external_gene_name <- factor(df$external_gene_name, levels = unique(df$external_gene_name ))
		df$colPlotCat2<- ifelse(df[,3]==levels(catDf$biotype_category)[1],"#00C5CD",
				ifelse(df[,3]==levels(catDf$biotype_category)[2],"#cccc7a",
				ifelse(df[,3]==levels(catDf$biotype_category)[3],"#993399",
				ifelse(df[,3]==levels(catDf$biotype_category)[4],"#CC0066",
				ifelse(df[,3]==levels(catDf$biotype_category)[5],"#FF9933",
					NA)))));

		df[,13] <-  as.numeric(gsub("Human","4.5",df[,13]));
		df[,4] <-  as.numeric(gsub("Gorilla","3.5",df[,4]));
		df[,5] <-  as.numeric(gsub("Chimpanzee","4",df[,5]));
		df[,6] <-  as.numeric(gsub("Orangutan","3",df[,6]));
		df[,7] <-  as.numeric(gsub("Mouse","2",df[,7]));
		df[,8] <-  as.numeric(gsub("Rat","2.5",df[,8]));
		df[,9] <-  as.numeric(gsub("Zebrafish","1.5",df[,9]));
		df[,10] <-  as.numeric(gsub("Drosophila","1",df[,10]));
		df[,11] <-  as.numeric(gsub("C.elegands","0.5",df[,11]));
		df <- df[!is.na(df[,14]),]
		bgin_posList[[i]] <- as.numeric(rownames(df)[1]);
		end_posList[[i]] <- as.numeric(length(rownames(df)));
		begin_pos_locList[[i]] <- df[bgin_pos,'start_position'];
		end_pos_locList[[i]] <- df[end_pos,'end_position'];
		groups_sumupList[[i]] <- factor(df$colPlotCat2, levels =c("#00C5CD","#cccc7a","#993399","#CC0066","#FF9933")); 
		genPLotList[[i]] <- df
		lol <- n +1;	
		i=i+1
	};

	assign("bgin_posList",bgin_posList,.GlobalEnv);
	assign("end_posList",end_posList,.GlobalEnv);
	assign("begin_pos_locList",begin_pos_locList,.GlobalEnv);
	assign("end_pos_locList",end_pos_locList,.GlobalEnv);
	assign("genPLotList",genPLotList,.GlobalEnv);
	assign("GenPlot",GenPlot,.GlobalEnv);
	assign("groups_sumupList",groups_sumupList,.GlobalEnv);

	#processing the data once more to produce a table to include in our results for the book
	GenPlotF[,7] <-  gsub("Gorilla","Yes",GenPlotF[,7]);
	GenPlotF[,8] <-  gsub("Chimpanzee","Yes",GenPlotF[,8]);
	GenPlotF[,9] <-  gsub("Orangutan","Yes",GenPlotF[,9]);
	GenPlotF[,10] <-  gsub("Mouse","Yes",GenPlotF[,10]);
	GenPlotF[,11] <-  gsub("Rat","Yes",GenPlotF[,11]);
	GenPlotF[,12] <-  gsub("Zebrafish","Yes",GenPlotF[,12]);
	GenPlotF[,13] <-  gsub("Drosophila","Yes",GenPlotF[,13]);
	GenPlotF[,14] <-  gsub("C.elegands","Yes",GenPlotF[,14]);
	assign("GenPlotF",GenPlotF,.GlobalEnv);
	write.xlsx(GenPlotF, file = paste0(wd,f_results,f_homologs,f_tables, specieReferenceName, "_specieRef_genes_homologs_found_BookTable_ensembl_",date, ".xlsx"));
};

formatting_plotData_Splitted.function(HomologsFound.NotFound,catDf,8,"310319")

#################################################
#  Plor R basic graphics, only BLACK COLOR!!  #########
#################################################

#now the strategy to represent each in one plot will be to give them an y value of 1 to n depending of the model
#data represented ordered by decrease number of homologs genes
#^
#9. human
#8. chimpanzee
#7. gorilla
#6. orangutan
#5 rat
#4 mouse
#3. zebrafish
#2.Fruitfly
#1. C.elegans
#
#y-axis
#
#x-axis Ordered gene names by genomic coordinates in chr21 >

library("extrafont")
font_import() # import all your fonts
fonts() #get a list of fonts
fonttable()
library("gridExtra")
library("cowplot")
library("ggpubr")
 library("grid")

#################################################
#  ggplot option #########
#################################################

#plot 1:  BY GENOMIC LCOATION  SHOWING GENE NAME 

legend_title <- "Biotype Category"
myTitlePlot <- "Homologs to human chr21 in DS animal models"

plotFC_bio1 <- ggplot(genPLotList[[1]], aes(x=genPLotList[[1]][,1], y=value,color =groups_sumupList[[1]])) +
	geom_point(aes(y=genPLotList[[1]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[1]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[1]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio1 <-plotFC_bio1 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) 

#2
plotFC_bio2 <- ggplot(genPLotList[[2]], aes(x=genPLotList[[2]][,1], y=value,color =groups_sumupList[[2]])) +
	geom_point(aes(y=genPLotList[[2]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[2]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[2]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio2 <-plotFC_bio2 + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) + theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")

#3
plotFC_bio3 <- ggplot(genPLotList[[3]], aes(x=genPLotList[[3]][,1], y=value,color =groups_sumupList[[3]])) +
	geom_point(aes(y=genPLotList[[3]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[3]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[3]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio3 <-plotFC_bio3 + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) + theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")
#4
plotFC_bio4 <- ggplot(genPLotList[[4]], aes(x=genPLotList[[4]][,1], y=value,color =groups_sumupList[[4]])) +
	geom_point(aes(y=genPLotList[[4]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[4]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[4]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio4 <-plotFC_bio4 + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) + theme(plot.title = element_text(colour="white"), axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")

#5

plotFC_bio5 <- ggplot(genPLotList[[5]], aes(x=genPLotList[[5]][,1], y=value,color =groups_sumupList[[5]])) +
	geom_point(aes(y=genPLotList[[5]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[5]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[5]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[5]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio5 <-plotFC_bio5 + theme(axis.text.y = element_text(size = 7, face = "bold"))  + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels))+ theme(legend.position="none") 

#6
plotFC_bio6 <- ggplot(genPLotList[[6]], aes(x=genPLotList[[6]][,1], y=value,color =groups_sumupList[[6]])) +
	geom_point(aes(y=genPLotList[[6]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[6]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[6]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[6]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio6 <-plotFC_bio6 + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) + theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")

#7
plotFC_bio7 <- ggplot(genPLotList[[7]], aes(x=genPLotList[[7]][,1], y=value,color =groups_sumupList[[7]])) +
	geom_point(aes(y=genPLotList[[7]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,5]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[7]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,6]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,8]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,7]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,9]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,10]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[7]][,11]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[7]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio7 <-plotFC_bio7 + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) + theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")
#8
plotFC_bio8 <- ggplot(genPLotList[[8]], aes(x=genPLotList[[8]][,1], y=value,color=groups_sumupList[[8]])) +
	geom_point(aes(y=genPLotList[[8]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,5]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +		
	geom_point(aes(y=genPLotList[[8]][,4]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,6]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,8]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,7]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]]))+
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,9]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]]))+
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,10]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList[[8]][,11]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[8]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 			
plotFC_bio8 <-plotFC_bio8 + scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5),labels = rev(labels)) + theme(plot.title = element_text(colour="white"), axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")




plotFC_bio9 <- ggplot(genPLotList[[1]], aes(x=genPLotList[[1]][,1], y=value,color=groups_sumupList[[1]])) +
	geom_point(aes(y=genPLotList[[1]][,13]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(0, 0.7, 0, 0.7), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=11), axis.text=element_text(), axis.text.x = element_text(),axis.text.y = element_text()) + guides(color = guide_legend(override.aes = list(size=5))) 

legend <- cowplot::get_legend(plotFC_bio9)


pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologsScatterplot_Book.pdf"), width=11.69, height=8.27)
ggarrange(plotFC_bio1, plotFC_bio2,plotFC_bio3, plotFC_bio4,
	ncol = 1, nrow = 4,  align = "hv" )
ggarrange(plotFC_bio5, plotFC_bio6, plotFC_bio7, plotFC_bio8,
	ncol = 1, nrow = 4,  align = "hv" )

grid.newpage()
grid.draw(legend)

dev.off()

#cuack <- ggarrange(plotFC_bio1, plotFC_bio2,plotFC_bio3, plotFC_bio4, common.legend = TRUE, legend = "bottom")

#plot 2:  BY GENOMIC LCOATION NOT SHOWING GENE NAME INSTEAD SHOWING GENOMIC COORDINATES
options(scipen=10000)
	plotFC_bio_loc <- ggplot(genPLotList[[1]], aes(x=genPLotList[[1]][,2], y=value,color =groups_sumupList[[1]])) +
		geom_point(aes(y=genPLotList[[1]][,13]), size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=7), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,5]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,4]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,6]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,8]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,7]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,9]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,10]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) +
		geom_point(aes(y=genPLotList[[1]][,11]),size = 1,shape = 15) + theme(plot.title = element_text(size = 30, face = "bold"), legend.title=element_text(size=17,face="bold"),axis.text=element_text(size=1), legend.text=element_text(size=9), axis.text.x = element_text(angle=60, hjust=0,size = 1, face = "plain")) + 
	labs( title= myTitlePlot)+ ylab("Species") + xlab("") + theme_gray(base_size = 14)+
	scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList[[1]])) +
	theme(plot.title = element_text(hjust = 0.5)) +ggpubr::rotate_x_text()+ theme(axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 4, face = "plain")) + guides(color = guide_legend(override.aes = list(size=5))) 
			


plotFC_bio_loc <-plotFC_bio_loc + scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9),labels = rev(labels)) + theme(axis.text.y = element_text(size = 12, face = "plain"),axis.text.x = element_text(size = 12, face = "plain"))



pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologsScatterplot_with_GeneNamesAll_cuack.pdf"), width=11.69, height=8.27)
print(plotFC_bio1)
dev.off()


pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologsScatterplot_with_GeneNames_byGenomicLocation.pdf"), width=18.69, height=9.27)
print(plotFC_bio_loc)
dev.off()


###############################################################
#
# Represntin the plot of homologous from the point of view of 
#    the genomic coordinates of mouse and rats chr
#
###############################################################

#1. creating the genesDf_Ref objects for mus and rats
# mouse: chr 10,16,17
# rat: chr 11,20


##################################################################################
#A. GENOME WIDE ANALYSIS
##################################################################################
speciesBiomaRt <- c("mmusculus_gene_ensembl","rnorvegicus_gene_ensembl");
speciesNames <- c("Mouse","Rat");
chrMouseSel <- c("10","16","17");
chrRattuSel <- c("11","20");
chrListSel <- list(chrMouse=chrMouseSel,chrRattu=chrRattuSel)

mouseHomologsHumanChr21 <- filter(GenPlotF,`homologs_found in Mouse` =="Yes")[,c(2,10)]
colnames(mouseHomologsHumanChr21)<-c("hsapiens_homolog_associated_gene_name", "homologs_found in Human")  
ratHomologsHumanChr21 <- filter(GenPlotF,`homologs_found in Rat` =="Yes")[,c(2,11)]
colnames(ratHomologsHumanChr21)<-c("hsapiens_homolog_associated_gene_name", "homologs_found in Human")  
chr21Homologs <- list(mouse=mouseHomologsHumanChr21, rat=ratHomologsHumanChr21)

chr_refGenes.function <- function(speciesBiomaRt,speciesNames,chrListSel,biotypeCat,chr21Homologs){
	##################################################################################
	# Function to annotate in the syntenic chrs of human chr21 in other species
	#  all the animal species genes, annotated with if the human homologues is found or not there.
	# its homologs in the species of interest defined as query
	# Its produce the stats. Input files are:
	#
	##################################################################################	

	#2. Quering ENSEMBL 
	######################
	#ensembl <- useEnsembl(biomart = "ensembl")
	#listDatasets(mart = ensembl)
	#for the reference specie:
	for (i in 1:length(chrListSel)){
		chrSelSpecie <- chrListSel[[i]]
		chr21Hom<- chr21Homologs[[i]]
		ensembl = useMart("ensembl", dataset = speciesBiomaRt[i]); 
		attributes = listAttributes(ensembl);
		attributesHom = attributes[which(attributes[,3]=="homologs"),];	
		rownames(attributesHom) <- 1:nrow(attributesHom);

		#to get basic info plus homologous in human no biotype is allowed
		attribRef <-c(as.numeric(rownames(attributesHom[which(attributesHom[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributesHom[which(attributesHom[,1]=='external_gene_name'),][1,])),
			as.numeric(rownames(attributesHom[which(attributesHom[,1]=='chromosome_name'),][1,])),as.numeric(rownames(attributesHom[which(attributesHom[,1]=='start_position'),][1,])),as.numeric(rownames(attributesHom[which(attributesHom[,1]=='end_position'),][1,])),
			as.numeric(rownames(attributesHom[which(attributesHom[,1]=='hsapiens_homolog_associated_gene_name'),][1,]))   );
		attribRef <-attributesHom[attribRef,1];  
		genesDf_ref = getBM(attributes=attribRef,filters="chromosome_name",values=chrSelSpecie, mart=ensembl); 
		#to get the biotype
		attribBio <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])), as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,]))   );
		attribBio <-attributes[attribBio,1];  
		genesDf_Bio = getBM(attributes=attribBio,filters="chromosome_name",values=chrSelSpecie, mart=ensembl); 

		genesDf_ref<- left_join(genesDf_ref,genesDf_Bio,by="ensembl_gene_id")
		genesDf_ref<-genesDf_ref[!duplicated(genesDf_ref$ensembl_gene_id), ]; 
		
		genesDf_Dupl<-genesDf_ref[duplicated(genesDf_ref$external_gene_name), ]; 
		genesDf_ref<-genesDf_ref[!duplicated(genesDf_ref$external_gene_name), ]; 
		#print(length(genesDf_Dupl));
		#print(genesDf_Dupl);
		genesDf_ref<- left_join(genesDf_ref,biotypeCat,by="gene_biotype")
		genesDf_ref <- genesDf_ref[order(match(genesDf_ref$chromosome_name,chrSelSpecie )),];
		genesDf_ref <- genesDf_ref[order(genesDf_ref$start_position),];
		genesDf_ref <-genesDf_ref[gtools::mixedorder(as.character(genesDf_ref$chromosome_name)),];	

		#joining with human ref data
		#humanData <- data.frame(hsapiens_homolog_associated_gene_name=unique(genesDf_ref[,c(2)]))
		#humanData$`homologs found in Human` <- "Human"
		genesDf_ref<- left_join(genesDf_ref,chr21Hom,by="hsapiens_homolog_associated_gene_name")
		genesDf_ref<-genesDf_ref[!duplicated(genesDf_ref$external_gene_name), ]; 
		genesDf_ref <-genesDf_ref[,c(1:5,7,8,9)]

		assign(paste0("genesChrSelDf_", speciesNames[i]), genesDf_ref,.GlobalEnv);
		assign(paste0("removed_dup_geneNames_genesChrSelDf_", speciesNames[i]), genesDf_Dupl,.GlobalEnv);
		i+i+1
	}
}
chr_refGenes.function(speciesBiomaRt,speciesNames,chrListSel,biotypeCat,chr21Homologs)


#preparing the data to be plotted:




####################################
# B. If we need to split the data
####################################
genesChrSelDf_Mouse_chr10<- filter(genesChrSelDf_Mouse,chromosome_name=="10" )
genesChrSelDf_Mouse_chr16<- filter(genesChrSelDf_Mouse, chromosome_name=="16")
genesChrSelDf_Mouse_chr17<- filter(genesChrSelDf_Mouse,chromosome_name=="17" )
genesChrSelDf_Rat_chr11 <- filter(genesChrSelDf_Rat, chromosome_name=="11")
genesChrSelDf_Rat_chr20 <- filter(genesChrSelDf_Rat, chromosome_name=="20")

dim(genesChrSelDf_Mouse_chr10)
dim(genesChrSelDf_Mouse_chr16)
dim(genesChrSelDf_Mouse_chr17)
dim(genesChrSelDf_Rat_chr11) 
dim(genesChrSelDf_Rat_chr20 )

# to plot depending on where are the homologus genes,
# as they are highly conserved we dont need to plot all mouse and rat chr
genesChrSelDf_Rat_chr11Plot <-genesChrSelDf_Rat_chr11[1:240,] #split in 2
rownames(genesChrSelDf_Rat_chr11Plot) <-1:nrow(genesChrSelDf_Rat_chr11Plot)
genesChrSelDf_Rat_chr20Plot <-genesChrSelDf_Rat_chr20[c(428:525,836:856),]
genesChrSelDf_Rat_chr20Plot[which(rownames(genesChrSelDf_Rat_chr20Plot)==836),2]<- " "
rownames(genesChrSelDf_Rat_chr20Plot) <-1:nrow(genesChrSelDf_Rat_chr20Plot)

genesChrSelDf_Mouse_chr10Plot <-genesChrSelDf_Mouse_chr10[c(1224:1336),]
rownames(genesChrSelDf_Mouse_chr10Plot) <-1:nrow(genesChrSelDf_Mouse_chr10Plot)

genesChrSelDf_Mouse_chr16Plot <-genesChrSelDf_Mouse_chr16[c(1200:1603),] #split in 2
rownames(genesChrSelDf_Mouse_chr16Plot) <-1:nrow(genesChrSelDf_Mouse_chr16Plot)

genesChrSelDf_Mouse_chr17Plot <-genesChrSelDf_Mouse_chr17[c(81:86,180:197,865:916,2017:2056),]
genesChrSelDf_Mouse_chr17Plot[which(rownames(genesChrSelDf_Mouse_chr17Plot)==2017),2]<- " "
genesChrSelDf_Mouse_chr17Plot[which(rownames(genesChrSelDf_Mouse_chr17Plot)==865),2]<- " "
genesChrSelDf_Mouse_chr17Plot[which(rownames(genesChrSelDf_Mouse_chr17Plot)==180),2]<- " "
genesChrSelDf_Mouse_chr17Plot[which(rownames(genesChrSelDf_Mouse_chr17Plot)==81),2]<- " "
rownames(genesChrSelDf_Mouse_chr17Plot) <-1:nrow(genesChrSelDf_Mouse_chr17Plot)


formatting_plotData_Splitted.function<- function(data,catDf,nbParts,specieName, chrName,date){
	# to decide in how many genes split the data for row, dont go over 105-115 genes
	GenPlotF <-unique(data[,c(1,2,4:8)]);
	GenPlot <-unique(data[,c(2,4,7:8)]);
	GenPlot<- GenPlot[!duplicated(GenPlot[,1]),];

	genPLotList <- list();
	bgin_posList <- list();
	end_posList <- list(); 
	begin_pos_locList <- list() ;  
	end_pos_locList <- list();
	groups_sumupList<- list();
	lol=1
	for (i in 1:nbParts){		
		n=round(length(GenPlot[,1])/nbParts)*i;
		df <- GenPlot[lol:n,];
		rownames(df) <- 1:nrow(df);
		number<- as.numeric(nrow(df))
		#adding colour code per biotype category as defined before:
		df[1,1] <- paste0(df[1,2]," bps: ", df[1,1])
		df[nrow(df),1] <- paste0(df[nrow(df),2], " bps: ",df[nrow(df),1] )
		df$external_gene_name <- factor(df$external_gene_name, levels = unique(df$external_gene_name) )
		df$biotype_category <- factor(df$biotype_category, levels = catDf$biotype_category )
		df$colPlotCat<- ifelse(df[,3]==levels(catDf$biotype_category)[1],"turquoise3",
				ifelse(df[,3]==levels(catDf$biotype_category)[2],"khaki",
				ifelse(df[,3]==levels(catDf$biotype_category)[3],"purple",
				ifelse(df[,3]==levels(catDf$biotype_category)[4],"deeppink3",
				ifelse(df[,3]==levels(catDf$biotype_category)[5],"orange",
					NA)))));
		#colDonutCat<- c("#00C5CD","#FFFF99","#993399","#CC0066","#FF9933")
		df$specie <- specieName;
		colnames(df)[6] <- specieName
		colnames(df)[4] <- "Human"
		#df$external_gene_name <- factor(df$external_gene_name, levels = unique(df$external_gene_name ))
		df$colPlotCat2<- ifelse(df[,3]==levels(catDf$biotype_category)[1],"#00C5CD",
				ifelse(df[,3]==levels(catDf$biotype_category)[2],"#cccc7a",
				ifelse(df[,3]==levels(catDf$biotype_category)[3],"#993399",
				ifelse(df[,3]==levels(catDf$biotype_category)[4],"#CC0066",
				ifelse(df[,3]==levels(catDf$biotype_category)[5],"#FF9933",
					NA)))));
		
		df[,4] <-  as.numeric(gsub("Yes","0.5",df[,4]));
		df[,6] <-  as.numeric(gsub(specieName,"2",df[,6]));
		#df <- df[!is.na(df[,14]),]
		bgin_posList[[i]] <- as.numeric(rownames(df)[1]);
		end_posList[[i]] <- as.numeric(length(rownames(df)));
		begin_pos_locList[[i]] <- df[bgin_pos,'start_position'];
		end_pos_locList[[i]] <- df[end_pos,'end_position'];
		groups_sumupList[[i]] <- factor(df$colPlotCat2, levels =c("#00C5CD","#cccc7a","#993399","#CC0066","#FF9933")); 
		genPLotList[[i]] <- df
		lol <- n +1;	
		i=i+1
	};

	assign(paste0("bgin_posList_",specieName,"_",chrName),bgin_posList,.GlobalEnv);
	assign(paste0("end_posList_",specieName,"_",chrName),end_posList,.GlobalEnv);
	assign(paste0("begin_pos_locList_",specieName,"_",chrName),begin_pos_locList,.GlobalEnv);
	assign(paste0("end_pos_locList_",specieName,"_",chrName),end_pos_locList,.GlobalEnv);
	assign(paste0("genPLotList_",specieName,"_",chrName),genPLotList,.GlobalEnv);
	assign(paste0("GenPlot_",specieName,"_",chrName),GenPlot,.GlobalEnv);
	assign(paste0("groups_sumupList_",specieName,"_",chrName),groups_sumupList,.GlobalEnv);

	#processing the data once more to produce a table to include in our results for the book
	GenPlotF[,4] <-  gsub("Human","Yes",GenPlotF[,4]);
	assign(paste0("GenPlotF_",specieName,"_",chrName),GenPlotF,.GlobalEnv);
	write.xlsx(GenPlotF, file = paste0(wd,f_results,f_homologs,f_tables, specieName,"_",chrName, "_syntenicChrs_to_human_chr21_ensembl_",date, ".xlsx"));
};



formatting_plotData_Splitted.function(genesChrSelDf_Mouse_chr16Plot,catDf,4,"Mouse","chr16","290419")
formatting_plotData_Splitted.function(genesChrSelDf_Mouse_chr10Plot,catDf,1,"Mouse","chr10","290419")
formatting_plotData_Splitted.function(genesChrSelDf_Mouse_chr17Plot,catDf,1,"Mouse","chr17","290419")
formatting_plotData_Splitted.function(genesChrSelDf_Rat_chr20Plot,catDf,1,"Rat","chr20","290419")
formatting_plotData_Splitted.function(genesChrSelDf_Rat_chr11Plot,catDf,2,"Rat","chr11","290419")


# further formatting the data of mouse chr17 and rat chr20 to show easily the genomic breaks in the plots:

genPLotList_Rat_chr20 <- genPLotList_Rat_chr20[[1]]
genPLotList_Rat_chr20[,1] <- as.character(genPLotList_Rat_chr20[,1])
genPLotList_Rat_chr20[98,1] <-paste0(genPLotList_Rat_chr20[98,2]," bps: ", genPLotList_Rat_chr20[98,1]) 
genPLotList_Rat_chr20[100,1] <- paste0(genPLotList_Rat_chr20[100,2]," bps: ", genPLotList_Rat_chr20[100,1])
genPLotList_Rat_chr20[99,c(2,3,4,5,6,7)] <- NA
genPLotList_Rat_chr20[99,1] <-" "
genPLotList_Rat_chr20$external_gene_name <- factor(genPLotList_Rat_chr20$external_gene_name, levels = unique(genPLotList_Rat_chr20$external_gene_name) )
genPLotList_Rat_chr20 <- list(genPLotList_Rat_chr20=genPLotList_Rat_chr20)





genPLotList_Mouse_chr17 <- genPLotList_Mouse_chr17[[1]]
genPLotList_Mouse_chr17[,1] <- as.character(genPLotList_Mouse_chr17[,1])
genPLotList_Mouse_chr17[74,1] <-paste0(genPLotList_Mouse_chr17[74,2]," bps: ", genPLotList_Mouse_chr17[74,1]) 
genPLotList_Mouse_chr17[76,1] <- paste0(genPLotList_Mouse_chr17[76,2]," bps: ", genPLotList_Mouse_chr17[76,1])
genPLotList_Mouse_chr17[25,1] <-paste0(genPLotList_Mouse_chr17[25,2]," bps: ", genPLotList_Mouse_chr17[25,1]) 
genPLotList_Mouse_chr17[23,1] <- paste0(genPLotList_Mouse_chr17[23,2]," bps: ", genPLotList_Mouse_chr17[23,1])
genPLotList_Mouse_chr17[6,1] <-paste0(genPLotList_Mouse_chr17[6,2]," bps: ", genPLotList_Mouse_chr17[6,1]) 
genPLotList_Mouse_chr17[8,1] <- paste0(genPLotList_Mouse_chr17[8,2]," bps: ", genPLotList_Mouse_chr17[8,1])
genPLotList_Mouse_chr17[c(75,24,7),c(2,3,4,5,6,7)] <- NA
genPLotList_Mouse_chr17[7,1] <-" "
genPLotList_Mouse_chr17[75,1] <-"  "
genPLotList_Mouse_chr17[24,1] <-"   "

genPLotList_Mouse_chr17$external_gene_name <- factor(genPLotList_Mouse_chr17$external_gene_name, levels = unique(genPLotList_Mouse_chr17$external_gene_name) )
genPLotList_Mouse_chr17 <- list(genPLotList_Mouse_chr17=genPLotList_Mouse_chr17)

##########################
# 1. plot for mouse
##########################

legend_title <- "Biotype Category"
myTitlePlot <- "Homologs to human chr21 in mouse chr16"
myTitlePlot2 <- "Homologs to human chr21 in mouse chr10"
myTitlePlot3 <- "Homologs to human chr21 in mouse chr17"

plotFC_bio1 <- ggplot(genPLotList_Mouse_chr16[[1]], aes(x=genPLotList_Mouse_chr16[[1]][,1], y=value,color =groups_sumupList_Mouse_chr16[[1]])) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[1]][,6]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[1]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio1 <-plotFC_bio1 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Mouse","Human")))  
#2
plotFC_bio2 <- ggplot(genPLotList_Mouse_chr16[[2]], aes(x=genPLotList_Mouse_chr16[[2]][,1], y=value,color =groups_sumupList_Mouse_chr16[[2]])) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[2]][,6]), size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[2]][,4]),size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio2 <-plotFC_bio2 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Mouse","Human")))+ theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") 
#3
plotFC_bio3 <- ggplot(genPLotList_Mouse_chr16[[3]], aes(x=genPLotList_Mouse_chr16[[3]][,1], y=value,color =groups_sumupList_Mouse_chr16[[3]])) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[3]][,6]), size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[3]][,4]),size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[3]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 

plotFC_bio3 <-plotFC_bio3 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Mouse","Human")))  + theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")
#4 
plotFC_bio4 <- ggplot(genPLotList_Mouse_chr16[[4]], aes(x=genPLotList_Mouse_chr16[[4]][,1], y=value,color =groups_sumupList_Mouse_chr16[[4]])) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[4]][,6]), size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Mouse_chr16[[4]][,4]),size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr16[[4]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio4 <-plotFC_bio4 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Mouse","Human")))  + theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none")

#chr10 mouse
plotFC_bio5 <- ggplot(genPLotList_Mouse_chr10[[1]], aes(x=genPLotList_Mouse_chr10[[1]][,1], y=value,color =groups_sumupList_Mouse_chr10[[1]])) +
	geom_point(aes(y=genPLotList_Mouse_chr10[[1]][,6]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot2) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr10[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Mouse_chr10[[1]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot2) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr10[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio5 <-plotFC_bio5 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Mouse","Human")))  
#chr17
plotFC_bio6 <- ggplot(genPLotList_Mouse_chr17[[1]], aes(x=genPLotList_Mouse_chr17[[1]][,1], y=value,color =groups_sumupList_Mouse_chr17[[1]])) +
	geom_point(aes(y=genPLotList_Mouse_chr17[[1]][,6]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot3) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr17[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Mouse_chr17[[1]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot3) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Mouse_chr17[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio6 <-plotFC_bio6 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Mouse","Human"))) + geom_vline(xintercept = 7,linetype = 2)+ geom_vline(xintercept = 24,linetype = 2)+ geom_vline(xintercept = 75,linetype = 2)
 


##########################
# 1. plot for Rat
##########################

legend_title <- "Biotype Category"
myTitlePlot <- "Homologs to human chr21 in rat chr11"
myTitlePlot2 <- "Homologs to human chr21 in rat chr20"

plotFC_bio_rat1 <- ggplot(genPLotList_Rat_chr11[[1]], aes(x=genPLotList_Rat_chr11[[1]][,1], y=value,color =groups_sumupList_Rat_chr11[[1]])) +
	geom_point(aes(y=genPLotList_Rat_chr11[[1]][,6]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Rat_chr11[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Rat_chr11[[1]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Rat_chr11[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio_rat1 <-plotFC_bio_rat1 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Rat","Human")))  
#2
plotFC_bio_rat2 <- ggplot(genPLotList_Rat_chr11[[2]], aes(x=genPLotList_Rat_chr11[[2]][,1], y=value,color =groups_sumupList_Rat_chr11[[2]])) +
	geom_point(aes(y=genPLotList_Rat_chr11[[2]][,6]), size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Rat_chr11[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Rat_chr11[[2]][,4]),size = 1.9,shape = 15) + 
xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Rat_chr11[[2]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio_rat2 <-plotFC_bio_rat2 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Rat","Human")))+ theme(plot.title = element_text(colour="white"),axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") 

#chr20
plotFC_bio_rat3 <- ggplot(genPLotList_Rat_chr20[[1]], aes(x=genPLotList_Rat_chr20[[1]][,1], y=value,color =groups_sumupList_Rat_chr20[[1]])) +
	geom_point(aes(y=genPLotList_Rat_chr20[[1]][,6]), size = 1.9,shape = 15) + 
labs( title= myTitlePlot2) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Rat_chr20[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) +
	geom_point(aes(y=genPLotList_Rat_chr20[[1]][,4]),size = 1.9,shape = 15) + 
labs( title= myTitlePlot2) + xlab("") + ylab("") + theme_gray(base_size = 14)+
scale_color_manual(legend_title,labels=levels(catDf$biotype_category),values=levels(groups_sumupList_Rat_chr20[[1]])) +
ggpubr::rotate_x_text()+ theme(plot.margin = unit(c(1, 0.7, 0, 0), "lines"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), legend.title=element_text(size=8,face="bold"), legend.text=element_text(size=7), axis.text=element_text(size=4), axis.text.x = element_text(angle=60,vjust=1,size = 6, face = "bold")) + guides(color = guide_legend(override.aes = list(size=5))) 
plotFC_bio_rat3 <-plotFC_bio_rat3 + theme(axis.text.y = element_text(size = 7, face = "bold")) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 1),limits = c(0, 2) ,breaks=c(0.5,2),labels = rev(c("Rat","Human"))) + geom_vline(xintercept = 99,linetype = 2)


#save plot
pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologschr21_RatSyntenicRegionScatterplot_Book.pdf"), width=11.69, height=9.27)
ggarrange(plotFC_bio_rat1, plotFC_bio_rat2,plotFC_bio_rat3,
	ncol = 1, nrow = 4,  align = "hv" )
plot(legend)
dev.off()

pdf(file=paste0(wd,f_results, f_homologs,f_figures,"homologschr21_MouseRatSyntenicRegionScatterplot_Book.pdf"), width=11.69, height=9.27)
ggarrange(plotFC_bio1, plotFC_bio2,plotFC_bio3, plotFC_bio4,
	ncol = 1, nrow = 4,  align = "hv" )
ggarrange(plotFC_bio5, plotFC_bio6,legend,
	ncol = 1, nrow = 4,  align = "hv" )
ggarrange(plotFC_bio_rat1, plotFC_bio_rat2,plotFC_bio_rat3,
	ncol = 1, nrow = 4,  align = "hv" )
#plot(legend)
dev.off()

df<-genesDf_specie_Mouse
dfRemote<- remoteHomologs_Mouse
dfPred<- predGenes_Mouse
speciesName<- c("Mouse","Rat", "Chimpanzee","Orangutan","Gorilla","Zebrafish", "C.elegans","Drosophila")
dfHigh<- list(Mouse=homologsHighConf_Mouse,Rat=homologsHighConf_Rat, Chimpanzee=homologsHighConf_Chimpanzee,Orangutan=homologsHighConf_Orangutan,Gorilla=homologsHighConf_Gorilla,Zebrafish=homologsHighConf_Zebrafish, Celegans= homologsHighConf_C.elegands,Drosophila= homologsHighConf_Drosophila)
dfPred<- list(Mouse=predGenes_Mouse,Rat=predGenes_Rat, Chimpanzee=predGenes_Chimpanzee,Orangutan=predGenes_Orangutan,Gorilla=predGenes_Gorilla,Zebrafish=predGenes_Zebrafish, Celegans= predGenes_C.elegands,Drosophila= predGenes_Drosophila)
dfRemote<- list(Mouse=remoteHomologs_Mouse,Rat=remoteHomologs_Rat, Chimpanzee=remoteHomologs_Chimpanzee,Orangutan=remoteHomologs_Orangutan,Gorilla=remoteHomologs_Gorilla,Zebrafish=remoteHomologs_Zebrafish, Celegans= remoteHomologs_C.elegands,Drosophila= remoteHomologs_Drosophila)

syntenic_regionDefinition_4Plot_function<- function(dfHigh,dfPred,dfRemote,speciesName,genesDf_ref,date){
	homologyResList<- list();
	homologyResFullList <- list();
	for (i in 1:length(speciesName)){
		specieName<- speciesName[[i]];
		df <- rbind(dfRemote[[i]],dfHigh[[i]],dfPred[[i]]);
		df1<- df[,c(1,2,5,6,3,8,9,10,21,13)];
		dfSel <- left_join(df1,genesDf_ref,by=c("ensembl_gene_id","external_gene_name"));
		colnames(dfSel)[c(11:15)] <- paste0("ref_", colnames(dfSel)[c(11:15)]);
		colnames(dfSel) <- gsub("\\.y","",gsub("\\.x","",colnames(dfSel)));
		dfSel <- dfSel[order(dfSel$chromosome),];
		dfSel <- dfSel[order(dfSel$ref_start_position),];
		homologyRes_fullTable<- dfSel;
		assign(paste0("homologyRes_fullTable",specieName),homologyRes_fullTable,.GlobalEnv);
		homologyResFullList[[i]] <- homologyRes_fullTable;
		names(homologyResFullList)[i] <- specieName;



		#only with high confidence homologs

		dfHighSel <- dfHigh[[i]][,c(1,2,5,6,3,8,9,10,21,13)];
		dfHighSel <- left_join(dfHighSel,genesDf_ref,by=c("ensembl_gene_id","external_gene_name"));
		colnames(dfHighSel)[c(11:15)] <- paste0("ref_", colnames(dfHighSel)[c(11:15)]);  
		colnames(dfHighSel) <- gsub("\\.y","",gsub("\\.x","",colnames(dfHighSel)));
		dfHighSel <- dfHighSel[order(dfHighSel$chromosome),];
		dfHighSel <- dfHighSel[order(dfHighSel$ref_start_position),];
		homologyResTable<- dfHighSel;
		assign(paste0("homologyResTable_",specieName),homologyResTable,.GlobalEnv);
		homologyResList[[i]] <- homologyResTable;
		names(homologyResList)[i] <- specieName;
		i=i+1
	};
	write.xlsx(homologyResList, file = paste0(wd,f_results,f_homologs,f_tables, "syntenicConservation_to_human_chr21_ensembl_HighHomologsOnly_",date, ".xlsx"));
	write.xlsx(homologyResFullList, file = paste0(wd,f_results,f_homologs,f_tables, "syntenicConservation_to_human_chr21_ensembl_allHomologs_",date, ".xlsx"));

};
syntenic_regionDefinition_4Plot_function(dfHigh,dfPred,dfRemote,speciesName,genesDf_ref,"020519")
		dfSel <- homologyResTable_Mouse[order(homologyResTable_Mouse$ref_start_position),];

#gnes_homolos_Mouse <-gnes_homolos_Mouse[gtools::mixedorder(as.character(gnes_homolos_Mouse$chromosome_name)),];	

#####################################################################################################
############################### END. Finished on the  0319 ########################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################





sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.5 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] openxlsx_4.1.0       gplots_3.0.1.1       ggplot2_3.1.0       
[4] gtable_0.2.0         RColorBrewer_1.1-2   dplyr_0.7.6         
[7] tidyr_0.8.1          biomaRt_2.34.2       BiocInstaller_1.28.0

loaded via a namespace (and not attached):
 [1] progress_1.2.0       gtools_3.8.1         tidyselect_0.2.5    
 [4] purrr_0.2.5          colorspace_1.4-0     stats4_3.4.4        
 [7] utf8_1.1.4           blob_1.1.1           XML_3.98-1.16       
[10] rlang_0.2.2          pillar_1.3.0         glue_1.3.0          
[13] withr_2.1.2          DBI_1.0.0            BiocGenerics_0.24.0 
[16] bit64_0.9-7          bindrcpp_0.2.2       bindr_0.1.1         
[19] plyr_1.8.4           stringr_1.3.1        munsell_0.5.0       
[22] zip_1.0.0            caTools_1.17.1.1     memoise_1.1.0       
[25] Biobase_2.38.0       IRanges_2.12.0       curl_3.2            
[28] parallel_3.4.4       fansi_0.4.0          AnnotationDbi_1.40.0
[31] Rcpp_1.0.0           KernSmooth_2.23-15   scales_1.0.0        
[34] gdata_2.18.0         S4Vectors_0.16.0     bit_1.1-14          
[37] hms_0.4.2            digest_0.6.17        stringi_1.2.4       
[40] grid_3.4.4           cli_1.0.1            tools_3.4.4         
[43] bitops_1.0-6         magrittr_1.5         RCurl_1.95-4.11     
[46] lazyeval_0.2.1       tibble_1.4.2         RSQLite_2.1.1       
[49] crayon_1.3.4         pkgconfig_2.0.2      prettyunits_1.0.2   
[52] assertthat_0.2.0     httr_1.4.0           R6_2.3.0            
[55] compiler_3.4.4  

