########################################################################
########################################################################
#QC: PCA plot, Volcano Plot
########################################################################
########################################################################

########################################################################
###################  PACKAGES NEEDED  ##################################
source("http://bioconductor.org/biocLite.R")

library("tidyr")
library("dplyr")
library("RColorBrewer")
library("dplyr")
library("VennDiagram")
library("gplots")
library("tidyr")

########################################################################
###################  PACKAGES NEEDED  ##################################


#Before begining:
set.seed(22); # need to set it to get always the same random results and plots

#wd:
#####################################################
setwd("/Users/marmmoreno/Desktop/onWorking/");  #mac
wd <- getwd();

f_results <- "/results/"
f_Rdata <- "RData/"
f_venn <- "vennDiagram/"

dir.create(file.path(getwd (), f_results), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_venn), showWarnings = F)


########################################################################
############################# vennDiagrams #############################
########################################################################
########################################################################

########################################
#1. DEGs vennDiagrams
########################################

##### A. Loading DEGs from models
########################################
########################################################################
#NOTE IMPORTANT!:
#  If the previous function is not run, instead the files are importing
#  or other models files are imported, as they will have the same variables
#  names use this specific way of importing in R where name is the model 
#  name to distinguish the results from each model.
########################################################################

nameArraysDEGsToImport <-c("Ts65Dn","Dp5Dp1","Dp1Yey","TgDyrk1a","Dp5Yah","Dp1Rhr","Dp3Yah")

f_InputarraysDEGs <- "/RData/arraysFcrosDEGs/regByModel/regFormated/"

namefolderRNASeqInfoImport <-c("TgDyrk1a","G8Dlx","G8Het","G7Dlx","adultTgDyrk1a",
	"RNASeq_Dp1Yey","Ts66","Ts66","rat_Rno20delDup_Del_Wt","rat_Rno20delDup_Dup_Wt",
	"rat_Rno20delDup_DelDup_Wt","rat_DupCBS","rat_Dyrk1a_cko")

nameDEGsRNASeqToImport <-c("TgDyrk1a","G8Dlx","G8Het","G7Dlx","adultTgDyrk1a","RNASeq_Dp1Yey",
	"hippo","ento","rat_Rno20delDup_Del_Wt","rat_Rno20delDup_Dup_Wt","rat_Rno20delDup_DelDup_Wt",
	"rat_DupCBS","rat_Dyrk1a_cko")

import_DEGs_for_venn.function <- function(nameArraysDEGsToImport,nameDEGsRNASeqToImport,f_InputarraysDEGs){
      ###############################################################
      # NOTE: this function will help to import all the DEGs data coming 
      # from different experiment (RNASeq/microarrays) if all the models
      # per technology are stored in the same directories (one for rnaseq, 
      # another for the arrays)
      # output: 2 list files will be created, 
      #  1. DEGsListReg: list  containing the DEGs name plus a
      #     second column with the regulation
      #  2. DEGsList: DEGs input ready to be feed to venn function
      ###############################################################
      DEGsListReg <-list()
      DEGsList <-list()
      for (i in 1:length(nameArraysDEGsToImport)){
            DEGs <- read.table(file=paste0(wd,f_InputarraysDEGs,nameArraysDEGsToImport[i],"_DEGs_regSense.txt"), header=T, sep = "\t")
            DEGs$model <- paste0(nameArraysDEGsToImport[i],"_Array")
            DEGsListReg[[i]]<- DEGs
            names(DEGsListReg)[i] <- paste0(nameArraysDEGsToImport[i],"_Array")
            DEGs.tmp <- as.data.frame(as.vector(unique(DEGs[,1]))) #comment="",quote=NULL
            colnames(DEGs.tmp) <- paste0(nameArraysDEGsToImport[i],"_Array")
            DEGsList[[i]]<- DEGs.tmp
            names(DEGsList)[i] <- paste0(nameArraysDEGsToImport[i],"_Array")
            i=i+1
      };

      n= length(DEGsList) +1
      for (lol in 1:length(nameDEGsRNASeqToImport)){
            
            if (nameDEGsRNASeqToImport[lol]=="TgDyrk1a"){
                  nameNew <- "TgDyrk1a_E15.5"
            } else if (nameDEGsRNASeqToImport[lol]=="G8Dlx"){
                  nameNew <- "G8Dlx_E15.5"
            } else if (nameDEGsRNASeqToImport[lol]=="G8Het"){
                  nameNew <- "Dyrk1a_het_E15.5_seRNASeq"
            } else if (nameDEGsRNASeqToImport[lol]=="G7Dlx"){
                  nameNew <- "G7Dlx_E15.5"
            } else{
                  nameNew <- nameDEGsRNASeqToImport[lol]
            }
            DEGs <- local({load(paste0(wd,"/results/gage/",namefolderRNASeqInfoImport[lol],"/DFA/RData/",nameDEGsRNASeqToImport[lol],"_DEGsALLInfo.RData")) 
                  stopifnot(length(ls())==1) 
            environment()[[ls()]]
                  })
            DEGs <- unique(DEGs[,c(2,6)])
            DEGs$model <- paste0(nameNew,"_seRNASeq")
            DEGsListReg[[n]]<- DEGs
            names(DEGsListReg)[n] <- paste0(nameNew,"_seRNASeq")
            DEGs.tmp <- as.data.frame(as.vector(unique(DEGs[,1]))) #comment="",quote=NULL
            colnames(DEGs.tmp) <- paste0(nameNew,"_seRNASeq")
            DEGsList[[n]]<- DEGs.tmp
            names(DEGsList)[n] <- paste0(nameNew,"_seRNASeq")
            n=n+1
            lol=lol+1
      }
      assign("DEGsList",DEGsList,.GlobalEnv)
      assign("DEGsListReg",DEGsListReg,.GlobalEnv)

};
import_DEGs_for_venn.function(nameArraysDEGsToImport,nameDEGsRNASeqToImport,f_InputarraysDEGs)


modelsDEGsNb<- t(as.data.frame(lapply(DEGsList, function(x) nrow(x))))
colnames(modelsDEGsNb)<- "modelsDEGsNb"

# names(DEGsList)
# [1] "Dp1Rhr_Array"             "Dp1Yey_Array"            
# [3] "Dp3Yah_Array"             "Dp5Dp1_Array"            
# [5] "Dp5Yah_Array"             "TgDyrk1a_Array"          
# [7] "Ts65Dn_Array"             "TgDyrk1a_E15.5_seRNASeq" 
# [9] "G8Dlx_E15.5_seRNASeq"     "Dyrk1a_E15.5_seRNASeq"
#[11] "G7Dlx_E15.5_seRNASeq"     "adultTgDyrk1a_seRNASeq"  
#[13] "Dp1Yey_seRNASeq"          "Ts66.hippo_seRNASeq"     
#[15] "Ts66.ento_seRNASeq" 

#annote with the homologs genes the rats to mouse
#####################################################
homologs_Annotation.function <- function(specieReferenceBiomaRt,specieReferenceName,inputGenes,speciesQueryBiomaRt,speciesQueryNames,nameDataDEGs,DEGsList,dateAnalysis){
	##################################################################################
	# Function to annotate all the genes of the chr of the reference specie with
	# its homologs in the species of interest defined as query
	# Its produce the stats. Input files are:
	#
	##################################################################################	
	f_biocatInput <- "/Users/marmmoreno/Documents/YH_Lab/collabs_outside/DS_models_book/results/RData/"
	load(file=paste0(f_biocatInput, "/biotypeCategory_pertenenceData.RData") ) #biotypeCat
	
	inputGenes <- unique(gsub("_.*","",inputGenes))
	#2. Quering ENSEMBL 
	######################
	#ensembl <- useEnsembl(biomart = "ensembl")
	#listDatasets(mart = ensembl)
	#for the reference specie:
	ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = specieReferenceBiomaRt,mirror="useast") 
	attributes = listAttributes(ensembl);
	attributesHom = attributes[which(attributes[,3]=="homologs"),];	
	rownames(attributesHom) <- 1:nrow(attributesHom);
	attribRef <-c(as.numeric(rownames(attributes[which(attributes[,1]=='ensembl_gene_id'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='external_gene_name'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='chromosome_name'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='start_position'),][1,])),as.numeric(rownames(attributes[which(attributes[,1]=='end_position'),][1,])),
		as.numeric(rownames(attributes[which(attributes[,1]=='gene_biotype'),][1,]))  );
	attribRef <-attributes[attribRef,1]; 


	genesDf_ref = getBM(attributes=attribRef,filters="external_gene_name",values=inputGenes, mart=ensembl); 
	values<-unique(genesDf_ref[,1]); 
	assign(paste0("genesDf_ref", "_", specieReferenceBiomaRt),genesDf_ref,.GlobalEnv);

	perchrHomologsFList <- list();
	for (i in 1:length(speciesQueryBiomaRt)){	
		f_specie <- paste0(gsub(" ","_",speciesQueryNames[i]),"/");
		name <- speciesQueryNames[i];
		specieAnalysis <- speciesQueryBiomaRt[i];

		#1. Creating the folders
		######################
		f_homologs<- "homologs/"
		f_tables <- "Tables/"
		f_figures <- "Plots/"
		f_stats <- "stats/"
		dir.create(file.path(getwd (), f_results,f_homologs), showWarnings = F);
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
		genesDf_specie<- genesDf_specie[,c(grep("ensembl_gene_id" ,colnames(genesDf_specie)),
			grep("external_gene_name" ,colnames(genesDf_specie)),grep("biotype_category" ,colnames(genesDf_specie)),
			grep("gene_biotype" ,colnames(genesDf_specie)),grep("homolo" ,colnames(genesDf_specie))  )];
		colnames(genesDf_specie) <- gsub(".*_homolog_","",colnames(genesDf_specie));
		genesDf_specie <-genesDf_specie[-which(genesDf_specie$ensembl_gene==""),];# not homolog found #all homologs		
		remoteHomologs <- genesDf_specie[which(genesDf_specie[,20]==0 & genesDf_specie[,14] >15),]; #remote, confidence =0, seq similarity > 15% i the query animal
		remoteHomologs$HomologyCategory <- "remote homologs";
		homologsHighConf<- genesDf_specie[which(genesDf_specie[,20]==1),]; # high confidence: conf score=1
		homologsHighConf$HomologyCategory <- "High confidence homologs";
		predGenes <- genesDf_specie[is.na(genesDf_specie[,20]),]; #creating this category because, the conf level of these genes is 0, but the % seq shared is high, 
		if (dim(predGenes)[1]==0) {
			genesDf_specieF<- rbind(homologsHighConf,remoteHomologs);
		} else {
			predGenes$HomologyCategory <- "predicted genes";
			genesDf_specieF<- rbind(homologsHighConf,remoteHomologs,predGenes);
		
		}
		
		#and it may be cause most of these genes are predicted in the specie of reference and we dont have 
		#enough proof of them yet that is why the conf score is 0 in the end
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
		write.xlsx(resultsHomologs, file = paste0(wd,f_results,f_homologs, f_specie, f_tables, specieReferenceName, "_specieRef_homologs_with_query_specie_",speciesQueryNames[i],"_",dateAnalysis,".xlsx"));
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

		#Name of the genes with homolog names
		#####################################
		homologs<- unique(df$HomologyCategory);
		genesWithHomologsList <-  list();
		nb <-length(DEGsList);
		for (lol in 1:length(homologs)){
			genesWithHomologsList[[lol]] <- as.data.frame(unique(df[which(df$HomologyCategory==unique(df$HomologyCategory)[lol]),'external_gene_name']));
			names(genesWithHomologsList)[[lol]] <- paste0(speciesQueryBiomaRt[i],"_homologs_of_",nameDataDEGs,"_",homologs[lol]);
			DEGsList[[nb+lol]] <- as.data.frame(unique(df[which(df$HomologyCategory==unique(df$HomologyCategory)[lol]),'associated_gene_name']));
			names(DEGsList)[[nb+lol]]<- gsub(" ","_",paste0(speciesQueryBiomaRt[i],"_homologs_of_",nameDataDEGs,"_",homologs[lol]));	
			lol=lol+1;
		};
		nb <-length(DEGsList);
		DEGsList[[nb+1]] <- as.data.frame(unique(df[,'associated_gene_name']));
		names(DEGsList)[[nb+1]]<- paste0(speciesQueryBiomaRt[i],"_homologs_of_",nameDataDEGs,"_AllHomologs");	

		

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
			HomologyCategory.tmp <-as.data.frame(unique(df[,c(2,21)]) %>% group_by(HomologyCategory) %>% count()); #nb of unique  genes with at least one homolog. 
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
			HomologyCategory$`Nb DEGs` <- length(inputGenes)
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
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[order(perchrHomologsNoOrtCat.tmp[,5],decreasing = TRUE ),];
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[order(match(perchrHomologsNoOrtCat.tmp[,2],unique(biotypeCat[,2]) )),];				
			perchrHomologsNoOrtCat.tmp <- perchrHomologsNoOrtCat.tmp[order(match(perchrHomologsNoOrtCat.tmp[,1],homology_typeOrder )),];	
			perchrHomologsNoOrtCat <- unique(full_join(perchrHomologsNoOrtCat,perchrHomologsNoOrtCat.tmp, by=c("HomologyCategory", "biotype_category","chromosome")));		
		}
		 
		i=i+1
	};
	#Saving the results
	resListStats <- list(HomologyCat=HomologyCategory,HomologyCatBio=HomologyCategoryBio,HomCatBio.queryHom=HomologyCategoryBio.queryHomologs,perchrHom=perchrHomologsNoOrtCat, perchrHomologsOrthoType=perchrHomologs)
	write.xlsx(resListStats, file = paste0(wd,f_results,f_homologs, f_stats, specieReferenceBiomaRt[i],"_homology_per_",nameDataDEGs,"_Stats","_",dateAnalysis,".xlsx"));
	assign(paste0("Homologs_Stats_",nameDataDEGs),resListStats,.GlobalEnv);
	assign(paste0("perchrHomologsFList",nameDataDEGs),perchrHomologsFList,.GlobalEnv);
	assign("DEGsListF",DEGsList,.GlobalEnv);
	assign(paste0("genesWithHomologsList_",nameDataDEGs),genesWithHomologsList,.GlobalEnv);
};
specieReferenceBiomaRt<- "rnorvegicus_gene_ensembl"
specieReferenceName <- "Rat"
speciesQueryNames <- c("Mouse","Human");
homologs_Annotation.function(specieReferenceBiomaRt,specieReferenceName,inputGenes,speciesQueryBiomaRt,speciesQueryNames,nameDataDEGs,DEGsList,"040719");

########################################
#
########################################
model_colors <- data.frame(matrix(nrow=length(names(DEGsList)), ncol=4))
colnames(model_colors) <- c("model","fill","cat.col","labels")
model_colors$model <- names(DEGsList)
model_colors$fill <- ifelse(model_colors$model=="Dp1Rhr_Array","darksalmon",
	ifelse(model_colors$model=="Dp1Yey_Array","darkorchid2",
	ifelse(model_colors$model=="Dp3Yah_Array","azure3",
	ifelse(model_colors$model=="Dp5Dp1_Array","deeppink",
	ifelse(model_colors$model=="Dp5Yah_Array","blue2",
	ifelse(model_colors$model=="TgDyrk1a_Array","lightgoldenrod2",
	ifelse(model_colors$model=="Ts65Dn_Array","aquamarine2",
	ifelse(model_colors$model=="TgDyrk1a_E15.5_seRNASeq","pink",
	ifelse(model_colors$model=="G8Dlx_E15.5_seRNASeq","palevioletred1",
	ifelse(model_colors$model=="Dyrk1a_het_E15.5_seRNASeq","lightsteelblue1",
	ifelse(model_colors$model=="G7Dlx_E15.5_seRNASeq","rosybrown1",
	ifelse(model_colors$model=="adultTgDyrk1a_seRNASeq","peachpuff4",
	ifelse(model_colors$model=="RNASeq_Dp1Yey_seRNASeq","thistle1",
	ifelse(model_colors$model=="Ts66.hippo_seRNASeq","indianred1",
	ifelse(model_colors$model=="Ts66.ento_seRNASeq","navajowhite2",
	ifelse(model_colors$model=="rat_Rno20delDup_Del_Wt","lightpink3",
	ifelse(model_colors$model=="rat_Rno20delDup_Dup_Wt","darkolivegreen2",
	ifelse(model_colors$model=="rat_Rno20delDup_DelDup_Wt","cornflowerblue",
	ifelse(model_colors$model=="rat_DupCBS","brown3",
	ifelse(model_colors$model=="rat_Dyrk1a_cko","khaki3",
	ifelse(model_colors$model=="rat_Rno11Rno20_DoubleDup_Wt","springgreen",
	ifelse(model_colors$model=="rat_Rno11Rno20_Rno11Dup_Wt","snow3",
	ifelse(model_colors$model=="rat_Rno11Rno20_Rno20Dup_Wt","darkolivegreen3",
		NA)))))))))))))))))))))))


model_colors$cat.col <- ifelse(model_colors$model=="Dp1Rhr_Array","darksalmon",
	ifelse(model_colors$model=="Dp1Yey_Array","darkorchid3",
	ifelse(model_colors$model=="Dp3Yah_Array","grey28",
	ifelse(model_colors$model=="Dp5Dp1_Array","deeppink",
	ifelse(model_colors$model=="Dp5Yah_Array","blue2",
	ifelse(model_colors$model=="TgDyrk1a_Array","goldenrod4",
	ifelse(model_colors$model=="Ts65Dn_Array","aquamarine4",
	ifelse(model_colors$model=="TgDyrk1a_E15.5_seRNASeq","pink4",
	ifelse(model_colors$model=="G8Dlx_E15.5_seRNASeq","palevioletred4",
	ifelse(model_colors$model=="Dyrk1a_het_E15.5_seRNASeq","lightsteelblue4",
	ifelse(model_colors$model=="G7Dlx_E15.5_seRNASeq","rosybrown3",
	ifelse(model_colors$model=="adultTgDyrk1a_seRNASeq","peachpuff4",
	ifelse(model_colors$model=="RNASeq_Dp1Yey_seRNASeq","thistle4",
	ifelse(model_colors$model=="Ts66.hippo_seRNASeq","indianred4",
	ifelse(model_colors$model=="Ts66.ento_seRNASeq","navajowhite4",
	ifelse(model_colors$model=="rat_Rno20delDup_Del_Wt","lightpink4",
	ifelse(model_colors$model=="rat_Rno20delDup_Dup_Wt","darkolivegreen4",
	ifelse(model_colors$model=="rat_Rno20delDup_DelDup_Wt","cornflowerblue",
	ifelse(model_colors$model=="rat_DupCBS","brown4",
	ifelse(model_colors$model=="rat_Dyrk1a_cko","khaki4",
	ifelse(model_colors$model=="rat_Rno11Rno20_DoubleDup_Wt","springgreen4",
	ifelse(model_colors$model=="rat_Rno11Rno20_Rno11Dup_Wt","snow4",
	ifelse(model_colors$model=="rat_Rno11Rno20_Rno20Dup_Wt","darkolivegreen4",
		NA)))))))))))))))))))))))

model_colors$labels <- gsub("6\\.","6 ",gsub("TgDyrk1a"," Tg(Dyrk1a)",gsub("_"," ",names(DEGsList))))

########################################

########## Folders #####################
library("VennDiagram")
f_venn <- "vennDiagram/"
f_vdata <- "vennTables/"
f_vplots <- "vennPlots/"
########################################


vennD2.function<- function(model_colors,model1,model2,label1,label2,DEGsList,DEGsListReg){
	library("VennDiagram");
	category <-c(label1,label2)
	f_category <- paste0(label1,"_",label2,"/")
	name_category <- paste0(label1,"_",label2)
	#creating the folders needed
	dir.create(file.path(getwd (), f_results, f_venn), showWarnings = F)
	dir.create(file.path(getwd (), f_results, f_venn,f_category), showWarnings = F)
	dir.create(file.path(getwd (), f_results, f_venn,f_category,f_vdata), showWarnings = F)
	dir.create(file.path(getwd (), f_results, f_venn,f_category,f_vplots), showWarnings = F)

	#overlaps
	file1 <- DEGsList[[grep(model1,names(DEGsList))]]
	file2 <- DEGsList[[grep(model2,names(DEGsList))]]
	area1 <- length(unique(unlist(file1)))
	area2 <- length(unique(unlist(file2)))
	n12 <- length(Reduce(intersect, list(t(file1),t(file2))))
	
	#saving data
	data <- data.frame(area1,area2,n12)
	colnames(data) <- c("area1", "area2", "Intercept_12")
	write.table(data, file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_sumup.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
	#saving each intercept list of DEGs
	write.table(data.frame(area_1=unique(unlist(file1))), file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_data_intercept_list_area1.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
	write.table(data.frame(area_2=unique(unlist(file2))), file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_data_intercept_list_area2.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
	write.table(data.frame(n12=Reduce(intersect, list(t(file1),t(file2)))), file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_data_intercept_list_n12.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
	
	#ploting the venn diagram
	fill<- 	c(model_colors[grep(model1,model_colors$model),2],model_colors[grep(model2,model_colors$model),2])
	cat.col <-c(model_colors[grep(model1,model_colors$model),3],model_colors[grep(model2,model_colors$model),3])
	
	pdf(file= paste0(wd,f_results, f_venn,f_category,f_vplots, name_category,".pdf"),height = 2.5, width = 2.5) 
	venn.plot1 <- draw.pairwise.venn(
	  area1 = area1,
	  area2 = area2,
	  cross.area = n12,
	  category = category,
	  fill = fill,
	  cat.col = cat.col,
	  margin = 0.3,
	  cex = 0.7,
	  alpha = 0.7,
	  cat.cex = 0.6,
	  cat.default.pos = "outer",
	  cat.just=list(c(0.5,-4.4), c(0.5,-4.4)), #1.hippo, 2 ento. c(a,b)  a:horizontl, b vertical. -b is up,-a drcha
	  ext.pos = 2.9,
	  ext.dist = -0.05,
	  ext.length = 0.85,
	  euler.d=FALSE,
	  scaled=FALSE,
	  ind=TRUE
	);
	assign(paste0("venn.plot1"),venn.plot1,.GlobalEnv)
	assign(paste0("n12",model1,model2),data.frame(n12=Reduce(intersect, list(t(file1),t(file2)))),.GlobalEnv)
	grid.draw(venn.plot1);
	grid.draw(venn.plot1)
	dev.off()
	grid.draw(venn.plot1)
	grid.draw(venn.plot1)
	dev.off()
	save(venn.plot1,file=paste0(wd,f_results, f_venn,f_category,f_vplots, name_category,".Rdata")) 
	
	#stats shared region of DEGs    
	n12Genes<- data.frame(n12=Reduce(intersect, list(t(file1),t(file2))))
	cond1 <- DEGsListReg[[model1]][,c(1,2)]
	cond1 <-as.data.frame(cond1[which(cond1[,1] %in% n12Genes[,1]),])
	cond1.Up <- cond1[which(cond1[,2]=="UpRegulated"),1]  
	cond1.Down <- cond1[which(cond1[,2]=="DownRegulated"),1]  
	
	cond2 <- DEGsListReg[[model2]][,c(1,2)]
	cond2 <-as.data.frame(cond2[which(cond2[,1] %in% n12Genes[,1]),])
	cond2.Up <- cond2[which(cond2[,2]=="UpRegulated"),1]  
	cond2.Down <- cond2[which(cond2[,2]=="DownRegulated"),1]  
	

	#total number hits
	############################################
	total.Shared <- length(unique(cond1[,1][cond1[,1] %in% cond2[,1]]))
	# shared, up in both
	############################################
	lshared.Up_Up <-length(unique(cond1.Up[cond1.Up %in% cond2.Up]))
	shared.Up_Up <-unique(cond1.Up[cond1.Up %in% cond2.Up])

	#shared, down in one cond, up in the other
	############################################
	lshared.1Up_2Down<- length(unique(cond1.Up[cond1.Up %in% cond2.Down]))
	shared.1Up_2Down<- unique(cond1.Up[cond1.Up %in% cond2.Down])
	
	lshared.2Up_1Down<- length(unique(cond2.Up[cond2.Up %in% cond1.Down]))
	shared.2Up_1Down<- unique(cond2.Up[cond2.Up %in% cond1.Down])

	# shared, down in both
	############################################
	lshared.Down_Down<-length(unique(cond1.Down[cond1.Down %in% cond2.Down])) 
	shared.Down_Down<-unique(cond1.Down[cond1.Down %in% cond2.Down]) 
	
	assign(paste0("lshared.Up_Up",model1,model2),lshared.Up_Up,.GlobalEnv)
	assign(paste0("shared.Up_Up",model1,model2),shared.Up_Up,.GlobalEnv)
	
	assign(paste0("lshared.1Up_2Down",model1,model2),lshared.1Up_2Down,.GlobalEnv)
	assign(paste0("shared.1Up_2Down",model1,model2),shared.1Up_2Down,.GlobalEnv)
	
	assign(paste0("lshared.2Up_1Down",model1,model2),lshared.2Up_1Down,.GlobalEnv)
	assign(paste0("shared.2Up_1Down",model1,model2),shared.2Up_1Down,.GlobalEnv)
	
	assign(paste0("lshared.Down_Down",model1,model2),lshared.Down_Down,.GlobalEnv)
	assign(paste0("shared.Down_Down",model1,model2),shared.Down_Down,.GlobalEnv)
	assign(paste0("f_category"),f_category,.GlobalEnv)
	assign(paste0("name_category"),name_category,.GlobalEnv)

	#stats shared/unique nb of genes and identity of the genes
	cond1$model <- model1
	colnames(cond1)[2] <- paste0(colnames(cond1)[2],".",model1)
	
	cond2$model <- model2
	colnames(cond2)[2] <- paste0(colnames(cond2)[2],".",model2)
	a <- full_join(cond1,cond2, by=colnames(cond1)[1])
	a$fmodel <- paste0(a[,3],": ",a[,5])
	a <- a[,-c(3,5,7)]
	assign(paste0("sharedBetween2models.",model1,"_",model2),a,.GlobalEnv)
	sharedList <- list(stats=data.frame(ltotal.Shared=ltotal.Shared,lshared.Up_Up=lshared.Up_Up,
		lshared.1Up_2Down=lshared.1Up_2Down,lshared.2Up_1Down=lshared.2Up_1Down,
		lshared.Down_Down=lshared.Down_Down),sharedDEGs=a,total.Shared=total.Shared,
		shared.Up_Up=shared.Up_Up,shared.1Up_2Down=shared.1Up_2Down,shared.2Up_1Down=shared.2Up_1Down,
		shared.Down_Down=shared.Down_Down)
	assign(paste0("sharedList.",model1,model2),sharedList,.GlobalEnv)
	write.xlsx(sharedList, file = paste0(wd,f_results, f_venn,f_category,f_vdata, "dataVenn_sharedList.xlsx"))


};
vennD2.function(model_colors,"Ts66.hippo_seRNASeq","adultTgDyrk1a_seRNASeq",model_colors[14,'labels'],model_colors[12,'labels'],DEGsList,DEGsListReg);

#####################################################################################################
############################### END. Finished on the  230219 ########################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################

