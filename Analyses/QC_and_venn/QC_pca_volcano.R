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
library("DESeq2")
require("vsn")
require("hexbin")
library("scatterplot3d")
library("gridGraphics");
library("cowplot");
library("ggpubr");
library("ggplot2");

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
f_pca<- "pca/"
f_dist <- "SamplesDistance/"
f_data <- "vennTables/"
f_plots <- "vennPlots/"
f_volcano <- "volcanoPlots/"

dir.create(file.path(getwd (), f_results), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_dist), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_pca), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_pca,f_Rdata), showWarnings = F)
dir.create(file.path(getwd (), f_results, f_volcano), showWarnings = F)


######################
#A PCA EGs & DEGs
######################

#For RNASeq data analysed with DESeq2 we should use the normalised counts ( head(cnts.norm) ) to do 
#the pca, or use the normalised counts after applying the variance stabilizing transformation 
# which is roughly similar to putting the data on the log2 scale, while also dealing with
# the sampling variability of low counts (assay(vsd) ). Better to use the vsd.
# Plotting pca using all EGs or selected DEGs

#1.calculating vsd
deseq2_normalizationEffect.function <- function(countData, metafile){
	#calculate the normalization vsd and effect of normalization
	dds <- DESeqDataSetFromMatrix(countData = countData, colData = metafile, design = ~ condition)
	dds <- dds[ rowSums(counts(dds)) > 1, ]  #minimal pre-filtering to remove rows that have only 0 or 1 read
	vsd <- vst(dds, blind = FALSE)
	ntd <- normTransform(dds) # this gives log2(n + 1)
	#head(assay(vsd), 3)
	#colData(vsd)
	rld <- rlog(dds, blind = FALSE)
	#head(assay(rld), 3)
	dds <- estimateSizeFactors(dds)
	df <- bind_rows(
	  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
	         mutate(transformation = "log2(x + 1)"),
	  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
	  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
	colnames(df)[1:2] <- c("x", "y")  
	assign("df_normEffect",df,.GlobalEnv)
	assign("dds",dds,.GlobalEnv)
	assign("vsd",vsd,.GlobalEnv)
	assign("rld",rld,.GlobalEnv)
	normEffect <- ggplot(df_normEffect, aes(x = x, y = y)) + geom_hex(bins = 80) +
 	 	coord_fixed() + facet_grid( . ~ transformation)  + theme_gray(base_size = 14)
	pdf(paste0(wd,f_results, f_norm,f_deseq2,f_name, "_ML_Scatterplot effect of normalization Deseq2.pdf", sep=""))
	print(normEffect)
	meanSdPlot(assay(ntd))
	meanSdPlot(assay(vsd))
	meanSdPlot(assay(rld))
	dev.off()
};

deseq2_normalizationEffect.function(RawcountData,colData);

#2. calculating pca components

pca_inputData.function <- function(vsd,DEGs,name){
	##########################################################
	#1. DESeq2 plotPCA function:only possible to use 
	# in the deseq2 object vsd (cannot be subsetted 
	# to select DEGs). Only give us the first 2 components
	##########################################################
	pcaData.all <- plotPCA(vsd, intgroup =  "condition", returnData = TRUE,ntop = dim(assay(vsd))[1])
	percentVar.all <- round(100 * attr(pcaData.all, "percentVar"))
	pcaData.all$labels<-factor(gsub("wt","WT",gsub("_.*_"," ",pcaData.all$name)),
		levels=gsub("wt","WT",gsub("_.*_"," ",pcaData.all$name)))
	assign("percentVar.all",percentVar.all,.GlobalEnv)
	assign("pcaData.all",pcaData.all,.GlobalEnv)
	##########################################################
	#2. calculating the PC using prcomp 
	# give us the n number of components, we can plot them all 
	# and prepare a table
	##########################################################
	###########
	#EGs ######
	###########	
	pca.EGs <- prcomp(t(assay(vsd)), scale = TRUE) #center = TRUE, scale. = TRUE
	PCA.EGs <- pca.EGs$x
	# Cumulative Proportion 
	vars.EGs <- apply(pca.EGs$x, 2, var)  
	props.EGs <- vars.EGs / sum(vars.EGs)
	percentagesacc.EGs <- cumsum(props.EGs)*100
	#Proportion of Variance #########################
	eig.EGs <- (pca.EGs$sdev)^2
	variances.EGs <- eig.EGs*100/sum(eig.EGs) # Variances in percentage, 
	all_variance.EGs <- data.frame(eig= signif(eig.EGs, digits = 6), 
		variance = signif(variances.EGs, digits = 4), cumvariance = round(percentagesacc.EGs, digits = 2))
	PCA1.EGs <- PCA.EGs[,"PC1"] 
	PCA2.EGs <- PCA.EGs[,"PC2"]
	PCA3.EGs <- PCA.EGs[,"PC3"]
	condition.EGs <- vsd$condition 
	PC1A.EGs <- signif(all_variance.EGs$variance[1], digits=4)
	PC2A.EGs <- signif(all_variance.EGs$variance[2], digits=4)
	PC3A.EGs <- signif(all_variance.EGs$variance[3], digits=4)
	###########
	#DEGs #####
	###########	
	vsdDEGs <- assay(vsd)[which(rownames(assay(vsd)) %in% DEGs$EnsembleGeneName),]
	DEGsNb<-dim(vsdDEGs)[1]
	pca <- prcomp(t(vsdDEGs), scale = TRUE) #center = TRUE, scale. = TRUE
	PCA <- pca$x
	# Cumulative Proportion 
	vars <- apply(pca$x, 2, var)  
	props <- vars / sum(vars)
	percentagesacc <- cumsum(props)*100
	#Proportion of Variance #########################
	eig <- (pca$sdev)^2
	variances <- eig*100/sum(eig) # Variances in percentage, 
	all_variance.DEGs <- data.frame(eig= signif(eig, digits = 6), 
		variance = signif(variances, digits = 4), cumvariance = round(percentagesacc, digits = 2))
	PCA1.DEGs <- PCA[,"PC1"] 
	PCA2.DEGs <- PCA[,"PC2"]
	PCA3.DEGs <- PCA[,"PC3"]
	condition <- vsd$condition 
	PC1A.DEGs <- signif(all_variance.DEGs$variance[1], digits=4)
	PC2A.DEGs <- signif(all_variance.DEGs$variance[2], digits=4)
	PC3A.DEGs <- signif(all_variance.DEGs$variance[3], digits=4)
	assign("PCA1.DEGs",PCA1.DEGs,.GlobalEnv)
	assign("PCA2.DEGs",PCA2.DEGs,.GlobalEnv)
	assign("PCA3.DEGs",PCA3.DEGs,.GlobalEnv)
	assign("PCA1.EGs",PCA1.EGs,.GlobalEnv)
	assign("PCA2.EGs",PCA2.EGs,.GlobalEnv)
	assign("PCA3.EGs",PCA3.EGs,.GlobalEnv)
	assign("PC1A.DEGs",PC1A.DEGs,.GlobalEnv)
	assign("PC2A.DEGs",PC2A.DEGs,.GlobalEnv)
	assign("PC3A.DEGs",PC3A.DEGs,.GlobalEnv)
	assign("PC1A.EGs",PC1A.EGs,.GlobalEnv)
	assign("PC2A.EGs",PC2A.EGs,.GlobalEnv)
	assign("PC3A.EGs",PC3A.EGs,.GlobalEnv)
	assign("all_variance.DEGs",all_variance.DEGs,.GlobalEnv)
	assign("all_variance.EGs",all_variance.EGs,.GlobalEnv)
	assign("DEGsNb",DEGsNb,.GlobalEnv)
	assign("vsdDEGs",vsdDEGs,.GlobalEnv)

};
pca_inputData.function(vsd,data.all.deseq.fcros,name);
save(vsdDEGs,all_variance.EGs,all_variance.EGs,DEGsNb,PCA1.DEGs,PCA2.DEGs,PCA3.DEGs,
	PCA1.EGs,PCA2.EGs,PCA3.EGs,PC1A.DEGs,PC2A.DEGs,PC3A.DEGs,PC1A.EGs,PC2A.EGs,PC3A.EGs, 
	file=paste0(wd, f_results, f_pca,f_Rdata,name,"pca_components.RData"));
##########################################################

#TgDyrk1a:6 samples
mycol2  <- colorRampPalette(rev(brewer.pal(9, "RdPu")))(11)[3:8] #MUT color
mycol1  <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(11)[2:7] #WT color
#mycol <- union(mycol2,mycol1)
pcaDataColNamesEGs <-data.frame(PCA1=PCA1.EGs,PCA2=PCA2.EGs,PCA3=PCA3.EGs)
pcaDataColNamesEGs$labels <- factor(rownames(pcaDataColNamesEGs),levels=rownames(pcaDataColNamesEGs));
pcaDataColNamesEGs$colPalette <- ifelse(names(PCA1.EGs)=="Wt1", mycol1[1],
	ifelse(names(PCA1.EGs)=="Wt2", mycol1[2],
	ifelse(names(PCA1.EGs)=="Wt3", mycol1[3],
	ifelse(names(PCA1.EGs)=="Wt4", mycol1[4],
	ifelse(names(PCA1.EGs)=="Wt5", mycol1[5],
	ifelse(names(PCA1.EGs)=="Wt6", mycol1[6],
	ifelse(names(PCA1.EGs)=="Dup1", mycol2[1],
	ifelse(names(PCA1.EGs)=="Dup2", mycol2[2],
	ifelse(names(PCA1.EGs)=="Dup3", mycol2[3],
	ifelse(names(PCA1.EGs)=="Dup4", mycol2[4],
	ifelse(names(PCA1.EGs)=="Dup5", mycol2[5],
	ifelse(names(PCA1.EGs)=="Dup6", mycol2[6],
		 NA))))))))))))
pcaDataColNamesDEGs <-data.frame(PCA1=PCA1.DEGs,PCA2=PCA2.DEGs,PCA3=PCA3.DEGs)
pcaDataColNamesDEGs$labels <- factor(rownames(pcaDataColNamesDEGs),levels=rownames(pcaDataColNamesDEGs));
pcaDataColNamesDEGs$colPalette <- ifelse(names(PCA1.DEGs)=="Wt1", mycol1[1],
	ifelse(names(PCA1.DEGs)=="Wt2", mycol1[2],
	ifelse(names(PCA1.DEGs)=="Wt3", mycol1[3],
	ifelse(names(PCA1.DEGs)=="Wt4", mycol1[4],
	ifelse(names(PCA1.DEGs)=="Wt5", mycol1[5],
	ifelse(names(PCA1.DEGs)=="Wt6", mycol1[6],
	ifelse(names(PCA1.DEGs)=="Dup1", mycol2[1],
	ifelse(names(PCA1.DEGs)=="Dup2", mycol2[2],
	ifelse(names(PCA1.DEGs)=="Dup3", mycol2[3],
	ifelse(names(PCA1.DEGs)=="Dup4", mycol2[4],
	ifelse(names(PCA1.DEGs)=="Dup5", mycol2[5],
	ifelse(names(PCA1.DEGs)=="Dup6", mycol2[6],
		 NA))))))))))))

##########################################################
#1 creating the 3D plots #################################
#library("scatterplot3d")
##########################################

nameLine <- "Dyrk1a"
order <- data.frame(labels=c(paste0(rep("Wt",times=5),1:5),paste0(rep("Dyrk1a+/- ",times=5),1:5)))
pcaDataColNamesEGs <- left_join(order,pcaDataColNamesEGs,by="labels")
pcaDataColNamesEGs <-pcaDataColNamesEGs[,c(2,3,4,1,5)]
rownames(pcaDataColNamesEGs) <- gsub(" ","_", pcaDataColNamesEGs$labels)

min <- -150
max<- 120

min1 <- -20
max1<- 20

x11(height=7, width=9)
m <- rbind(c(1,1,3,3,5),
           c(2,2,4,4,6))
layout(m)

par(mar=c(3.5, 2, 3.5, 2))
scatterplot3d(pcaDataColNamesEGs$PCA1,pcaDataColNamesEGs$PCA2,pcaDataColNamesEGs$PCA3,       
                color=pcaDataColNamesEGs$colPalette, pch=16, 
                type="p",cex.symbols=3.5,angle=35, col.axis="grey28", bg="grey81",lty.hide=2, lab = c(2,2),
                 lab.z = 2,cex.axis=0.8,cex.lab=0.9,y.margin.add=0.1,
                xlim=c(min,max), ylim=c(min,max), zlim=c(min,max),
                col.grid="lightblue", 
                main=paste0(nameLine," \n (",dim(assay(vsd))[1]," EGs)"),
                xlab=paste0("PC1 (",PC1A.EGs,"%)"),
                ylab=paste0("                           PC2 (",PC2A.EGs,"%)"),
                zlab=paste0("PC3 (",PC3A.EGs,"%)"))

frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
# Table grob
grid.draw(table1a)
popViewport(3)


scatterplot3d(pcaDataColNamesDEGs$PCA1,pcaDataColNamesDEGs$PCA2,pcaDataColNamesDEGs$PCA3,       
                color=pcaDataColNamesDEGs$colPalette, pch=16, 
                type="p",cex.symbols=3.5,angle=35, col.axis="grey28", bg="grey81",lty.hide=2, 
                lab = c(2,2), lab.z = 2,cex.axis=0.8,cex.lab=0.9,y.margin.add=0.1,
                xlim=c(min1,max1), ylim=c(min1,max1), zlim=c(min1,max1),
                col.grid="lightblue", 
                main=paste0(nameLine,"\n (",DEGsNb," DEGs)"),
                xlab=paste0("PC1 (",PC1A.DEGs,"%)"),
                ylab=paste0("                           PC2 (",PC2A.DEGs,"%)"),
                zlab=paste0("PC3 (",PC3A.DEGs,"%)"))

frame()
# Grid regions of current base plot (ie from frame)
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
# Table grob
grid.draw(table2a)
popViewport(3)

frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.draw(NULL)
legend("center", legend =gsub("wt","WT",gsub(".*_"," ",pcaDataColNamesEGs$labels)), pch = 16,
 col = pcaDataColNamesEGs$colPalette,pt.cex=2.5, cex=1, inset=c(0.3))
popViewport(3)
dev.copy2pdf(device = quartz, file = paste0( wd,f_results,f_pca, name, "_3dpca.pdf"),onefile=TRUE)
dev.off()

##########################################################
#2. creating 2D PCA plots ################
#library("gridGraphics");
#library("cowplot");
#library("ggpubr");
#library("ggplot2");
##########################################
######################################################################
# NOTE: be CAREFUL, the variable labels in color=labels, 
# have to be a factor with the groups/categories/genoypes 
# ordered in the same order as the rows in the df 
# pcaDataColNamesDEGsif not the colors are asigned not oganised!!!!!
######################################################################

pc_qp1<- qplot(x=PCA1, y=PCA2, data=as.data.frame(pcaDataColNamesEGs), color=labels)+ 
theme_gray(base_size = 11) + theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",
legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + 
geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
scale_color_manual(breaks = rownames(pcaDataColNamesEGs), values=pcaDataColNamesEGs$colPalette)+ 
labs(x = paste0("PCA1 (",PC1A.EGs, "% )"), y=paste0("PCA2 (",PC2A.EGs, "% )"));  

pc_qp2<-qplot(x=PCA1, y=PCA3, data=as.data.frame(pcaDataColNamesEGs), color=labels) + 
theme_gray(base_size = 11)+ theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
	legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",
	legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) +
	 geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+
	  scale_color_manual(breaks = rownames(pcaDataColNamesEGs), values=pcaDataColNamesEGs$colPalette)+ labs(x = paste0("PCA1 (",PC1A.EGs, "% )"), y=paste0("PCA3 (",PC3A.EGs, "% )")); 

pc_qp3<-qplot(x=PCA2, y=PCA3, data=as.data.frame(pcaDataColNamesEGs), color=labels) + 
theme_gray(base_size = 11)+ theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
	legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",
	legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + 
geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
scale_color_manual(breaks = rownames(pcaDataColNamesEGs), values=pcaDataColNamesEGs$colPalette)+  
labs(x = paste0("PCA2 (",PC2A.EGs, "% )"), y=paste0("PCA3 (",PC3A.EGs, "% )")); 

pc_qp4<-qplot(x=PCA1, y=PCA2, data=as.data.frame(pcaDataColNamesDEGs), color=labels) + 
theme_gray(base_size = 11)+ theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
	legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",
	legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + 
geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
scale_color_manual(breaks = rownames(pcaDataColNamesDEGs), values=pcaDataColNamesDEGs$colPalette)+  
labs(x = paste0("PCA1 (",PC1A.DEGs, "% )"), y=paste0("PCA2 (",PC2A.DEGs, "% )")); 

pc_qp5<-qplot(x=PCA1, y=PCA3, data=as.data.frame(pcaDataColNamesDEGs), color=labels) + 
theme_gray(base_size = 11)+ theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
	legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",
	legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + 
geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
scale_color_manual(breaks = rownames(pcaDataColNamesDEGs), values=pcaDataColNamesDEGs$colPalette)+ 
labs(x = paste0("PCA1 (",PC1A.DEGs, "% )"), y=paste0("PCA3 (",PC3A.DEGs, "% )")); 

pc_qp6<-qplot(data=pcaDataColNamesEGs, x=PCA2, y=PCA3,color=labels )+
	 	 geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") +
	 	  geom_vline(xintercept = 0, colour = "gray65")+ 
	 	  scale_color_manual(values=pcaDataColNamesEGs$colPalette,labels=pcaDataColNamesEGs$labels)+
	 	    labs(x = paste0("PCA2 (",PC2A.DEGs, "% )"), 
	 	    	y=paste0("PCA3 (",PC3A.DEGs, "% )")) +
	 	    guides(colour = guide_legend("Condition",nrow = 2, byrow = TRUE))+ 
	 	    theme_gray(base_size = 12)+
	 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
	 	legend.title = element_text(size=11, face="bold"),legend.position="bottom", 
	 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
	 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));


#labels
upText <-paste0("   rat Dyrk1a PCA (",dim(assay(vsd))[1],") EGs");
upText.p <- ggparagraph(text = upText, face = "bold", size = 11, color = "black");
midText <-paste0("   rat Dyrk1a  PCA (",DEGsNb,") DEGs");
midText.p <- ggparagraph(text = midText, face = "bold", size = 11, color = "black");
# Arrange the plots on the same page
x11(height=8, width=8) #good dimensions
legend <- cowplot::get_legend(pc_qp6); 
blank<-rectGrob(gp=gpar(col="white")); # make a white spacer grob
p1 <- ggarrange(upText.p, ncol = 1, nrow = 1);
p2 <- ggarrange(pc_qp1+  theme(legend.position="none"), pc_qp2+  theme(legend.position="none"), pc_qp3+  
	theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));
p12 <- ggarrange(p1, p2, ncol = 1, nrow = 2,heights = c(0.1, 1));
p3 <- ggarrange(midText.p, ncol = 1, nrow = 1);
p4 <- ggarrange(pc_qp4+  theme(legend.position="none"), pc_qp5+  theme(legend.position="none"), pc_qp6 +  
	theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));
p34 <- ggarrange(p3, p4, ncol = 1, nrow = 2,heights = c(0.1, 1));
ggarrange(blank,p12,blank,p34,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5));
dev.copy2pdf(device = quartz, file = paste0( wd,f_results,f_pca, name, "_2dpca.pdf"),onefile=TRUE,height=5.5, width=8);

################################################################################################
################################################################################################
#
# Distance plots: 
# get the list of genes associated to each disregulated pathway
#
################################################################################################

######################
#2. Distance similarity matrix
######################
# A useful first step in an RNA-seq analysis is often to assess overall similarity between samples:
# Which samples are similar to each other, which are different? Does this fit to the expectation 
# from the experiment’s design?

# We use the R function dist to calculate the Euclidean distance between samples. To ensure we 
# have a roughly equal contribution from all genes, we use it on the VST data. We need to transpose
# the matrix of values using t, because the dist function expects the different samples to be rows 
#of its argument, and different dimensions (here, genes) to be columns.

#folder for output
f_dist <- "SamplesDistance/"
dir.create(file.path(getwd (), f_results, f_dist), showWarnings = F)

## For deseq2 data
#####################################
#####################################

sampleDistance_deseq2.function <-  function(vsd, metafile,vsdDEGs,dds) {	
	# Euclidean distances 
	#####################################
	#####################################
	sampleDistsEGs <- dist(t(assay(vsd)))
	sampleDistsEGsMatrix <- as.matrix(sampleDistsEGs) # Formatting the results of the previous step
	rownames(sampleDistsEGsMatrix) <- colnames(sampleDistsEGsMatrix)

	sampleDistsDEGs <- dist(t(vsdDEGs))
	sampleDistsDEGsMatrix <- as.matrix(sampleDistsDEGs) # Formatting the results of the previous step
	rownames(sampleDistsDEGsMatrix) <- colnames(sampleDistsDEGsMatrix)

	# Poisson Distance
	#####################################
	#####################################
	#is to use the Poisson Distance (Witten 2011), implemented 
	# in the PoiClaClu package. This measure of dissimilarity 
	# between counts also takes the inherent variance structure
	# of counts into consideration when calculating the distances 
	# between samples. The PoissonDistance function takes the 
	# original count matrix (not normalized) with samples as rows
	# instead of columns, so we need to transpose the counts in dds
	# library("PoiClaClu")

	poisd.EGs <- PoissonDistance(t(counts(dds)),type="deseq")
	samplePoisDistMatrix.EGs <- as.matrix( poisd.EGs$dd )
	colnames(samplePoisDistMatrix.EGs) <- colnames(counts(dds))
	rownames(samplePoisDistMatrix.EGs) <- colnames(samplePoisDistMatrix.EGs)

	ddsDEGs <- counts(dds)[which(rownames(counts(dds)) %in% rownames(vsdDEGs)),]
	ddsDEGsNb<-dim(ddsDEGs)[1]
	poisd.DEGs <- PoissonDistance(t(ddsDEGs),type="deseq")
	samplePoisDistMatrix.DEGs <- as.matrix( poisd.DEGs$dd )
	colnames(samplePoisDistMatrix.DEGs) <- colnames(counts(dds))
	rownames(samplePoisDistMatrix.DEGs) <- colnames(samplePoisDistMatrix.DEGs)

	assign("samplePoisDistMatrix.EGs",samplePoisDistMatrix.EGs,.GlobalEnv)
	assign("samplePoisDistMatrix.DEGs",samplePoisDistMatrix.DEGs,.GlobalEnv)
	assign("ddsDEGsNb",ddsDEGsNb,.GlobalEnv)
	assign("sampleEuclDistsEGsMatrix",sampleDistsEGsMatrix,.GlobalEnv)
	assign("sampleEuclDistsDEGsMatrix",sampleDistsDEGsMatrix,.GlobalEnv)
};

#run the function
sampleDistance_deseq2.function(vsd,colData,vsdDEGs,dds); 

#####################################
#####################################
# Plot Distance similarity bw samples
#####################################

col1Dis <- colorRampPalette(rev(brewer.pal(9, "RdPu")))(100)  

#for 2 genotypes 
dist.heatmapsDESeq2.plot.function <- function(sampleEuclDistsEGsMatrix,samplePoisDistMatrix.EGs,
	sampleEuclDistsDEGsMatrix,samplePoisDistMatrix.DEGs,col1) {
	#Plotting the sample similarity heatmaps with the euclidean 
	# and poisson distance computed with type "deseq".  
	x11(width=3, height=4)
	par(cex.main=0.8)
	grab_grob.function <- function(){
  		grid.echo()
  		grid.grab()
	};

	gplots::heatmap.2(sampleEuclDistsEGsMatrix, trace="none",key=FALSE,cexCol=0.8,cexRow=0.8, 
		dendrogram="column", scale="row",  col = col1, main= " Euclidean EGs",margin=c(6,6))
	p1 <- grab_grob.function()
	gplots::heatmap.2(samplePoisDistMatrix.EGs, trace="none", key=FALSE,cexCol=0.8,cexRow=0.8,
		dendrogram="column", scale="row",  col = col1,  main= " Poisson EGs",margin=c(6,6))
	p2 <- grab_grob.function()
	gplots::heatmap.2(sampleEuclDistsDEGsMatrix, trace="none",key=FALSE,cexCol=0.8,cexRow=0.8, 
		dendrogram="column", scale="row",  col = col1, main= " Euclidean DEGs",margin=c(6,6))
	p3 <- grab_grob.function()
	gplots::heatmap.2(samplePoisDistMatrix.DEGs, trace="none",key=FALSE,cexCol=0.8,cexRow=0.8, 
		dendrogram="column", scale="row",  col = col1, main= " Poisson DEGs",margin=c(6,6))
	p4 <- grab_grob.function()
	pAll<- grid.arrange(p1,p2,p3,p4, ncol=2,clip=TRUE)
	pdf(file=paste0(wd, f_results, f_dist, name,"_SampleDistance_blindF.pdf", sep=""))
	grid.draw(pAll)
	gplots::heatmap.2(sampleEuclDistsEGsMatrix, trace="none",  dendrogram="column", scale="row",  
		col = col1, main= " Euclidean EGs",margin=c(12,12))
	gplots::heatmap.2(samplePoisDistMatrix.EGs, trace="none",  dendrogram="column", scale="row",  
		col = col1,  main= " Poisson EGs",margin=c(12,12))
	gplots::heatmap.2(sampleEuclDistsDEGsMatrix, trace="none",  dendrogram="column", scale="row", 
	 col = col1, main= " Euclidean DEGs",margin=c(12,12))
	gplots::heatmap.2(samplePoisDistMatrix.DEGs, trace="none",  dendrogram="column", scale="row", 
	 col = col1, main= " Poisson DEGs",margin=c(12,12))
	dev.off()
};
dist.heatmapsDESeq2.plot.function(sampleEuclDistsEGsMatrix,samplePoisDistMatrix.EGs,sampleEuclDistsDEGsMatrix,
	samplePoisDistMatrix.DEGs,col1Dis)


################################################################################################
########################################################################################################################################
#
# Volcano plots: 
# get the list of genes associated to each disregulated pathway
#
########################################################################################################################################
################################### colours and groups
#palette:
mypalVolcan <- c("azure4", "purple1", "darkgreen", "deeppink")
mypalVolcan <- setNames(mypalVolcan, c("NotSignificant", "Significant", "FC", "Significant&FC"))

OI <- c("NotSignificant"="azure4", "Significant"="purple1","FC"="darkgreen","Significant&FC"="deeppink")
legend_title <- "Legend"
library(ggrepel)


# GENERATING THE PLOTS
###############################
###############################
# We need 3 types of  code, one to get the individual plots of each ML. ANd another two to get the plots necesary for the paper
# To plot all of them together wo innecessary things and adding the legend "apart" we need two plots:
# a) plotVolc : for all ML except the last one as it will carry the legend
# b) plotEND: is a plot as the others but carrying the legend on the right
# c) white: the thrid plot will be an empty plot

# STATS UNDER THE LABELING
#######################################
#Not significant - mean that fcros did not select these genes as DEGs based on the alpha cut off
#Significant:  fcros selected these genes based on the alpha cut off, and all of these genes have a fdr< 0.05
#FC: of all the expressed genes, the ones that have a abs(FC)> 1.4
#######################################

volcanoPlot.function<- function(name, dataEGs,dataDEGsFcros,width,n){
	#README: function to create volcano plots.
	# The function is practically same for human Hgu133a and RNASEQ and 
	# our own Gen 1st microarrays. The only change for the Hgu133a data is that label=SYMBOL
	# instead for all RNASeq data and our own microarrays label=swissprot.
	f_volcano<- "volcanoPlot/"
	dir.create(file.path(getwd (), f_results, f_volcano), showWarnings = FALSE)


	op <- par( cex=0.5)
	condFC2 <- abs(dataEGs$FC2) > 1.4
	#condFDR <- dataDEGsFcros[,4] < 0.05 #not needed the values are alway smaller
	num <- dim(dataEGs)[2] +1
	dataEGs$Group <- "NotSignificant"
	dataEGs[which(dataEGs[,1] %in% dataDEGsFcros[,1]),num] <- "Significant" 
	dataEGs[which(dataEGs[,4]> min(dataDEGsFcros$fdr_fcros)  &  condFC2),num] <- "FC" 
	dataEGs[which(dataEGs[,1] %in% dataDEGsFcros[,1]),][which(abs(dataEGs[which(dataEGs[,1] %in% dataDEGsFcros[,1]),'FC2'])>1.4),num]   <- "Significant&FC" 
	dataEGs <- dataEGs[order(dataEGs$FC2, dataEGs$fdr, decreasing=T),] #first fc and then fdr
	dataEGs <- dataEGs[!duplicated(dataEGs[,2]), ]

	#label top genes
	
	top_dataEGs <- dataEGs[order(dataEGs$FC2, dataEGs$fdr, decreasing=T),][1:n,] #first fc and then fdr
	top_dataEGs <- rbind(top_dataEGs, dataEGs[order(-dataEGs$FC2, dataEGs$fdr, decreasing=T),][1:n,])
	top_dataEGs <- top_dataEGs[!duplicated(top_dataEGs[,2]), ]


	#1. individual plot
	#######################
	IndivPlot = ggplot(data=dataEGs, aes(FC2, -log10(fdr), colour=Group)) +
	  geom_point(aes(size=3),alpha=0.4) + 
	  scale_color_manual(values=mypalVolcan) + scale_size(guide="none")+
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
	  	legend.title = element_text(size=14, face="bold"),legend.position="bottom", 
	  	legend.box = "horizontal",legend.text=element_text(size=13),legend.key.width=unit(1.6,"line"),
	  	legend.key.height=unit(1.1,"line")) +
	  labs(y="-log10 FDR", 
	       x="log2 Fold change", 
	       title=paste("Distribution of differentially expressed genes in ",name ),
	       caption = "Source: Y.Herault team @IGBMC") + guides(fill=FALSE) + 
	  guides(colour = guide_legend(override.aes = list(size=4))) 
	IndivPlot1<- IndivPlot + geom_text_repel(data=top_dataEGs[1:n,], aes(label=swissprot), 
		colour="navy",show.legend = TRUE)

	#IndivPlot
	pdf(file=paste0(wd, f_results, f_volcano, name, "_volcano.pdf"),height=8.27, width=width)
	plot(IndivPlot1)
	dev.off();

	#to combine plots
	####################
	plotVolc = ggplot(data=dataEGs, aes(FC2, -log10(fdr), colour=Group)) +
	  geom_point(aes(size=8),alpha=0.4) + 
	  scale_color_manual(values=mypalVolcan) + scale_size(guide="none")+
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"),
	  	legend.title = element_text(size=14, face="bold"),legend.position="right", legend.box = "horizontal",
	  	legend.text=element_text(size=13),legend.key.width=unit(1.6,"line"),legend.key.height=unit(1.1,"line")) +
	  labs(y="-log10 FDR", 
	       x="log2 Fold change", 
	       title=name) + guides(fill=FALSE) + guides(colour = guide_legend(override.aes = list(size=4)))+ 
	  theme(legend.position="none") + 
	  scale_size(range = c(1, 1))
	#plotVolc
	#plotVolc+geom_text_repel(data=top_dataEGs[1:n,], aes(label=swissprot), colour="navy",show.legend = FALSE)

	Dp3PLot <- plotVolc+geom_text_repel(data=top_dataEGs[1:n,], aes(label=swissprot), colour="navy",show.legend = TRUE)

	Dp3PLotL <-Dp3PLot + xlim(-1.6, 9)+ ylim(-1, 40)
	assign("Dp3PLotL",Dp3PLotL,.GlobalEnv);
	assign("Dp3PLot",Dp3PLot,.GlobalEnv);

	pdf(file=paste0(wd, f_results, f_volcano, name, "_volcanob.pdf"))
	plot(Dp3PLotL)
	dev.off();
};

n <- 15
volcanoPlot.function("Dyrk1a", afDf,data.all.deseq.fcros,12,n)

#####################################################################################################
############################### END. Finished on the  190218 ########################################
###################  Mar Muniz. PhDstudent Yann Herault lab. @IGBMC #################################
#####################################################################################################

