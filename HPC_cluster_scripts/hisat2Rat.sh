#!/bin/sh
#SBATCH -p igbmc                     # partition
#SBATCH -o slurm.%N.%j.out           # STDOUT
#SBATCH -e slurm.%N.%j.err           # STDERR
#SBATCH --mail-type=ALL              # 
#SBATCH --mail-user=user@igbmc.fr

#################################################################################################
###################################### README  ##################################################
# 
# Script to do all the steps of the HISAT2 pipeline for genome aligning 
# using HISAT2 and generation of ensembl ID gene counts using HTSeqCount
# Done and maintained by Mar Muniz Moreno PhDstudent YHerault team @IGBMC
#################################################################################################
#################################################################################################
#to run it on the cluster:                                                                      #
########################### 																	#
#module load htseq/0.9.1																		#
#module load samtools/1.9                                                                       #
#to run: 																						#
# sbatch --mem=64GB --cpus=16 -J hst2  -v -v hisat2Rat.sh DNML37_S1 180918_MaD02/S18030 S18030	# 
#################################################################################################
#################################################################################################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
name=$1
name2=$2
name3=$3

echo "=========================================================="
echo "raw data of :" ${name2} "  ," ${name} ;
echo "output in : " ${name3} ;
echo "Starting on :"  ${current_date_time} ;
echo "=========================================================="
echo "=========================================================="

#input files
genome=/shared/misc/data5/herault_eq/munizmom/Ensembl/EMBL/ftp.ensembl.org/pub/release-96/embl/rattus_norvegicus/primary_assembly
hisat2Index=/shared/misc/data5/herault_eq/munizmom/Ensembl/EMBL/ftp.ensembl.org/pub/release-96/embl/rattus_norvegicus/primary_assembly/hisat2Index
results=/shared/space2/herault/OmicsAnalyses/omicsAnalyses/ongoing/projects/rattus_norvegicus/${name3}/
rawData=/shared/space2/herault/RNASeq/rattus_norvegicus/${name2}/raw/gz
gtf=/shared/misc/data5/herault_eq/munizmom/Ensembl/EMBL/ftp.ensembl.org/pub/release-96/embl/rattus_norvegicus/genes



#mkdir ${results}/newJunctions
##alignment_hisat2
########################
#index building
#hisat2-build -p 16 ${wd}/${model}/primary_assembly/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa  ${wd}/${model}/primary_assembly/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa


## to run all the files of a folder consecutively
################################################################################################
#for i in `ls ${rawData}/*.fastq.gz |awk '{print $NF}' FS=/|sed -e 's/\.fastq.gz//g' |  sort -u`
#do
#echo ${i}.fastq.gz
#echo ${i}.bam
#echo #hisat2 -p --phred33  --dta -t --summary-file --novel-splicesite-outfile ${results}/newJunctions/${i}.junctions  -x ${genome}/Mus_musculus.GRCm38.dna.primary_assembly.fa  -U ${rawData}/${i}.fastq.qz -S ${results}/sam/${i}.sam
#hisat2 -p 4 --phred33  --dta -t --summary-file  -x ${genome}/Mus_musculus.GRCm38.dna.primary_assembly.fa  -U ${rawData}/${i}.fastq.gz -S ${results}/sam/${i}.sam
#done

################################################################################################
mkdir ${results}
mkdir ${results}/hisat2/
mkdir ${results}/hisat2/sam/
hisat2 -p 8 --phred33  --dta -t --summary-file --no-unal -x ${hisat2Index}/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
-U ${rawData}/${name}.fastq.gz -S ${results}/hisat2/sam/${name}.sam  

################################################################################################
#samtools faidx  ${wd}/${model}/primary_assembly/Mus_musculus.GRCm38.dna.primary_assembly.fa  #done in the genome folder
mkdir ${results}/hisat2/bam
samtools view -bS ${results}/hisat2/sam/${name}.sam > ${results}/hisat2/bam/${name}.bam
samtools view -h -F 4 -b ${results}/hisat2/bam/${name}.bam > ${results}/hisat2/bam/${name}_only_mapped.bam
samtools sort -m 2G -@ 8 -O bam -o ${results}/hisat2/bam/${name}.sorted.bam  ${results}/hisat2/bam/${name}_only_mapped.bam
samtools index  ${results}/hisat2/bam/${name}.sorted.bam 
################################################################################################

#running HTSeq-count
#gunzip  ${gtf}/Rnor_6.0.96.chr.gtf.gz
mkdir ${results}/hisat2/HTSeq/
samtools view ${results}/hisat2/bam/${name}.sorted.bam | htseq-count --mode=intersection-nonempty \
--stranded=no --type=exon --idattr=gene_id --additional-attr=gene_name -  ${gtf}/Rnor_6.0.96.chr.gtf \
> ${results}/hisat2/HTSeq/${name}.htseq.txt

################################################################################################

echo ${name2} ${name};
echo "************************"
echo "Finished on :" $(date);

