#!/bin/sh
#SBATCH -p igbmc                     # partition
#SBATCH -o slurm.%N.%j.out           # STDOUT
#SBATCH -e slurm.%N.%j.err           # STDERR
##SBATCH --mail-type=ALL              # can be BEGIN, END, FAIL or REQUEUE
#SBATCH --mail-user=munizmom@igbmc.fr

#################################################################################################
###################################### README  ##################################################
#
# Script to generate kallisto Transcripts counts for RNASeq data.
# Done and maintained by Mar Muniz Moreno PhDstudent YHerault team @IGBMC
#################################################################################################
#################################################################################################
#################################################################################################
##to run it on the cluster:
##########################
#
#sbatch --mem=80GB --cpus=16 -J klstT  -v -v kallisto_Transcripts.sh 180603 180603/S18057 mus_musculus grcm38.93 Mus_musculus.GRCm38.cdna.all.fa
#################################################################################################
#
#README:
#################################################################################################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";

echo "=========================================================="
echo "Starting on : ${current_date_time}"
echo "Current directory : $(pwd)"
echo "=========================================================="


name3=$3
name1=$1
name2=$2
PrefixNameIndex=$4
TranscriptomeFile=$5
transcriptome=/shared/misc/data5/herault_eq/munizmom/Ensembl/EMBL/ftp.ensembl.org/pub/release-96/embl/${name3}/transcriptome
results=/shared/space2/herault/OmicsAnalyses/omicsAnalyses/ongoing/projects/${name3}/${name1}/
rawData=/shared/space2/herault/RNASeq/${name3}/${name2}/raw/gz

echo "raw data of : ${name3} , ${name2} " ;
echo "output in : ${name3}, ${name1}" ;
echo "=========================================================="

mkdir ${results}/
mkdir ${results}/kallisto_Trancripts

# Quantification of the  abundances of the transcripts using  Kallisto (wo alignment)
######################################################################################## 

#echo kallisto index -i ${rawData}/${i}.idx ${rawData}/${i}.fastq.gz 
#kallisto index --make-unique -i ${transcriptome}/${PrefixNameIndex}_kallisto.idx ${transcriptome}/${TranscriptomeFile}

for i in `ls ${rawData}/*.fastq.gz |awk '{print $NF}' FS=/|sed -e 's/\.fastq.gz//g' |  sort -u`
do
kallisto quant --single --bootstrap-samples=100 --sd=50 --fragment-length=200  -o ${results}kallisto_Trancripts/${i} -i ${transcriptome}/${PrefixNameIndex}_kallisto.idx  ${rawData}/${i}.fastq.gz  
done
######################################################################################## 

echo ${name}
echo "************************"
echo "Finished on : $(date)"

