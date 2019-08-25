#!/bin/sh
#SBATCH -p igbmc                     # partition
#SBATCH -o slurm.%N.%j.out           # STDOUT
#SBATCH -e slurm.%N.%j.err           # STDERR
#SBATCH --mail-type=ALL              # can be BEGIN, END, FAIL or REQUEUE
#SBATCH --mail-user=user@igbmc.fr
#sbatch --mem=40GB --cpus=16 -J ftpGet -v ftp.sh ftpaddress userRaw pswdRaw userAnalysed pswdAnalysed f_output
#################################################################################


current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";

echo "=========================================================="
echo "Starting on : ${current_date_time}" ;
echo "Current directory : $(pwd)"
echo "=========================================================="


f_output=$6
userRaw=$2
pswdRaw=$3
f_RawData=$1
userAnalysed=$4
pswdAnalysed=$5
###########################
#Defining the experiment directory
wd=/shared/space/herault/RNASeq/${f_output}

#creating the folders structure #######
mkdir ${wd}/raw
mkdir ${wd}/raw/GEOSub
mkdir ${wd}/raw/gz
mkdir ${wd}/raw/QC

mkdir ${wd}/SeqPlatFiles
mkdir ${wd}/SeqPlatFiles/reports
mkdir ${wd}/SeqPlatFiles/HTSeq-counts
mkdir ${wd}/SeqPlatFiles/bam
mkdir ${wd}/stouts

mkdir ${wd1}/raw
mkdir ${wd1}/raw/GEOSub
mkdir ${wd1}/raw/gz
mkdir ${wd1}/raw/QC

mkdir ${wd1}/SeqPlatFiles
mkdir ${wd1}/SeqPlatFiles/reports
mkdir ${wd1}/SeqPlatFiles/HTSeq-counts
mkdir ${wd1}/SeqPlatFiles/bam
mkdir ${wd1}/stouts

#######################################
# 1. Copying the raw data
######################################
cd ${wd}/raw/GEOSub/
wget -A txt,xls --random-wait  --ftp-user=${userRaw} --ftp-password=${pswdRaw} -m -p -E -k -K -np -nH --cut-dirs=100 ftp://ngs.igbmc.fr/${f_RawData}

cd ${wd}/raw/QC/
wget -A pdf --random-wait  --ftp-user=${userRaw} --ftp-password=${pswdRaw} -m -p -E -k -K -np -nH --cut-dirs=100 ftp://ngs.igbmc.fr/${f_RawData}

cd ${wd}/raw/gz/
wget -A gz --random-wait  --ftp-user=${userRaw} --ftp-password=${pswdRaw} -m -p -E -k -K -np -nH --cut-dirs=100 ftp://ngs.igbmc.fr/${f_RawData}


echo "************************"
echo "Finished on : $(date)"

#######################################
# 2. Copying the Analysed files
######################################

cd ${wd}/SeqPlatFiles/bam/
wget -A bam,bai --random-wait  --ftp-user=${userAnalysed} --ftp-password=${pswdAnalysed} -m -p -E -k -K -np -nH --cut-dirs=100  ftp://ngs.igbmc.fr/analyzeddata/

cd ${wd1}/SeqPlatFiles/HTSeq-counts/
wget -A tsv --random-wait  --ftp-user=${userAnalysed} --ftp-password=${pswdAnalysed} -m -p -E -k -K -np -nH --cut-dirs=100  ftp://ngs.igbmc.fr/analyzeddata/

cd ${wd1}/SeqPlatFiles/reports/
wget -A pdf --random-wait  --ftp-user=${userAnalysed} --ftp-password=${pswdAnalysed} -m -p -E -k -K -np -nH --cut-dirs=100  ftp://ngs.igbmc.fr/analyzeddata/

######################################
######################################

echo "************************"
echo "Finished on : $(date)"
