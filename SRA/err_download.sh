#!/bin/bash
#$ -S /bin/bash
#$ -q regular.q@cbsubscb01,regular.q@cbsubscb02,regular.q@cbsubscb03,regular.q@cbsubscb04,regular.q@cbsubscb05,regular.q@cbsubscb06,regular.q@cbsubscb07,regular.q@cbsubscb08,regular.q@cbsubscb10,regular.q@cbsubscb11,regular.q@cbsubscb12,regular.q@cbsubscb13
#$ -j y
#$ -cwd
#$ -pe bscb 2
#$ -t 1-197:1
#$ -N fa2_download
#$ -l h_rt=06:00:00

d1=$(date +%s)

newdir=$JOB_ID${SGE_TASK_ID}

echo $HOSTNAME
echo ${SGE_TASK_ID}
echo $1
echo $newdir

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir

### change the path to this folder accordingly:
cp $HOME/home/ya76/GitHub_Repositories/FlyAtlas2.RNAseq/SRA/accesions.list .

SRA=$(awk "NR==$SGE_TASK_ID" accesions.list)

fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $SRA

### change the path to this folder accordingly:
mv fastq/* /fs/cbsufsrv5/data2/ya76/FlyAtlas2/READS/ERR
cd ..
rm -r ./$newdir

