#!/bin/bash
#PBS -q q512G
#PBS -l mem=100gb,nodes=1:ppn=3,walltime=1000:00:00 
#HSCHED -s Singlecell+cellranger+human
#PPN limit 3

sample="t1219"
mkdir /asnas/ciwm_group/zhangqi/scRNA/Cellranger/${sample}
cd /asnas/ciwm_group/zhangqi/scRNA/Cellranger/${sample}

bin=/xtdisk/ciwm_group/tanyzh/Mouse/software/cellranger-6.0.0/cellranger
ref=/asnas/ciwm_group/zhangqi/ref
opt=/asnas/ciwm_group/zhangqi/scRNA/

${bin} count --id=t1219 \
           --transcriptome=${ref}/refdata-gex-GRCh38-2020-A \
           --fastqs=${opt}/t1219 \
           --sample=t1219 \
           --localcores=16 \
           --localmem=64