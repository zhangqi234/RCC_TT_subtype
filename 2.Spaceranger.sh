#!/bin/bash
#PBS -q q512G 
#PBS -l mem=150gb,nodes=1:ppn=6,walltime=1000:00:00 
#HSCHED -s Thrombus+spaceranger+human
#PPN limit 6

sample="2620-8"

bin=/home/tanyzh/xtdisk/ThrombusST/software/spaceranger-1.3.0/spaceranger
opt=/asnas/ciwm_group/zhangqi/ThrombusST
cd /home/tanyzh/xtdisk/ThrombusST

${bin} count --id=P26208 \
                --transcriptome=/home/tanyzh/xtdisk/ThrombusST/software/refdata-gex-GRCh38-2020-A \
                --probe-set=/home/tanyzh/xtdisk/ThrombusST/software/spaceranger-1.3.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \
                --fastqs=${opt}/2620-8 \
                --sample=2620-8 \
                --image=/asnas/ciwm_group/zhangqi/ThrombusST/image/2620-8/C-2620-8.tif \
                --slide=V11L12-013 \
                --area=C1 \
                --loupe-alignment=/asnas/ciwm_group/zhangqi/ThrombusST/image/2620-8/V11L12-013-2620-8-C1.json \
                --localcores=16 \
                --localmem=64 
