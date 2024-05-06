dir=/data5/zhangq/scRNA/Molecular/Validation/Revised/Tumor/pyscenic
tfs=/data5/zhangq/scRNA/MultiSample/NCTMerge9/pySCENICtest/hs_hgnc_tfs.txt
feather=/data5/zhangq/scRNA/MultiSample/NCTMerge9/pySCENICtest/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=/data5/zhangq/scRNA/MultiSample/NCTMerge9/pySCENICtest/motifs-v9-nr.hgnc-m0.001-o0.0.tbl 

input_loom=/data5/zhangq/scRNA/Molecular/Validation/Revised/Tumor/pyscenic/PTTT_DTTEpi_Annotation.loom
ls $tfs  $feather  $tbl  

# pyscenic grn
/data5/zhangq/miniconda3/envs/pyscenic/bin/pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
$input_loom  \
$tfs 

#pyscenic cistarget
/data5/zhangq/miniconda3/envs/pyscenic/bin/pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20  \
--mask_dropouts

#pyscenic AUCell
/data5/zhangq/miniconda3/envs/pyscenic/bin/pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 20 