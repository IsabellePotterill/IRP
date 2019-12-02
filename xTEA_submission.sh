#!/bin/sh 

#xTEA example submission scripts

#to make sample names
ls 24*/DNA/final.bam_GRCh38/*.bam | sed -e 's:\.final.bam::' | sed -e 's:[0-9]*/DNA/final.bam_GRCh38/::' > /lustre/scratch119/casm/team176tv/ip10/xTEA/hg38/sample.name.txt

#to make bam list
ls /lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/24*/DNA/final.bam_GRCh38/*.bam |sed -r -e 's:[^$]*:& &:' | sed -e 's:\.final.bam::' -e 's:/lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/[0-9]*/DNA/final.bam_GRCh38/::' > /lustre/scratch119/casm/team176tv/ip10/xTEA/hg38/illumina_bam_list.txt

#make sh files for each cell
python /nfs/users/nfs_i/ip10/xTEA/gnrt_pipeline_local_v38.pyc -i sample.name.txt -b illumina_bam_list.txt -x null -p /lustre/scratch119/casm/team176tv/ip10/xTEA -o xtea.sh -n 20 -l /nfs/users/nfs_h/hj6/software/xTEA/rep_lib_annotation/ -r /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_ERCC92/genome.fa --xtea /nfs/users/nfs_h/hj6/software/xTEA/ --flklen 3000 -y 1 -g /nfs/users/nfs_h/hj6/software/xTEA/hg19_gene.txt --cr 1 --nclip 1 --nd 1 --nfclip 1 --nfdisc 1


# run each sh file 
for sample in SCGC*;do bsub -q normal -G team176 -o n_log3.txt -n 8 -e n_err3.txt -R"select[mem>60000] rusage[mem=60000] span[hosts=1]" -M60000 "sh /lustre/scratch119/casm/team176tv/ip10/xTEA/high_cov_2/$sample/L1/run_xTEA_pipeline.sh";done