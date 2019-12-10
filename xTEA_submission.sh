#!/bin/sh 

#xTEA example submission scripts
#put sample.name.txt and illumina_bam_list.txt into directory that you will then run python command

#to put sample names to sample.name.txt file
#removing .final.bam as well as path to leave identifier e.g. NA12878

ls /lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/24*/DNA/final.bam_GRCh38/*.bam | sed -e 's:\.final.bam::' | sed -e 's:/lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/[0-9]*/DNA/final.bam_GRCh38/::' > /lustre/scratch119/casm/team176tv/ip10/xTEA/hg38/sample.name.txt

#to put bam list into illumina_bam_list.txt
#put sample ID and path into bam file e.g. NA12878 /path/na12878_illumina_1_sorted.bam

ls /lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/24*/DNA/final.bam_GRCh38/*.bam |sed -r -e 's:[^$]*:& &:' | sed -e 's:\.final.bam::' -e 's:/lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/[0-9]*/DNA/final.bam_GRCh38/::' > /lustre/scratch119/casm/team176tv/ip10/xTEA/hg38/illumina_bam_list.txt

#make sh files for each cell 
#run in directory you have put sample.name.txt and illumina_bam_list.txt

python /nfs/users/nfs_i/ip10/xTEA/gnrt_pipeline_local_v38.pyc -i sample.name.txt -b illumina_bam_list.txt -x null -p /lustre/scratch119/casm/team176tv/ip10/xTEA -o xtea.sh -n 20 -l /nfs/users/nfs_h/hj6/software/xTEA/rep_lib_annotation/ -r /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_ERCC92/genome.fa --xtea /nfs/users/nfs_h/hj6/software/xTEA/ --flklen 3000 -y 1 -g /nfs/users/nfs_h/hj6/software/xTEA/hg19_gene.txt --cr 1 --nclip 1 --nd 1 --nfclip 1 --nfdisc 1


# run each sh file 

for sample in SCGC*;do bsub -q normal -G team176 -o n_log3.txt -n 8 -e n_err3.txt -R"select[mem>60000] rusage[mem=60000] span[hosts=1]" -M60000 "sh /lustre/scratch119/casm/team176tv/ip10/xTEA/high_cov_2/$sample/L1/run_xTEA_pipeline.sh";done

#look at output in candidate_disc_filtered_cns.txt in file/L1/ for each bam file
#column headings:
#chrm refined-pos lclip-pos rclip-pos TSD nalclip narclip naldisc nardisc nalpolyA narpolyA lcov rcov nlclip nrclip nldisc nrdisc nlpolyA nrpolyA lclip-cns-start:end rclip_cns-start:end ldisc-cns-start:end rdisc-cns-start:end Transduction-info ntrsdct-clip ntrsdct-disc ntrsdct-polyA ldisc-same ldisc-diff rdisc-same rdisc-diff 3mer-inversion confidential insertion-length
