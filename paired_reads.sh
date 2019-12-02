#!/bin/sh 
# a file to get unmapped and mapped pairs from bam files, compare to L1 non-redundant list and blast unmap pairs to L1 database


samtools view -u -f 8 -F 260 *.sorted.bam  > oneEndMapped.bam #-f 8 mate unmapped, -F 260 skip when itself is unmapped

samtools view -u  -f 4 -F 264 *sorted.bam  > oneEndUnmapped.bam # -f 4 any read unmapped, -F 264 skip entry where mate is unmapped

#intersect the mapped reads with L1 flank regions
bedtools intersect -a oneEndMapped.bam -b /nfs/users/nfs_i/ip10/reference/L1_final.bed -wa > mapped_to_flank.bam

#get read IDs from mathces to L1 flank
samtools view mapped_to_flank.bam | cut -f 1 > read_IDs

#get unmapped pairs to the matches using read IDS
samtools view oneEndUnmapped.bam | grep -f read_IDs > L1_unmaped_reads.bam

#get fasta sequence
cut -f -1,10 L1_unmaped_reads.bam > fasta_att

#blast reads against blast database
blastall -d /lustre/scratch119/casm/team176tv/ip10/local_align/l1_ref/L1_DATABASE/L1_database.fasta -i final.fasta -p blastn -m 8 -o output.blast