#!/bin/sh 
#a script to run BWA MEM on cram files

## run one at the time in their relevant directories
run_lanes="27128_1 27128_2 27128_3 27128_4 27128_5 27128_6" 
for run_lane in $run_lanes 
do
grep -A1 sample_supplier_name$ $run_lane/*.imeta | grep value | sed -e 's/.imeta-value: /.cram:/' \
> $run_lane.samples.txt
done


run_lanes="27128_1" # run this in the following path /lustre/scratch117/casm/team176/rr11/Human-brain-Alzheimer/cnv_validation
run_lanes="27128_2" #
run_lanes="27128_3" #
run_lanes="27128_4" #
run_lanes="27128_5" #
run_lanes="27128_6" 

## converting CRAM to FASTQ 
for run_lane in $run_lanes
do
mkdir $run_lane.fastq
for line in `cat $run_lane.samples.txt`
do
cram=`echo $line | sed -e 's/:.*//'`
sample=`echo $line | sed -e 's/.*://'`
bsub -o $run_lane.fastq/$sample.o -e $run_lane.fastq/$sample.e \
-G team176 \
-R'select[mem>=16000] rusage[mem=16000]' -M16000 \
/software/hpag/biobambam/latest/bin/bamtofastq \
inputformat=cram \
F=$run_lane.fastq/$sample.1.fastq \
F2=$run_lane.fastq/$sample.2.fastq \
filename=$cram
done
done

## Check jobs ran OK
ls */*.cram | wc -l && \
grep -l 'Successfully completed' *.fastq/*.o | wc -l && \
grep -l 'Exited' *.fastq/*.o | wc -l

## merge Fastq across lanes
mkdir fastq
for sample in `sed -e 's/.*://' *.samples.txt  | sort -u` #1.samples.txt is the list of samples except the runnig jobs above
do
for read in 1 2
do
bsub -o fastq/$sample.$read.merge.o -e fastq/$sample.$read.merge.e \
-G team176 \
"cat `ls *.fastq/$sample.$read.fastq | sort | tr '\n' ' '` > fastq/$sample.$read.fastq"
done
done

## check that my jobs run OK
ls */*.cram | wc -l && \
grep -l 'Successfully completed' fastq/*.o | wc -l && \
grep -l 'Exited' fastq/*.o | wc -l

## delete intermediate Fastq file
for run_lane in $run_lanes
do
rm -rf $run_lane.fastq
done

## BWA-MEM mapping 
mkdir BWA_mem_xtea_37
for sample in `ls fastq/*.1.fastq  | sed -e 's/\.1.fastq//' | sed -e 's/^fastq\///'`
do
bsub \
-n 14 -R"span[hosts=1]" \
-o BWA_mem_xtea_37/$sample.o -e BWA_mem_xtea_37/$sample.e \
-q long \
-G team176 \
-R 'select[mem>=30000] rusage[mem=30000]' -M30000 \
"/nfs/users/nfs_r/rr11/Tools/bwa/bwa mem -t 14 \
/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/genome.fa \
fastq/$sample.1.fastq \
fastq/$sample.2.fastq \
| /software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools view -bSho ./BWA_mem_xtea_37/$sample.xtea.mem.bam"
done

## sort bam
mkdir BWA_mem_xtea_sorted_37
for sample in `ls BWA_mem_xtea_37/*.bam | sed -e 's/\.xtea.mem.bam//' | sed -e 's/^BWA_mem_xtea_37\///'`
do
bsub \
-o BWA_mem_xtea_sorted_37/$sample.o -e BWA_mem_xtea_sorted_37/$sample.e \
-q normal \
-G team176 \
-R'select[mem>6000] rusage[mem=6000]' -M6000 \
"/software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools sort \
-o ./BWA_mem_xtea_sorted_37/$sample.sorted.bam ./BWA_mem_xtea_37/$sample.xtea.mem.bam"
done


## samtools index after sorting
for sample in `ls *.bam | sed -e 's/\.xtea.mem.bam//'`
do
bsub \
-o $sample.o -e $sample.e \
-q normal \
-G team176 \
-R'select[mem>4000] rusage[mem=4000]' -M4000 \
"/software/hgi/pkglocal/samtools-1.3.1-htslib-1.3.2/bin/samtools index $sample"
done

