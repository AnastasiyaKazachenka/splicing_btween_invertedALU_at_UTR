#!/bin/sh

#SBATCH --job-name=splice site
#SBATCH --ntasks=8 #number of cores
#SBATCH --nodes=1 # Ensure that all cores are on one machine 
#SBATCH --time=1-00:05 #time limit of zero requests no limit. Other formats are allowed 
#SBATCH --mem-per-cpu=4096 #memory per core in MB


#ICMT 
#SINE1 chr1:6222139-6222443
#SINE2 chr1:6223099-6223406
#addded 10bp to each side of sine
#3' UTR chr1:6221193-6225082


fullpath=$1   #collects string of input file from 1st argument (*.fastq)
fastq="${fullpath##*/}"
fname="${fastq%.bam}"

indirectory_full="${fullpath%/*}"
indirectory="${indirectory_full##*/}"

if [ ! -d "$indirectory" ]; then
	mkdir -p "$indirectory"
fi

cd "$indirectory"

module load SAMtools/1.8-foss-2016b
module load BEDTools/2.26.0-foss-2016b

#generate bam for locus
samtools view -b "$1" "$3" > sine1_"$fname"_"$5".bam
samtools index sine1_"$fname"_"$5".bam
samtools view -b "$1" "$4" > sine2_"$fname"_"$5".bam
samtools index sine2_"$fname"_"$5".bam

bedtools bamtobed -split -i sine1_"$fname"_"$5".bam > sine1_split_"$fname"_"$5".bed
bedtools bamtobed -i sine1_"$fname"_"$5".bam > sine1_nosplit_"$fname"_"$5".bed

echo "$3" > SINE1_coord_"$fname"_"$5".tmp.bed
echo "$4" >> SINE2_coord_"$fname"_"$5".tmp.bed
sed 's/:/\t/g' SINE1_coord_"$fname"_"$5".tmp.bed | sed 's/-/\t/g' > SINE1_coord_"$fname"_"$5".bed
rm SINE1_coord_"$fname"_"$5".tmp.bed
sed 's/:/\t/g' SINE2_coord_"$fname"_"$5".tmp.bed | sed 's/-/\t/g' > SINE2_coord_"$fname"_"$5".bed
rm SINE2_coord_"$fname"_"$5".tmp.bed
bedtools intersect -a sine1_split_"$fname"_"$5".bed -b SINE1_coord_"$fname"_"$5".bed -wa > sine1_bits_"$fname"_"$5".tmp.bed
bedtools intersect -a sine1_split_"$fname"_"$5".bed -b SINE2_coord_"$fname"_"$5".bed -wa > sine2_bits_"$fname"_"$5".tmp.bed
bedtools intersect -a sine1_bits_"$fname"_"$5".tmp.bed -b sine2_bits_"$fname"_"$5".tmp.bed -v | awk -F'\t' '$3-$2 >2' > sine1_bits_"$fname"_"$5".bed
bedtools intersect -a sine2_bits_"$fname"_"$5".tmp.bed -b sine1_bits_"$fname"_"$5".tmp.bed -v | awk -F'\t' '$3-$2 >2' > sine2_bits_"$fname"_"$5".bed
rm sine1_bits_"$fname"_"$5".tmp.bed
rm sine2_bits_"$fname"_"$5".tmp.bed
awk '{print $4}' sine2_bits_"$fname"_"$5".bed | grep -f - sine1_bits_"$fname"_"$5".bed | awk '{print $4}' | grep -f - sine1_nosplit_"$fname"_"$5".bed > sine1andsine2_"$fname"_"$5".tmp.bed
samtools view sine1_"$fname"_"$5".bam | awk '$6 ~ /N/ {print $1}' | grep -f - sine1andsine2_"$fname"_"$5".tmp.bed > sine1andsine2_"$fname"_"$5".bed
rm sine1andsine2_"$fname"_"$5".tmp.bed

samtools view -b "$1" "$2" > UTR_"$fname"_"$5".bam
samtools index UTR_"$fname"_"$5".bam
bedtools bamtobed -i UTR_"$fname"_"$5".bam > UTR_"$fname"_"$5".bed

echo "$2" > UTR_coord_"$fname"_"$5".tmp.bed
sed 's/:/\t/g' UTR_coord_"$fname"_"$5".tmp.bed | sed 's/-/\t/g' > UTR_coord_"$fname"_"$5".bed
rm UTR_coord_"$fname"_"$5".tmp.bed
bedtools intersect -a UTR_"$fname"_"$5".bed -b UTR_coord_"$fname"_"$5".bed -f 1 -wa > UTR_final_"$fname"_"$5".bed

echo "$fname" > counts_"$fname"_"$5".txt
echo "$5" >> counts_"$fname"_"$5".txt
awk '{print $4}' sine1andsine2_"$fname"_"$5".bed | sort | uniq | wc -l | cut -d' ' -f1 >> counts_"$fname"_"$5".txt
wc -l UTR_final_"$fname"_"$5".bed | cut -d' ' -f1 >> counts_"$fname"_"$5".txt

rm sine1_split_"$fname"_"$5".bed
rm sine1_nosplit_"$fname"_"$5".bed
rm sine1andsine2_"$fname"_"$5".bed
rm UTR_"$fname"_"$5".bed
rm UTR_final_"$fname"_"$5".bed
rm sine1_"$fname"_"$5".bam
rm sine1_"$fname"_"$5".bam.bai
rm sine2_"$fname"_"$5".bam
rm sine2_"$fname"_"$5".bam.bai
rm UTR_"$fname"_"$5".bam.bai
rm UTR_"$fname"_"$5".bam
rm UTR_coord_"$fname"_"$5".bed
rm SINE1_coord_"$fname"_"$5".bed
rm SINE2_coord_"$fname"_"$5".bed
rm sine1_bits_"$fname"_"$5".bed
rm sine2_bits_"$fname"_"$5".bed
