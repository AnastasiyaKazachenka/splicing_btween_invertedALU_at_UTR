#!/bin/sh

#SBATCH --job-name=ALU_splicing
#SBATCH --ntasks=8 #number of cores
#SBATCH --nodes=1 # Ensure that all cores are on one machine 
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4096 #memory per core in MB



fullpath=$1   #collects string of input file from 1st argument (*.fastq)
fastq="${fullpath##*/}"
fname="${fastq%.bam}"

indirectory_full="${fullpath%/*}"
indirectory="${indirectory_full##*/}"

module load SAMtools/1.8-foss-2016b
module load BEDTools/2.26.0-foss-2016b


while read x y z n
do
    ../../shellscripts/peptise_SINE_work/inverted_ALU_splicing_onepair.sh $1 $x $y $z $n
done < UTRs_with_ALUinvertedpairs_short_plus10and10bp.csv

paste ./"$indirectory"/counts_"$fname"_*.txt > ./"$indirectory"/all_inverted_ALU_"$fname".txt
rm ./"$indirectory"/counts_"$fname"_*.txt