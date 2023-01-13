#STEP 1 - get all UTRs

cat gencode.v37.basic.annotation.gtf | awk 'OFS="\t" {if ($3=="UTR") {print $1,$4,$5,$10,$6,$7,$3,$16}}' | tr -d '";' > allUTRs_gencovev37_basic.bed

#STEP 2 - get all ALUs

grep "Alu" GRCh38.78_Dfam2_full.bed > onlyALU_GRCh38.78_Dfam2.bed

#STEP 3 - all ALUs that overlap UTRs

module load BEDTools/2.26.0-foss-2016b
bedtools intersect -a ../../gencode.v37/allUTRs_gencovev37_basic.bed -b ../../GRCh38.78_Dfam2_full.gtf/onlyALU_GRCh38.78_Dfam2.bed -wo > UTR_and_overlaping_ALUs.bed

#Remove UTRs that overlap only 1 Alu

cat UTR_and_overlaping_ALUs.bed | sort | uniq > UTR_and_overlaping_ALUs.uniq.bed
awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' UTR_and_overlaping_ALUs.uniq.bed | sort | uniq -c | awk '$1>1 {print $2,$3,$4,$5,$6,$7,$8}' | sed 's/ /\t/g' | grep -f - UTR_and_overlaping_ALUs.uniq.bed > UTR_overlaping_2andMORE_ALUs.bed

#remove UTRs that dont have inverted ALUs

awk '$14=="-"' UTR_overlaping_2andMORE_ALUs.bed > UTR_overlaping_2andMORE_ALUs.minus_strand.bed
awk '$14=="+" {print $1,$2,$3,$4,$5,$6,$7,$8}' UTR_overlaping_2andMORE_ALUs.bed | sort | uniq | sed 's/ /\t/g' | grep -f - UTR_overlaping_2andMORE_ALUs.minus_strand.bed | awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' | sort | uniq | sed 's/ /\t/g' | grep -f - UTR_overlaping_2andMORE_ALUs.bed > UTR_overlaping_2andMORE_ALUs.invertedALUs.bed
rm UTR_overlaping_2andMORE_ALUs.minus_strand.bed

#make ALU pairs

awk '$14=="-" {print $9,$10,$11,$12,$13,$14}' UTR_overlaping_2andMORE_ALUs.invertedALUs.bed |  sed 's/ /\t/g' | sort | uniq > invertedALUs.minus_strand.bed
awk '$14=="+"' UTR_overlaping_2andMORE_ALUs.invertedALUs.bed > UTR_overlaping_2andMORE_ALUs.invertedALUs.plus_strand.bed

bedtools intersect -a UTR_overlaping_2andMORE_ALUs.invertedALUs.plus_strand.bed -b invertedALUs.minus_strand.bed -wa -wb > UTRs_with_ALUinvertedpairs.bed

rm UTR_overlaping_2andMORE_ALUs.invertedALUs.plus_strand.bed
rm invertedALUs.minus_strand.bed