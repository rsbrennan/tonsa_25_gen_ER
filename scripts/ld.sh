#!/bin/bash

sort -k1,1 -k2,2n /users/r/b/rbrennan/tonsa_genomics/analysis/variants.vcf > /users/r/b/rbrennan/tonsa_genomics/analysis/variants.sort.vcf

samtools=~/bin/samtools-1.6/samtools

# generate sam files for analysis

cd ~/tonsa_genomics/analysis/ldx

$samtools view -H ~/tonsa_genomics/data/aligned/merged/AA_F00_Rep1.bam > ~/tonsa_genomics/analysis/ldx/sam.header

cat ~/tonsa_genomics/analysis/snp_all_out | tail -n +2| cut -f 1 | sort | uniq > ~/tonsa_genomics/analysis/ldx/scaffolds.txt

vcf_i=~/tonsa_genomics/analysis/variants.sort.vcf

cat $vcf_i | grep '^#' > ~/tonsa_genomics/analysis/ldx/vcf_header.txt

#######################
# loop over samples
#######################

for rep in $(ls ~/tonsa_genomics/data/aligned/merged/*.bam | cut -f 10 -d "/" | cut -f 1 -d "."); do

  echo $rep

  echo "sorting ${rep}"

  $samtools sort ~/tonsa_genomics/data/aligned/merged/${rep}.bam -o ~/tonsa_genomics/analysis/ldx/${rep}.sorted.bam
  $samtools index ~/tonsa_genomics/analysis/ldx/${rep}.sorted.bam
  echo "done sorting"

  bam_in=~/tonsa_genomics/analysis/ldx/${rep}.sorted.bam

#  # cycling over each scaffold, to be memory efficient.

  echo "begin calculating LD for every scaffold"

  while read scaf; do

            #subset to correct scaffold
            $samtools view -h $bam_in $scaf > tmp.${rep}.scaf.sam

            #subset vcf
            cat $vcf_i | grep -v '^#' | grep -w $scaf > tmp1.${rep}.vcf
            cat  ~/tonsa_genomics/analysis/ldx/vcf_header.txt tmp1.${rep}.vcf > tmp.${rep}.vcf

            #calculate ld on the scaffold
            ~/bin/LDx.pl -l 10 -h 1000 -s 500 -q 20 -a 0.05 -i 5 tmp.${rep}.scaf.sam tmp.${rep}.vcf \
              > ~/tonsa_genomics/analysis/ldx/ld.${rep}.tmp

            # add scaffold identifier

            sed -i "s/^/$scaf\t/g" ~/tonsa_genomics/analysis/ldx/ld.${rep}.tmp

            cat  ~/tonsa_genomics/analysis/ldx/ld.${rep}.tmp >> ~/tonsa_genomics/analysis/ldx/${rep}.ld.out

    done<~/tonsa_genomics/analysis/ldx/scaffolds.txt

    rm ld.${rep}.tmp
    rm tmp.${rep}.vcf
    rm tmp.${rep}.scaf.sam
#    # rm ~/tonsa_genomics/analysis/ldx/${rep}.sorted.bam
    echo $rep
    echo "done"

done
