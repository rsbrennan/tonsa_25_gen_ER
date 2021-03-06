cd ~/tonsa_genomics/data/aligned/merged/

samtools mpileup -Q 20 -B --max-depth 10000 --skip-indels -f /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa AA_F00_Rep1.bam AA_F00_Rep2.bam AA_F00_Rep3.bam AA_F00_Rep4.bam AA_F25_Rep1.bam AA_F25_Rep2.bam AA_F25_Rep3.bam AA_F25_Rep4.bam AH_F25_Rep1.bam AH_F25_Rep2.bam AH_F25_Rep3.bam AH_F25_Rep4.bam HA_F25_Rep1.bam HA_F25_Rep2.bam HA_F25_Rep3.bam HA_F25_Rep4.bam HH_F00_Rep1.bam HH_F00_Rep2.bam HH_F00_Rep3.bam HH_F00_Rep4.bam HH_F03_Rep1.bam HH_F03_Rep2.bam HH_F03_Rep3.bam HH_F03_Rep4.bam HH_F25_Rep1.bam HH_F25_Rep2.bam HH_F25_Rep3.bam HH_F25_Rep4.bam| java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 30 --min-reads 2 --min-avg-qual 20 --min-var-freq 0.01 --variants --p-value 0.1 > ~/tonsa_genomics/analysis/snp_all_out


# remove variants with only missing data
cat ~/tonsa_genomics/analysis/snp_all_out | wc -l
# 10,368,816

awk '(NR==1) || ($10 < 1 ) ' ~/tonsa_genomics/analysis/snp_all_out > ~/tonsa_genomics/analysis/snp_all_out.nomiss

cat ~/tonsa_genomics/analysis/snp_all_out.nomiss | wc -l
#833,910

