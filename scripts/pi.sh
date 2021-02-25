cd ~/tonsa_genomics/analysis/pi

for i in $(ls ~/tonsa_genomics/data/aligned/merged | cut -f 1 -d '.'); do

  echo "Status: starting $i";

  samtools mpileup -Q 20 -B --max-depth 6000 --skip-indels     -f /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa     ~/tonsa_genomics/data/aligned/merged/${i}.bam     > ~/tonsa_genomics/analysis/pi/${i}.mpileup;

  echo "Status: $i pileup done; starting popoolation";

  perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 50 --min-qual 20     --min-coverage 20 --min-count 2 --max-coverage 3000 --min-covered-fraction 0.5     --input ~/tonsa_genomics/analysis/pi/${i}.mpileup     --window-size 100 --step-size 100     --output ~/tonsa_genomics/analysis/pi/${i}.l2.pi --measure pi;

  rm ~/tonsa_genomics/analysis/pi/${i}.mpileup;

  echo "Status: $i popoolation done, pileup removed";

done
