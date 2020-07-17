my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=/data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fna.gz
my_samblstr=~/bin/samblaster/samblaster

cd ~/tonsa_genomics/data/trimmed/lane2/

for sample in `ls ~/tonsa_genomics/data/trimmed/lane2 | grep '.fq.gz' | cut -f 1 -d "."| uniq`

do

    echo "starting sample ${sample}"
    #starting with only name of rep. need to pull out files

    rep_1=$(ls ~/tonsa_genomics/data/trimmed | grep ${sample} | grep 'R1')
    rep_2=$(ls ~/tonsa_genomics/data/trimmed | grep ${sample} | grep 'R2')

    echo $rep_1
    echo $rep_2

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$sample\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $rep_1 $rep_2 | \
    $my_samblstr |\
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o ~/tonsa_genomics/data/aligned/lane2/${sample}_lane2.bam

done
