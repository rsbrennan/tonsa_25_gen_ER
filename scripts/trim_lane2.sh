
cd ~/tonsa_genomics/data/trimmed/lane2/

indir=/data/copepods/tonsa_rapid_genomics/lane2/
outdir=~/tonsa_genomics/data/trimmed/lane2/
for i in $(ls /data/copepods/tonsa_rapid_genomics/lane2 | grep 'fastq' | cut -f 1 -d '.' | sort | uniq)

do {

  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 2 \
    ${indir}${i}.R1.fastq.gz \
    ${indir}${i}.R2.fastq.gz \
    ${outdir}${i}.R1.qc.fq.gz s1_se \
    ${outdir}${i}.R2.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:31

echo $i done

  }

done
