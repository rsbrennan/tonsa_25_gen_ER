
output=~/tonsa_genomics/analysis/fastqc/pre_trim/

cd /data/copepods/tonsa_rapid_genomics

~/bin/FastQC/fastqc *fastq.gz -t 4 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

multiqc ~/tonsa_genomics/analysis/fastqc/pre_trim/

