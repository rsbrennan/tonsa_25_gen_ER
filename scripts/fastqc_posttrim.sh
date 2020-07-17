
output=~/tonsa_genomics/analysis/fastqc/post_trim/

cd ~/tonsa_genomics/data/trimmed/

~/bin/FastQC/fastqc *fq.gz -t 4 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

multiqc ~/tonsa_genomics/analysis/fastqc/post_trim/

