# detecting mitochondrial haplotypes


This markdown file runs through the steps taken to pull out mitochondrial haplotypes from our pooled data. Much of the actual analyses (in terms of the phylogenetics) are based on Figueroa et al., 2020. 

Figueroa, N.J., Figueroa, D.F. and Hicks, D., 2020. Phylogeography of Acartia tonsa Dana, 1849 (Calanoida: Copepoda) and phylogenetic reconstruction of the genus Acartia Dana, 1846. Marine Biodiversity, 50(2), pp.1-20.


## What lineage is this reference genome?


Pulled down reference sequence from the figuera paper. 

clade X
>KC287411.1 Acartia tonsa voucher ActoME2406 cytochrome c oxidase subunit I (COI) gene, partial cds; mitochondrial
AGTAGGAATGTGATCAGGAATGGTTGGAACCGGTTTAAGAATAATTATTCGAATAGAGCTTGGCCAGGCT
GGAAGACTTATCGGGGATGATCAAATTTATAATGTAGTAGTCACCGCCCACGCTTTTATTATAATTTTTT
TCATAGTTATACCAATTTTAATTGGAGGATTTGGAAATTGACTAGTACCTTTAATACTAGGGGCTGCAGA
TATAGCTTTTCCTCGAATAAATAACATAAGATTTTGACTGTTATTACCCGCATTAATTATACTGTTATCT
AGGTCGCTAGTAGAGAGAGGGGCGGGAACAGGATGAACAGTATACCCCCCTTTATCAAGCAATATTGCCC
ATGCTGGAGCGTCAGTCGACTTTGCTATTTTTTCTCTCCATCTTGCAGGTGCAAGATCAATTTTAGGAGC
AGTAAATTTTATTTCCACAATTGGGAATCTTCGATCTTTTGGGATAGTACTTGATTTAATACCATTATTT
GCATGAGCAGTCTTAATTACAGCAGTTTTACTATTATTGTCTTTACCTGTATTAGCTGGTGCTATTACTA
TATTATGGACTGACCGAAATTTAAACTCTTCTTTTTATGAC

clade F
>JF304084.1 Acartia tonsa isolate ActoH087 cytochrome c oxidase subunit I (COI) gene, partial cds; mitochondrial
AACCGGGTTAAGAATAATTATTCGAATAGAGCTTGGCCAAGCAGGGAGATTAATTGGGGATGATCAAATT
TATAACGTTGTAGTAACAGCTCACGCTTTTATTATAATTTTTTTTATAGTTATGCCAATTTTAATTGGAG
GATTTGGAAACTGGCTAGTACCTTTAATACTAGGGGCCGCAGATATGGCTTTCCCACGAATAAATAACAT
GAGATTTTGATTATTATTACCCGCTTTAATTATATTATTGTCTAGGTCTTTAGTAGAAAGAGGAGCTGGA
ACAGGATGAACAGTTTACCCTCCCTTATCAAGTAATATTGCACATGCTGGAGCATCTGTTGACTTCGCTA
TTTTTTCTCTGCACCTCGCAGGTGCAAGATCAATTTTAGGAGCAGTAAATTTTATTTCTACAATTGGAAA
CCTTCGAGCTTTTGGAATGGTTCTTGATTTAATACCTTTGTTTGCATGAGCAGTGTTAATTACTGCAGTT
TTACTTTTATTATCGCTACCTGTATTGGCTGGGGCTATTACTATGTTGTTAACAGACCGAAATTTAAACT
CTTC



Pull down COI to align against.

```bash

/data/programs/ncbi-blast/makeblastdb -in /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa -parse_seqids -dbtype nucl -out reference_blast_db


cd /data/copepods/tonsa_genome

blastn -db /data/copepods/tonsa_genome/reference_blast_db -query ~/tonsa_genomics/analysis/mitotypes/coI_X_clade.fasta -out results.out

# this will output the match, or the subject fasta
# sseqid means Subject Seq-id
# sseq means Aligned part of subject sequence

blastn -db /data/copepods/tonsa_genome/reference_blast_db -query ~/tonsa_genomics/analysis/mitotypes/coI_X_clade.fasta -outfmt '6 sseqid sseq'

# LS052587.1    AGCAGGAATGTGGTCAGGAATGGTTGGAACCGGTTTAAGAATAATTATTCGAATAGAGCTTGGCCAGGCTGGAAGACTTATCGGGGATGATCAAATTTATAATGTAGTAGTCACCGCCCACGCTTTTATTATAATTTTTTTCATAGTTATACCAATTTTAATTGGAGGATTTGGAAATTGACTAGTACCTTTAATACTAGGGGCTGCAGATATAGCTTTTCCTCGAATAAATAACATAAGATTTTGACTGTTATTACCCGCATTAATTATACTGTTATCTAGATCGCTAGTAGAGAGAGGGGCGGGAACAGGATGAACAGTATATCCCCCTTTATCAAGCAATATTGCCCATGCTGGAGCGTCAGTCGACTTTGCTATTTTTTCTCTCCATCTTGCAGGTGCAAGATCAATTTTAGGAGCAGTAAATTTTATTTCCACAATTGGGAATCTTCGATCTTTTGGAATAGTACTTGATTTAATACCATTATTTGCATGAGCAGTCTTAATTACAGCAGTTTTACTATTATTGTCTTTACCTGTATTAGCTGGTGCTATTACTATATTGTTGACTGACCGAAATTTAAACTCTTCTTTTTATGA

```

then add this output to the fastas to be aligned, below.


### align my data to the coI fasta

```bash
my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=~/tonsa_genomics/analysis/mitotypes/coI_X_clade.fasta
my_samblstr=~/bin/samblaster/samblaster


cd ~/tonsa_genomics/data/trimmed/

for sample in `ls ~/tonsa_genomics/data/trimmed | grep '.fq.gz'  | cut -f 1 -d "."| uniq`

do

    echo "starting sample ${sample}"
    rep_1=$(ls ~/tonsa_genomics/data/trimmed | grep ${sample} | grep 'R1')
    rep_2=$(ls ~/tonsa_genomics/data/trimmed | grep ${sample} | grep 'R2')

    echo $rep_1
    echo $rep_2

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$sample\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $rep_1 $rep_2 | \
    $my_samblstr |\
    $my_samtools view -h -u -F 4 - | \
    $my_samtools sort - -O bam -o ~/tonsa_genomics/analysis/mitotypes/${sample}.coI.bam

done



```

# call snps

want to call snps, then get consensus. After consensus, could use opposite snp to see if another haplotype is present


```bash

cd ~/tonsa_genomics/analysis/mitotypes

for sample in `ls ~/tonsa_genomics/analysis/mitotypes | grep 'coI.bam' | cut -f 1 -d "."`

do

samtools mpileup -Q 20 -B --max-depth 14000 --skip-indels -f coI_X_clade.fasta ${sample}.coI.bam |\
java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 30 --min-reads 5 --min-avg-qual 20 --min-var-freq 0.05 --variants --p-value 0.1 --output-vcf 1 --strand-filter 0 | bgzip > ~/tonsa_genomics/analysis/mitotypes/${sample}.varscan.vcf.gz

done

# renames vcf files

mv AA_F00_Rep1.varscan.vcf.gz Ambient_F00_Rep1.varscan.vcf.gz 
mv AA_F00_Rep2.varscan.vcf.gz Ambient_F00_Rep2.varscan.vcf.gz 
mv AA_F00_Rep3.varscan.vcf.gz Ambient_F00_Rep3.varscan.vcf.gz 
mv AA_F00_Rep4.varscan.vcf.gz Ambient_F00_Rep4.varscan.vcf.gz 
mv AA_F25_Rep1.varscan.vcf.gz Ambient_F25_Rep1.varscan.vcf.gz 
mv AA_F25_Rep2.varscan.vcf.gz Ambient_F25_Rep2.varscan.vcf.gz 
mv AA_F25_Rep3.varscan.vcf.gz Ambient_F25_Rep3.varscan.vcf.gz 
mv AA_F25_Rep4.varscan.vcf.gz Ambient_F25_Rep4.varscan.vcf.gz 
mv AH_F25_Rep1.varscan.vcf.gz Acidification_F25_Rep1.varscan.vcf.gz 
mv AH_F25_Rep2.varscan.vcf.gz Acidification_F25_Rep2.varscan.vcf.gz 
mv AH_F25_Rep3.varscan.vcf.gz Acidification_F25_Rep3.varscan.vcf.gz 
mv AH_F25_Rep4.varscan.vcf.gz Acidification_F25_Rep4.varscan.vcf.gz 
mv HA_F25_Rep1.varscan.vcf.gz Warming_F25_Rep1.varscan.vcf.gz 
mv HA_F25_Rep2.varscan.vcf.gz Warming_F25_Rep2.varscan.vcf.gz 
mv HA_F25_Rep3.varscan.vcf.gz Warming_F25_Rep3.varscan.vcf.gz 
mv HA_F25_Rep4.varscan.vcf.gz Warming_F25_Rep4.varscan.vcf.gz 
mv HH_F03_Rep1.varscan.vcf.gz Greenhouse_F03_Rep1.varscan.vcf.gz 
mv HH_F03_Rep2.varscan.vcf.gz Greenhouse_F03_Rep2.varscan.vcf.gz 
mv HH_F03_Rep3.varscan.vcf.gz Greenhouse_F03_Rep3.varscan.vcf.gz 
mv HH_F03_Rep4.varscan.vcf.gz Greenhouse_F03_Rep4.varscan.vcf.gz 
mv HH_F25_Rep1.varscan.vcf.gz Greenhouse_F25_Rep1.varscan.vcf.gz 
mv HH_F25_Rep2.varscan.vcf.gz Greenhouse_F25_Rep2.varscan.vcf.gz 
mv HH_F25_Rep3.varscan.vcf.gz Greenhouse_F25_Rep3.varscan.vcf.gz 
mv HH_F25_Rep4.varscan.vcf.gz Greenhouse_F25_Rep4.varscan.vcf.gz 


```


## create consensus fasta

from bcftools: Create consensus sequence by applying VCF variants to a reference fasta file. By default, the program will apply all ALT variants. Using the --sample (and, optionally, --haplotype) option will apply genotype (or haplotype) calls from FORMAT/GT. The program ignores allelic depth


So to get the snps, we can set any alternate allele as homozygous alt (1/1). and can remove anything we want as the reference. The alternate alleles will be called in the consensus for any snp present.

Below, we find that the allele frequencies are almost exclusively 0 or 1 within each sample (or very close). So we can set these as fixed ref or fixed alternate, respectively, to determine the dominant haplotype. If we flip these, we can get the minor haplotype, if present.


```r

library(data.table)
setwd("~/tonsa_genomics/analysis/mitotypes")
files <- list.files(path="~/tonsa_genomics/analysis/mitotypes",pattern="varscan.vcf.gz")

sp_nm <- gsub(".varscan.vcf.gz", "", files)

dl = lapply(files, function(x)fread(x, header=T)) 
names(dl) <- sp_nm

# lets set the genotypes to get the major allele as the genotype.
# but if it is called snp, but at low freq, lets flag it to remove.
## we don't want these very low freq sites to be called at alt alleles in the consensus,
    #

flip_f = function(x){
    if(x > 50){
        geno <- "1/1"
    }
    else{
        geno <- "REMOVE_GENO"
    }
    return(geno)
}

# polarize for major AF to get consensus.
pdf("~/tonsa_genomics/figures/haplotype_X_af.pdf", h=9, w=11)
par(mfrow=c(5,5)) # rows, cols

for(lst_idx in 1:length(dl)){

    # get allele frequency.
    d2 <- as.data.frame(do.call('rbind', strsplit(as.character(dl[[lst_idx]]$Sample1),':',fixed=TRUE)))
    freqs <- as.numeric(gsub("%", "", d2$V7))
    hist(freqs, main=names(dl)[lst_idx], col="gray50",
        xlab="allele frequency", xlim=c(0,100), breaks=seq(0,100,2.5))
    # get the major allele genotype
    d2$V1 <- unlist(lapply(freqs, flip_f))

    dat2 <- dl[[lst_idx]]
    dat2$Sample1 <- apply( d2 , 1 , paste , collapse = ":" )

    # remove rows with the low freq variants:
    dat.out <- dat2[which(d2$V1 != "REMOVE_GENO"),]

    write.table(file=paste("~/tonsa_genomics/analysis/mitotypes/",names(dl)[lst_idx],".major.vcf", sep=""),
            dat.out,
            sep="\t", quote=F, row.names=F, col.names=F)

}

dev.off()

########
# get the minor haplotype:


flip_minor = function(x){
    if(x < 50){
        geno <- "1/1"
    }
    else{
        geno <- "REMOVE_GENO"
    }
    return(geno)
}

for(lst_idx in 1:length(dl)){

    # get allele frequency.
    d2 <- as.data.frame(do.call('rbind', strsplit(as.character(dl[[lst_idx]]$Sample1),':',fixed=TRUE)))
    freqs <- as.numeric(gsub("%", "", d2$V7))
    # get the minor allele genotype
    d2$V1 <- unlist(lapply(freqs, flip_minor))

    dat2 <- dl[[lst_idx]]
    dat2$Sample1 <- apply( d2 , 1 , paste , collapse = ":" )

    # remove rows with the low freq variants:
    dat.out <- dat2[which(d2$V1 != "REMOVE_GENO"),]

    write.table(file=paste("~/tonsa_genomics/analysis/mitotypes/",names(dl)[lst_idx],".minor.vcf", sep=""),
            dat.out,
            sep="\t", quote=F, row.names=F, col.names=F)

}

```


Looking at the histogram. For F0 samples, there are a bunch of low frequency variants. Our reference here is Clade X. This means that the majority of the reads match clade X reference.

By F25, we get very low frequency of fixed alternate alleles (6 or so per sample), but the low frequency variants drop out. This suggests that at the beginning of the experiment, we had mixed lines. Majority clade X, but a minority of another clade. This minor mitotype is lost in all lines by F25, and even by F3. 

We can get the consensus for the minor and major haplotypes below and check to see what their haplotypes actually are.



### get consensus

```bash

cd ~/tonsa_genomics/analysis/mitotypes


export BCFTOOLS_PLUGINS=~/bin/bcftools/plugins

for sample in `ls ~/tonsa_genomics/analysis/mitotypes | grep 'major.vcf$' | cut -f 1 -d "."`

do

# add header to vcf:
    zcat ${sample}.varscan.vcf.gz | grep '^#' | cat -  ~/tonsa_genomics/analysis/mitotypes/${sample}.major.vcf | bgzip >~/tonsa_genomics/analysis/mitotypes/${sample}.major.vcf.gz

#index
    tabix -p vcf -f ~/tonsa_genomics/analysis/mitotypes/${sample}.major.vcf.gz


    cat coI_X_clade.fasta | /data/programs/bcftools-1.9/bcftools consensus --sample Sample1 ~/tonsa_genomics/analysis/mitotypes/${sample}.major.vcf.gz > consensus.major.${sample}.fa

# make name in consensus fa more informative.

    sed -i "s/^>.*$/>${sample}_major/g" consensus.major.${sample}.fa

done


### minor haplotype:

for sample in `ls ~/tonsa_genomics/analysis/mitotypes | grep 'minor.vcf$' | grep 'AA' | cut -f 1 -d "."`

do

# add header to vcf:
    zcat ${sample}.varscan.vcf.gz | grep '^#' | cat -  ~/tonsa_genomics/analysis/mitotypes/${sample}.minor.vcf | bgzip >~/tonsa_genomics/analysis/mitotypes/${sample}.minor.vcf.gz

#index
    tabix -p vcf -f ~/tonsa_genomics/analysis/mitotypes/${sample}.minor.vcf.gz

    cat coI_X_clade.fasta | /data/programs/bcftools-1.9/bcftools consensus --sample Sample1 ~/tonsa_genomics/analysis/mitotypes/${sample}.minor.vcf.gz > consensus.minor.${sample}.fa
# make name in consensus fa more informative.
    sed -i "s/^>.*$/>${sample}_minor/g" consensus.minor.${sample}.fa

done

#cat coI.fasta | /data/programs/bcftools-1.9/bcftools consensus --sample AA_F00_Rep1 AA_F00_Rep1.alt.vcf.gz > consensus.alt.fa

```

### run mr bayes

merge fasta files:

downloaded fasta files from NCBI: searched using list from supp of figuerra. Then added A. dana and my samples to this.

Aligned files using muscle.

raw fasta is: figueroa_aligned.fasta

`./mitotype_pipeline.sh` aligns the reads, makes input for mr bayes, and outputs a rough plot.

see help associated with this script for details. `./mitotype_pipeline.sh --help`


```bash

cd ~/tonsa_genomics/analysis/mitotypes

cat consensus*.fa  > /data/copepods/mitotypes/EE_tonsa.fasta

# then run mitotyping pipeline:

cd /data/copepods/mitotypes

./mitotype_pipeline.sh EE_tonsa.fasta EE_tonsa_multigen

```

