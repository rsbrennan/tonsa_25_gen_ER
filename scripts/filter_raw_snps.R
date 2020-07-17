
## filtering raw variants

library(stringr)
library(data.table)

dat <- read.table("~/tonsa_genomics/analysis/snp_all_out.nomiss", stringsAsFactors=FALSE, skip=1)
datnm <- read.table("~/tonsa_genomics/analysis/snp_all_out.nomiss", stringsAsFactors=FALSE, nrows=1)

pops <- c(
            "AA_F00_Rep1","AA_F00_Rep2","AA_F00_Rep3","AA_F00_Rep4",
            "AA_F25_Rep1","AA_F25_Rep2","AA_F25_Rep3","AA_F25_Rep4",
            "AH_F25_Rep1","AH_F25_Rep2","AH_F25_Rep3","AH_F25_Rep4",
            "HA_F25_Rep1","HA_F25_Rep2","HA_F25_Rep3","HA_F25_Rep4",
            "HH_F00_Rep1","HH_F00_Rep2","HH_F00_Rep3","HH_F00_Rep4",
            "HH_F03_Rep1","HH_F03_Rep2","HH_F03_Rep3","HH_F03_Rep4",
            "HH_F25_Rep1","HH_F25_Rep2","HH_F25_Rep3","HH_F25_Rep4")

colnames(dat) <- c(datnm[1,1:10], pops)

dat2 <- dat

# filter by coverage:

depthout <- as.data.frame(matrix(nrow=nrow(dat2), ncol=length(pops)))
colnames(depthout) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat2[,grep(i_pop, colnames(dat2))]
        depth <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,2])
        # sum up reads
        depthout[,grep(i_pop, colnames(depthout))] <- depth

    }

head(depthout)
colMeans(depthout)
hist(as.matrix(depthout))
quantile(as.matrix(depthout), c(0.5, 0.975, 0.99))

hi_cv <- (apply(depthout, 1, function(x) {ifelse((length(which(x > 912)) > 0), FALSE, TRUE)}))
sum(hi_cv)
#[1] 780195
dat3 <- dat2[hi_cv,]

#many of these are skewed by indels. only keep reads where depth of actual bialleleic snps > 40
# from the manual: Also, VarScan reports variants on a biallelic basis.
    #That is, for a given SNP call, the "reads1" column is the number of
    #reference-supporting reads (RD), and the "reads2" column is the number of
    #variant-supporting reads (AD).
    #There may be additional reads at that position showing other bases (SNP or indel variants).
    #If these other variants meet the calling criteria, they will be reported in
    #their own line. If not, it may look like you have "missing" reads.
# columns for each call are:
    #consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
keep <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(keep) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])
        # sum up reads
        keep[,grep(i_pop, colnames(keep))] <- (maj+ minor)

    }

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 40)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 505635
dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# [1] 505635

# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]
nrow(dat4)
# [1] 451522
dat3 <- dat4

# here calculate allele freqs
# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
af <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(af) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate af
        af[,grep(i_pop, colnames(af))] <- maj/(maj+ minor)

    }

# get rid of invariant sites
dat4 <- dat3[(which(rowSums(af) > 0)),]
dat3 <- dat4
af <- af[(which(rowSums(af) > 0)),]
nrow(dat3)
#[1] 451356

af.out <- (cbind(paste(dat3$Chrom, dat3$Position, sep=":"),af))

afct.maf <- (sapply(af,function(x)
          ifelse(x > 0.5, (1-x), x)))

# low maf cut off of < 0.05 in at least 4 groups.
## note that this corresponds to 2 reads at 40x, which seems reasonable.
low_maf <- (apply(afct.maf, 1, function(x) {ifelse((length(which(x > 0.05)) < 4), FALSE, TRUE)}))
sum(low_maf)
# [1] 166090

dat4 <- dat3[low_maf,]
nrow(dat4)
#[1] 166090

af_f <- af[low_maf,]

af.out <- (cbind(paste(dat4$Chrom, dat4$Position, sep=":"),af_f))

colnames(af.out) <- c("SNP", colnames(af_f))

###########
######
# save filtered genotypes
#####
##########

write.table(file="~/tonsa_genomics/analysis/filtered_variants.txt", dat4, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/tonsa_genomics/analysis/filtered_allele_freqs.txt", af.out, sep="\t",
              row.names=FALSE, quote=FALSE)
