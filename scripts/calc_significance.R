
library(qvalue)
library(data.table)
library(tidyverse)

# convert afs to mafs

af <- read.table("~/tonsa_genomics/analysis/filtered_allele_freqs.txt", header=T)

rowvals <- rowMeans(af[,2:ncol(af)])

afonly <- as.matrix(af[,2:ncol(af)])

indices <- (rowvals > 0.5)

afonly[indices] <- (1- afonly)[indices,]

out <- cbind(af$SNP, as.data.frame(afonly))
colnames(out) <- c("SNP", colnames(out[2:ncol(out)]))

write.table(out, file="~/tonsa_genomics/analysis/filtered_maf_freqs.txt", sep="\t", row.names=F, quote=F)

# calc p-value thresholds for each locus. base this on wf sims from starting allele freqs, where do we see shifts larger than drift.

#library(purrr)

cmh <- read.csv("~/tonsa_genomics/analysis/cmh_results.txt", sep="\t", )

filenames <- list.files("~/tonsa_genomics/analysis/cmh_simulation", pattern="*.pval", full.names=TRUE)

dflist <- lapply(filenames[1:500],fread)

colnames <- c("SNP","pvalue") 
dflist <- lapply(dflist, setNames, colnames)

df <- reduce(dflist, full_join, by = "SNP") %>% replace(., is.na(.), 1)

#rename columns:

colnames(df)[2:ncol(df)] <- paste0("pval_", seq(1,ncol(df)-1,1))

# merge the two data frames
newdata <- -log10(df[,2:ncol(df)])
df2 <- cbind(df[,1], newdata)
alldf <- left_join(cmh, df2, by="SNP") %>% replace(., is.na(.), 0)

# make df to hold output:
df.out <- as.data.frame(matrix(ncol=11, nrow=nrow(alldf)))
colnames(df.out) <- c("CHR", "POS","SNP", "aa_pval", "ah_pval", "ha_pval", "hh_pval",
                             "aa_fdr", "ah_fdr", "ha_fdr", "hh_fdr")
df.out[1:7] <- cmh
aain <- -log10(df.out$aa_pval)
ahin <- -log10(df.out$ah_pval)
hain <- -log10(df.out$ha_pval)
hhin <- -log10(df.out$hh_pval)

# calc emp. pvals with the simulated dataset:
df.out$aa_fdr <- empPvals(stat = aain, stat0 = as.matrix(alldf[,8:ncol(alldf)]))
df.out$ah_fdr <- empPvals(stat = ahin, stat0 = as.matrix(alldf[,8:ncol(alldf)]))
df.out$ha_fdr <- empPvals(stat = hain, stat0 = as.matrix(alldf[,8:ncol(alldf)]))
df.out$hh_fdr <- empPvals(stat = hhin, stat0 = as.matrix(alldf[,8:ncol(alldf)]))

# significant or not?
df.out$aa_sig <- FALSE
df.out$ah_sig <- FALSE
df.out$ha_sig <- FALSE
df.out$hh_sig <- FALSE

df.out$aa_sig[which(df.out$aa_fdr < 0.05)] <- TRUE
df.out$ah_sig[which(df.out$ah_fdr < 0.05)] <- TRUE
df.out$ha_sig[which(df.out$ha_fdr < 0.05)] <- TRUE
df.out$hh_sig[which(df.out$hh_fdr < 0.05)] <- TRUE

df.out$ah_sig_nolab <- FALSE
df.out$ha_sig_nolab <- FALSE
df.out$hh_sig_nolab <- FALSE

df.out$ah_sig_nolab[which(df.out$ah_fdr < 0.05 & df.out$aa_sig == FALSE)] <- TRUE
df.out$ha_sig_nolab[which(df.out$ha_fdr < 0.05 & df.out$aa_sig == FALSE)] <- TRUE
df.out$hh_sig_nolab[which(df.out$hh_fdr < 0.05 & df.out$aa_sig == FALSE)] <- TRUE

table(df.out$aa_sig)
#  FALSE   TRUE
 # 382541  12126
table(df.out$ah_sig)
# FALSE   TRUE
#387151   7516
table(df.out$ha_sig)
# FALSE   TRUE
#376343  18324
table(df.out$hh_sig)
# FALSE   TRUE
#370908  23759

table(df.out$ah_sig_nolab)
# FALSE   TRUE
#391267   3400

table(df.out$ha_sig_nolab)
# FALSE   TRUE
#381358  13309
table(df.out$hh_sig_nolab)
# FALSE   TRUE
#375570  19097

write.table(file="~/tonsa_genomics/analysis/cmh.fdr.txt", df.out,
              sep="\t", quote=FALSE, row.names=FALSE)


df.out <- read.table("~/tonsa_genomics/analysis/cmh.fdr.txt", header=T)
