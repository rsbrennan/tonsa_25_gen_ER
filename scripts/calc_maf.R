# convert afs to mafs

af <- read.table("~/tonsa_genomics/analysis/filtered_allele_freqs.txt", header=T)

rowvals <- rowMeans(af[,2:ncol(af)])

afonly <- as.matrix(af[,2:ncol(af)])

indices <- (rowvals > 0.5)

afonly[indices] <- (1- afonly)[indices,]

out <- cbind(af$SNP, as.data.frame(afonly))
colnames(out) <- c("SNP", colnames(out[2:ncol(out)]))

write.table(out, file="~/tonsa_genomics/analysis/filtered_maf_freqs.txt", sep="\t", row.names=F, quote=F)
