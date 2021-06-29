
library(data.table)
library(scales)
library(ggplot2)

aaCmh <- fread('~/tonsa_genomics/analysis/AA.cmh',header=FALSE)
hhCmh <- fread('~/tonsa_genomics/analysis/HH.cmh',header=FALSE)
ahCmh <- fread('~/tonsa_genomics/analysis/AH.cmh',header=FALSE)
haCmh <- fread('~/tonsa_genomics/analysis/HA.cmh',header=FALSE)

d1 <- data.frame(CHR = aaCmh$V1, POS = aaCmh$V2, SNP =paste(aaCmh$V1, aaCmh$V2, sep=":"),
                  aa = aaCmh$V32, ah = ahCmh$V32, ha = haCmh$V32, hh = hhCmh$V32)


write.table(d1, "~/tonsa_genomics/analysis/cmh_results.txt", sep="\t", quote=FALSE, row.names=FALSE)
