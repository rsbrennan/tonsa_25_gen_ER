# pi and selection.
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggbeeswarm)

#locate the directory containing the files.
dir <- "~/tonsa_genomics/analysis/pi"
files <- file.path(dir, list.files(dir))
files <- files[grep("params", files, invert=TRUE)]
files <- files[grep("F03", files, invert=TRUE)]
files <- files[grep("HH_F00", files, invert=TRUE)]
files <- files[grep("l2.pi", files, invert=FALSE)]

# read in files
d <- lapply(files, fread)

# rename without path, etc.
names(d) <- (str_split_fixed(files, "/", 5)[,5] %>% str_split_fixed( "[.]", 3))[,1]

# assign column names
d <- lapply(d, function(x) {
  colnames(x) <- c("gene", "position", "n_variants", "prop_covered", "pi")
  x
})

# convert pi to numeric
d <- lapply(d, function(x) {
  x$pi <- as.numeric(as.character(x$pi))
  x
})

# make snp name column
d <- lapply(d, function(x) {
  x$snp <- paste(x$gene, x$position, sep="_")
  x
})

# assign group id to each
for(i in 1:length(d)){
 d[[i]]$gp <- names(d)[i]

}

pops <- names(d)

for(i in 1:length(d)){
 print(sum(!is.na(d[[i]]$pi)))

}

# drop na's from the dataset
d2 <- vector(mode = "list", length = 16)

for(i in 1:length(d)){
 d2[[i]] <- d[[i]][!is.na(d[[i]]$pi),]
}


###############################################
########
######## reformat for stats
########
################################################


v1 <- vector(mode = "list", length = 16)

for(i in 1:length(d2)){
 v1[[i]] <- d2[[i]]$snp
}

overlap <- Reduce(intersect, v1)
length(overlap)
#[1] 1940
# then parse down to only the intersect
d3 <- vector(mode = "list", length = 16)
for(i in 1:length(d2)){
 d3[[i]] <- d2[[i]][d2[[i]]$snp %in% overlap,]
}

for(i in 1:length(d3)){
 print(nrow(d3[[i]]))
}

#merge entire dist.
mdat <- bind_rows(d3)
nrow(mdat)

mdat$treatment <- substr(mdat$gp,1,6)

#mdat$treatment[which(mdat$treatment == "AA_F00" | mdat$treatment == "HH_F00")] <- "founding"

#mdat$gp <- ordered(mdat$gp, 
 #                     levels = c("founding", "AA_F25", "AH_F25", "HA_F25", "HH_F25"))


forstats <- cbind(mdat$gp, as.character(mdat$treatment), mdat$pi)
colnames(forstats) <- c("sample", "treatment", "pi")

write.table(forstats,"~/tonsa_genomics/analysis/pi/forstats.txt", quote=FALSE, sep="\t", row.names=F)


###############################################
########
######## run stats
########
################################################

dat <- read.csv("~/tonsa_genomics/analysis/pi/forstats.txt", header=T, sep="\t")
dat$treatment <- as.character(dat$treatment)
dat$rep <- substr(dat$sample, 8,11)

dat$treatment[grep("AA_F00", dat$sample)] <- "AA_F00"
dat$treatment[grep("HH_F00", dat$sample)] <- "HH_F00"

out.p <- pairwise.wilcox.test(dat$pi,
                          dat$sample,
                          p.adjust.method="bonferroni")$p.value

out.p2 <- pairwise.wilcox.test(dat$pi,
                          dat$treatment,
                          p.adjust.method="bonferroni")$p.value


write.table(out.p,"~/tonsa_genomics/analysis/pi/pi_stats_all.txt", quote=FALSE, sep="\t", row.names=T)
write.table(out.p2,"~/tonsa_genomics/analysis/pi/pi_stats_trts.txt", quote=FALSE, sep="\t", row.names=T)

dat %>% group_by(treatment) %>%
    summarise(
        n = n(),
        mean = round(mean(pi, na.rm=T),4),
        sd = sd(pi, na.rm=T))

#  treatment     n   mean     sd
#  <chr>     <int>  <dbl>  <dbl>
#1 AA_F00     7760 0.0148 0.0111
#2 AA_F25     7760 0.0133 0.0111
#3 AH_F25     7760 0.0127 0.0111
#4 HA_F25     7760 0.0133 0.0110
#5 HH_F25     7760 0.0138 0.0112


dat %>% group_by(sample, rep) %>%
    summarise(
        n = n(),
        mean = round(mean(pi, na.rm=T),4),
        sd = sd(pi, na.rm=T)) %>%
    as.data.frame()
#         sample  rep    n   mean         sd
#1  AA_F00_Rep1 Rep1 1940 0.0144 0.01095262
#2  AA_F00_Rep2 Rep2 1940 0.0153 0.01127395
#3  AA_F00_Rep3 Rep3 1940 0.0149 0.01101943
#4  AA_F00_Rep4 Rep4 1940 0.0147 0.01095264
#5  AA_F25_Rep1 Rep1 1940 0.0138 0.01118990
#6  AA_F25_Rep2 Rep2 1940 0.0131 0.01120862
#7  AA_F25_Rep3 Rep3 1940 0.0135 0.01084632
#8  AA_F25_Rep4 Rep4 1940 0.0128 0.01108379
#9  AH_F25_Rep1 Rep1 1940 0.0130 0.01099934
#10 AH_F25_Rep2 Rep2 1940 0.0130 0.01087125
#11 AH_F25_Rep3 Rep3 1940 0.0121 0.01143521
#12 AH_F25_Rep4 Rep4 1940 0.0128 0.01099498
#13 HA_F25_Rep1 Rep1 1940 0.0132 0.01111625
#14 HA_F25_Rep2 Rep2 1940 0.0135 0.01091141
#15 HA_F25_Rep3 Rep3 1940 0.0131 0.01099260
#16 HA_F25_Rep4 Rep4 1940 0.0132 0.01093459
#17 HH_F25_Rep1 Rep1 1940 0.0137 0.01101892
#18 HH_F25_Rep2 Rep2 1940 0.0134 0.01109227
#19 HH_F25_Rep3 Rep3 1940 0.0143 0.01140782
#20 HH_F25_Rep4 Rep4 1940 0.0138 0.01112597


################
#
# figure
#
################

#ylim1 = boxplot.stats(mdat$pi)$stats[c(1, 5)]

pc <- ggplot(data=mdat, aes(x=gp, y=pi, fill=treatment)) +
  #geom_quasirandom(alpha=0.2, size=0.5, shape=21, show.legend=FALSE) +
  geom_violin(draw_quantiles = c(.5), show.legend = FALSE) +
#  geom_boxplot(outlier.shape=NA, show.legend = FALSE, alpha=0.2) +
 # scale_fill_manual(values=c("gray96", "gray57"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
#  coord_cartesian(ylim = ylim1*1.05)+
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Founding population","Ambient", "Acidic", 
                               "Warming", "Greenhouse"))


ggsave("~/tonsa_genomics/figures/pi_fig.pdf", pc, width = 6, height = 4)
