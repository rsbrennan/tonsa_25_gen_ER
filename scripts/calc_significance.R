
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

# plot num sig, num lab removed, etc.
num_sig <- data.frame(

            group=c("AA", "AH", "HA", "HH"),
            total_sig = c(table(df.out$aa_sig)[2], table(df.out$ah_sig)[2],table(df.out$ha_sig)[2],table(df.out$hh_sig)[2]),
            lab_removed = c(NA, table(df.out$ah_sig_nolab)[2], table(df.out$ha_sig_nolab)[2], table(df.out$hh_sig_nolab)[2] )            )

num_sig$percent_sig_all <- round(num_sig$total_sig/nrow(df.out), 4)*100
num_sig$percent_sig_lab_rm <- round(num_sig$lab_removed/nrow(df.out), 4)*100
num_sig$percent_lab_rm <- 100-(round(num_sig$lab_removed/num_sig$total_sig, 4)*100)
num_sig$num_lab_rm <- num_sig$total_sig - num_sig$lab_removed

num_sig_melt <- data.frame(

            group=rep(c("AA", "AH", "HA", "HH"),2),
            sig_count = c(table(df.out$aa_sig)[2], table(df.out$ah_sig)[2],table(df.out$ha_sig)[2],table(df.out$hh_sig)[2],
                      c(NA, table(df.out$ah_sig_nolab)[2], table(df.out$ha_sig_nolab)[2], table(df.out$hh_sig_nolab)[2] )),
            sig_group = c(rep("total", 4), rep("lab removed", 4)),
            percent = c(num_sig$percent_sig_all, num_sig$percent_sig_lab_rm)
            )

num_sig_melt$sig_group <- factor(num_sig_melt$sig_group, levels=c("total", "lab removed"))

num_sig_plot <- num_sig[2:4,]

labs <- c("Acidic", "Warm", "Green-\nhouse")

p1 <- ggplot(num_sig_melt, aes(x=sig_group, y=sig_count, group=group, fill=group, shape=group, color=group)) +
  geom_line(show.legend=F)  +
  geom_point(size=5, color="black") + theme_bw() +
  ggtitle("Number of loci significant")+
  scale_y_continuous(breaks = round(seq(0, max(num_sig_melt$sig_count, na.rm=T), by = 2500),1),limits = c(0, NA)) +
    scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidic", 
                               "Warming", "Greenhouse")) +
    scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidic", 
                               "Warming", "Greenhouse"))+
    scale_shape_manual(values=c(21,22, 23, 24)) +
    ylab("Number of significant loci") +
    xlab("") +
    guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24), 
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)+ theme(legend.title = element_blank())



p2 <- ggplot(num_sig_melt, aes(x=sig_group, y=percent, group=group, fill=group,shape=group,  color=group)) +
  geom_line()  +
  geom_point(size=5, color="black") + theme_bw() +
  ggtitle("Percent of total loci") + ylim(0, 6.2) +
    scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidic", 
                               "Warming", "Greenhouse")) +
    scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidic", 
                               "Warming", "Greenhouse"))+
    scale_shape_manual(values=c(21,22, 23, 24)) +
    ylab("Percent of loci significant") +
    xlab("") +
    guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24), 
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)+ theme(legend.title = element_blank())


p3 <- ggplot(num_sig_plot, aes(x=group, y=percent_lab_rm, group=group, fill=group,shape=group,  color=group)) +
  geom_point(size=5, color="black") + theme_bw() +
  ggtitle("Percent of loci removed\ndue to lab adaptation") +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0, NA)) +
    scale_fill_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c( "Acidic",
                               "Warming", "Greenhouse")) +
    scale_color_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c( "Acidic",
                               "Warming", "Greenhouse"))+
    scale_shape_manual(values=c(22, 23, 24)) +
    ylab("Percent of loci removed") +
    xlab("") +
    guides(fill=guide_legend(override.aes=list(
        shape=c(22, 23, 24), 
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)+ theme(legend.title = element_blank())+ scale_x_discrete(labels= labs)


p4 <- ggplot(num_sig_plot, aes(x=group, y=num_lab_rm, group=group, fill=group,shape=group,  color=group)) +
  geom_point(size=5, color="black") + theme_bw() +
  ggtitle("Number of loci removed\ndue to lab adaptation") +
  scale_y_continuous(breaks = c(0, 2500, 5000, 6000), limits = c(0, NA)) +
    scale_fill_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Acidic", 
                               "Warming", "Greenhouse")) +
    scale_color_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c( "Acidic", 
                               "Warming", "Greenhouse"))+
    scale_shape_manual(values=c(22, 23, 24)) +
    ylab("Number of loci removed") +
    xlab("") +
    guides(fill=guide_legend(override.aes=list(
        shape=c(22, 23, 24), 
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)+ theme(legend.title = element_blank()) + scale_x_discrete(labels= labs)


pall <- ggarrange(p1, p2, p3, p4, nrow=1, ncol=4, common.legend=T, legend="bottom")
ggsave("~/tonsa_genomics/figures/num_sig_loci.pdf",pall, w=13, h=5)


df.out <- read.table("~/tonsa_genomics/analysis/cmh.fdr.txt", header=T)


################################################
### overlap of different groups:
################################################

ah_lab <- rep(FALSE, nrow(df.out))
ha_lab <- rep(FALSE, nrow(df.out))
hh_lab <- rep(FALSE, nrow(df.out))

ah_lab[which(df.out$ah_fdr < 0.05 & df.out$aa_sig == TRUE)] <- TRUE
ha_lab[which(df.out$ha_fdr < 0.05 & df.out$aa_sig == TRUE)] <- TRUE
hh_lab[which(df.out$hh_fdr < 0.05 & df.out$aa_sig == TRUE)] <- TRUE

ah <- sum(ah_lab == TRUE & ha_lab == FALSE & hh_lab == FALSE)
ha <- sum(ah_lab == FALSE & ha_lab == TRUE & hh_lab == FALSE)
hh <- sum(ah_lab == FALSE & ha_lab == FALSE & hh_lab == TRUE)

ah_ha <- sum(ah_lab == TRUE & ha_lab == TRUE & hh_lab == FALSE)
ah_hh <- sum(ah_lab == TRUE & ha_lab == FALSE & hh_lab == TRUE)
hh_ha <- sum(ah_lab == FALSE & ha_lab == TRUE & hh_lab == TRUE)

all <- sum(ah_lab == TRUE & ha_lab == TRUE & hh_lab == TRUE)

library(eulerr)

lab_overlaps <- euler(c("Acidic" = ah,
                        "Warm" = ha,
                        "Greenhouse" = hh,
                "Acidic&Warm" = ah_ha,
                "Acidic&Greenhouse" = ah_hh,
                "Greenhouse&Warm" = hh_ha,
                "Acidic&Warm&Greenhouse" = all))

pdf("~/tonsa_genomics/figures/lab_adaptation_venn.pdf",
    height = 3, width = 3)
plot(lab_overlaps,
     fills = list(fill = c("#F2AD00", "#00A08A", "#CC3333")),
     #edges = list(lty = 1:3),
     labels = list(font = 2),
     quantities = TRUE)

dev.off()


##########################################################
##### adaptation venn
##########################################################

ah_lab <- rep(FALSE, nrow(df.out))
ha_lab <- rep(FALSE, nrow(df.out))
hh_lab <- rep(FALSE, nrow(df.out))

ah_lab[which(df.out$ah_fdr < 0.05 & df.out$aa_sig == FALSE)] <- TRUE
ha_lab[which(df.out$ha_fdr < 0.05 & df.out$aa_sig == FALSE)] <- TRUE
hh_lab[which(df.out$hh_fdr < 0.05 & df.out$aa_sig == FALSE)] <- TRUE

ah <- sum(ah_lab == TRUE & ha_lab == FALSE & hh_lab == FALSE)
ha <- sum(ah_lab == FALSE & ha_lab == TRUE & hh_lab == FALSE)
hh <- sum(ah_lab == FALSE & ha_lab == FALSE & hh_lab == TRUE)

ah_ha <- sum(ah_lab == TRUE & ha_lab == TRUE & hh_lab == FALSE)
ah_hh <- sum(ah_lab == TRUE & ha_lab == FALSE & hh_lab == TRUE)
hh_ha <- sum(ah_lab == FALSE & ha_lab == TRUE & hh_lab == TRUE)

all <- sum(ah_lab == TRUE & ha_lab == TRUE & hh_lab == TRUE)


lab_overlaps <- euler(c("Acidic" = ah,
                        "Warm" = ha,
                        "Greenhouse" = hh,
                "Acidic&Warm" = ah_ha,
                "Acidic&Greenhouse" = ah_hh,
                "Greenhouse&Warm" = hh_ha,
                "Acidic&Warm&Greenhouse" = all))

pdf("~/tonsa_genomics/figures/adaptation_venn.pdf",
    height = 3, width = 3)
plot(lab_overlaps,
     fills = list(fill = c("#F2AD00", "#00A08A", "#CC3333")),
     #edges = list(lty = 1:3),
     labels = list(font = 2),
     quantities = TRUE)

dev.off()


#plot p-value vs new p-value

df.out$hh_sig_nolab <- factor(df.merge$hh_sig_nolab, levels = c("TRUE", "FALSE"))

p1 <- ggplot(df.out, aes(x=-log10(aa_fdr), y=-log10(aa_pval), color=aa_sig)) +
  geom_point(alpha=0.2) + theme_bw() +
  ggtitle("AA CMH vs simulation") +
  xlab("-log10(simulation p-value)") + ylab("-log10(empirical p-value)")+ 
  labs(color='Significant') +
  guides(colour = guide_legend(override.aes = list(alpha=1)))

p2 <- ggplot(df.out, aes(x=-log10(aa_fdr), y=-log10(aa), color=aa_sig)) +
  geom_point(alpha=0.2) + theme_bw() + ylim(0,15) +
    ggtitle("AA CMH vs simulation")+
  xlab("-log10(simulation p-value)") + ylab("-log10(empirical p-value)")+ 
  labs(color='Significant') +
  guides(colour = guide_legend(override.aes = list(alpha=1)))


p3 <- ggplot(ddf.out, aes(x=-log10(hh_fdr), y=-log10(hh_pval), color=hh_sig_nolab)) +
    geom_point(alpha=0.2, data = subset(df.merge, hh_sig_nolab == 'TRUE'),
             aes(x=-log10(hh_fdr), y=-log10(hh), color = hh_sig_nolab)) +
    geom_point(alpha=0.2, data = subset(df.merge, hh_sig_nolab == 'FALSE'),
             aes(x=-log10(hh_fdr), y=-log10(hh), color = hh_sig_nolab)) +
    theme_bw()+
      ggtitle("HH CMH vs simulation, lab adaptation removed")+
  xlab("-log10(simulation p-value)") + ylab("-log10(empirical p-value)")+ 
  labs(color='Significant') +
  guides(colour = guide_legend(override.aes = list(alpha=1)))


p4 <- ggplot(df.merge, aes(x=-log10(hh_fdr), y=-log10(hh), color=hh_sig_nolab)) +
  #geom_point(alpha=1) + theme_bw() + ylim(0,15) +
    geom_point(alpha=0.2, data = subset(df.merge, hh_sig_nolab == 'TRUE'),
             aes(x=-log10(hh_fdr), y=-log10(hh), color = hh_sig_nolab)) +
    geom_point(alpha=0.2, data = subset(df.merge, hh_sig_nolab == 'FALSE'),
             aes(x=-log10(hh_fdr), y=-log10(hh), color = hh_sig_nolab)) +
    theme_bw()+ ylim(0,15) +
          ggtitle("HH CMH vs simulation, lab adaptation removed")+
  xlab("-log10(simulation p-value)") + ylab("-log10(empirical p-value)")+ 
  labs(color='Significant') +
  guides(colour = guide_legend(override.aes = list(alpha=1)))


ggsave("~/tonsa_genomics/figures/sim_vs_cmh.png",ggarrange(p1,p2,p3,p4, nrow=2, ncol=2, common.legend=T), w=10, h=10)

