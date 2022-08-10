
library(scales)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(eulerr)
library(UpSetR)

af <- read.csv("~/tonsa_genomics/analysis/filtered_maf_freqs.txt", sep="\t", header=T)
annot <- read.csv("~/tonsa_genomics/analysis/results_with_gene_ids.txt", header=T, sep="\t")
df.out <- read.table("~/tonsa_genomics/analysis/cmh.fdr.txt", header=T)


# note that my allele freq calcs are the same as given by the poolseq package, when the alleles are polarized by the reference.
  # this is always good to see
aa.f0 <- rowMeans(af[,grep("AA_F00", colnames(af))])
hh.f0 <- rowMeans(af[,grep("HH_F00", colnames(af))])
hh.aa.f0 <- rowMeans(af[,grep("HH_F00|AA_F00", colnames(af))])
aa.f25 <- rowMeans(af[,grep("AA_F25", colnames(af))])
ah.f25 <- rowMeans(af[,grep("AH_F25", colnames(af))])
ha.f25 <- rowMeans(af[,grep("HA_F25", colnames(af))])
hh.f25 <- rowMeans(af[,grep("HH_F25", colnames(af))])
hh.f3 <- rowMeans(af[,grep("HH_F03", colnames(af))])

aa_hh <- as.data.frame(cbind(af$SNP,aa.f0, hh.f25))
aa_hh$aa.f0 <- round(as.numeric(aa_hh$aa.f0), 4)
aa_hh$hh.f25 <- round(as.numeric(aa_hh$hh.f25), 4)

df <- data.frame(SNP=df.out$SNP,
                  aa.f0.af = aa.f0, hh.f0.af = hh.f0,
                  aa.f25.af = aa.f25, ah.f25.af = ah.f25, ha.f25.af = ha.f25, hh.f25.af = hh.f25,
                  hh.f3.af = hh.f3)

df.af <- full_join(df, df.out, by="SNP")
df.af$delta_aa <- df.af$aa.f25.af - df.af$aa.f0.af
df.af$delta_ah <- df.af$ah.f25.af - df.af$aa.f0.af
df.af$delta_ha <- df.af$ha.f25.af - df.af$aa.f0.af
df.af$delta_hh <- df.af$hh.f25.af - df.af$aa.f0.af

sum(abs(df.af$delta_hh[df.af$hh_fdr < 0.01]) < 0.05)/(sum(df.af$hh_fdr < 0.01))
sum(abs(df.af$delta_aa[df.af$aa_fdr < 0.01]) < 0.05)/(sum(df.af$aa_fdr < 0.01))


####################################################
# num sig loci, lab removed, etc.
####################################################


# plot num sig, num lab removed, etc.
num_sig <- data.frame(

            group=c("AA", "AH", "HA", "HH"),
            total_sig = c(table(df.out$aa_sig)[2], table(df.out$ah_sig)[2],table(df.out$ha_sig)[2],table(df.out$hh_sig)[2]),
            lab_removed = c(NA, table(df.out$ah_sig_nolab)[2], table(df.out$ha_sig_nolab)[2], table(df.out$hh_sig_nolab)[2] )            )

num_sig$percent_sig_all <- round(num_sig$total_sig/nrow(df.out), 4)*100
num_sig$percent_sig_lab_rm <- round(num_sig$lab_removed/nrow(df.out), 4)*100
num_sig$percent_lab_rm <- 100-(round(num_sig$lab_removed/num_sig$total_sig, 4)*100)
num_sig$num_lab_rm <- num_sig$total_sig - num_sig$lab_removed


write.table(file="~/tonsa_genomics/analysis/cmh.sigNums.txt", num_sig,
              sep="\t", quote=FALSE, row.names=FALSE)

num_sig_melt <- data.frame(

            group=rep(c("AA", "AH", "HA", "HH"),2),
            sig_count = c(table(df.out$aa_sig)[2], table(df.out$ah_sig)[2],table(df.out$ha_sig)[2],table(df.out$hh_sig)[2],
                      c(NA, table(df.out$ah_sig_nolab)[2], table(df.out$ha_sig_nolab)[2], table(df.out$hh_sig_nolab)[2] )),
            sig_group = c(rep("total", 4), rep("lab removed", 4)),
            percent = c(num_sig$percent_sig_all, num_sig$percent_sig_lab_rm)
            )

num_sig_melt$sig_group <- factor(num_sig_melt$sig_group, levels=c("total", "lab removed"))

num_sig_plot <- num_sig[2:4,]

labs <- c("Acidification", "Warming", "OWA")

p1 <- ggplot(num_sig_melt, aes(x=sig_group, y=sig_count, group=group, fill=group, shape=group, color=group)) +
  geom_line(show.legend=F)  +
  geom_point(size=5, color="black") + theme_bw() +
  ggtitle("Number of loci significant")+
  scale_y_continuous(breaks = round(seq(0, max(num_sig_melt$sig_count, na.rm=T), by = 2500),1),limits = c(0, NA)) +
    scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "OWA")) +
    scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "OWA"))+
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
  ggtitle("Percent of total loci") + 
  ylim(0, 2.7) +
    scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "OWA")) +
    scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "OWA"))+
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
                    labels = c( "Acidification",
                               "Warming", "OWA")) +
    scale_color_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c( "Acidification",
                               "Warming", "OWA"))+
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
                    labels = c("Acidification",
                               "Warming", "OWA")) +
    scale_color_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c( "Acidification",
                               "Warming", "OWA"))+
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


## prop sig alone.
num_sig_melt_aarm <- num_sig_melt[grep("AA",num_sig_melt$group, invert =T),]

p2 <- ggplot(num_sig_melt_aarm, aes(x=sig_group, y=percent, group=group, fill=group,shape=group,  color=group)) +
  geom_line()  +
  geom_point(size=5, color="black") + theme_bw() +
  ggtitle("Percent of total loci") + ylim(0, 2.7) +
    scale_fill_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Acidification",
                               "Warming", "OWA")) +
    scale_color_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Acidification",
                               "Warming", "OWA"))+
    scale_shape_manual(values=c(22, 23, 24)) +
    ylab("Percent of loci significant") +
    xlab("") +
    guides(fill=guide_legend(override.aes=list(
        shape=c(22, 23, 24),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)+ theme(legend.title = element_blank())


ggsave("~/tonsa_genomics/figures/num_sig_loci_only.pdf",p2,w=4, h=3.4)


################################################################################################
### overlap of different groups:
################################################################################################

aa_sig <- rep(FALSE, nrow(df.out))
ah_sig <- rep(FALSE, nrow(df.out))
ha_sig <- rep(FALSE, nrow(df.out))
hh_sig <- rep(FALSE, nrow(df.out))

aa_sig[which(df.out$aa_fdr < 0.01)] <- TRUE
ah_sig[which(df.out$ah_fdr < 0.01)] <- TRUE
ha_sig[which(df.out$ha_fdr < 0.01)] <- TRUE
hh_sig[which(df.out$hh_fdr < 0.01)] <- TRUE

# find unique significant
aa <- sum(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == FALSE & hh_sig == FALSE)
ah <- sum(aa_sig == FALSE  &ah_sig == TRUE  & ha_sig == FALSE & hh_sig == FALSE)
ha <- sum(aa_sig == FALSE  &ah_sig == FALSE & ha_sig == TRUE  & hh_sig == FALSE)
hh <- sum(aa_sig == FALSE  &ah_sig == FALSE & ha_sig == FALSE & hh_sig == TRUE)

# pairwise overlaps

aa_ah <- sum(aa_sig == TRUE  & ah_sig == TRUE  & ha_sig == FALSE & hh_sig == FALSE)
aa_ha <- sum(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == TRUE & hh_sig == FALSE)
aa_hh <- sum(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == FALSE & hh_sig == TRUE)
ah_ha <- sum(aa_sig == FALSE  & ah_sig == TRUE  & ha_sig == TRUE & hh_sig == FALSE)
ah_hh <- sum(aa_sig == FALSE  & ah_sig == TRUE  & ha_sig == FALSE & hh_sig == TRUE)
ha_hh <- sum(aa_sig == FALSE  & ah_sig == FALSE  & ha_sig == TRUE & hh_sig == TRUE)

# 3-way overlaps
aa_ah_hh <- sum(aa_sig == TRUE  & ah_sig == TRUE  & ha_sig == FALSE & hh_sig == TRUE)
aa_ha_hh <- sum(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == TRUE & hh_sig == TRUE)
aa_ah_ha <- sum(aa_sig == TRUE  & ah_sig == TRUE  & ha_sig == TRUE & hh_sig == FALSE)
ah_ha_hh <- sum(aa_sig == FALSE  & ah_sig == TRUE  & ha_sig == TRUE & hh_sig == TRUE)


all <- sum(aa_sig == TRUE & ah_sig == TRUE & ha_sig == TRUE & hh_sig == TRUE)

all_overlaps <- euler(c(
                        "Ambient" = aa,
                        "Acidic" = ah,
                        "Warm" = ha,
                        "OWA" = hh,
                "Ambient&Acidic" = aa_ah,
                "Ambient&Warm" = aa_ha,
                "Ambient&OWA" = aa_hh,
                "Acidic&Warm" = ah_ha,
                "Acidic&OWA" = ah_hh,
                "Warm&OWA" = ha_hh,
                "Ambient&Acidic&OWA" = aa_ah_hh,
                "Ambient&Warm&OWA" = aa_ha_hh,
                "Ambient&Acidic&Warm" = aa_ah_ha,
                "Acidic&Warm&OWA" = ah_ha_hh,
                "Ambient&Acidic&Warm&OWA" = all))

pdf("~/tonsa_genomics/figures/all_adaptation_venn.pdf",
    height = 3, width = 3)
plot(all_overlaps,
     fills = list(fill = c("#F2AD00", "#00A08A", "#CC3333")),
     #edges = list(lty = 1:3),
     labels = list(font = 2),
     quantities = TRUE)

dev.off()



################
# different plot

all_input <- all_overlaps <- c(
                        "Ambient" = aa,
                        "Acidic" = ah,
                        "Warm" = ha,
                        "OWA" = hh,
                "Ambient&Acidic" = aa_ah,
                "Ambient&Warm" = aa_ha,
                "Ambient&OWA" = aa_hh,
                "Acidic&Warm" = ah_ha,
                "Acidic&OWA" = ah_hh,
                "Warm&OWA" = ha_hh,
                "Ambient&Acidic&OWA" = aa_ah_hh,
                "Ambient&Warm&OWA" = aa_ha_hh,
                "Ambient&Acidic&Warm" = aa_ah_ha,
                "Acidic&Warm&OWA" = ah_ha_hh,
                "Ambient&Acidic&Warm&OWA" = all)


pdf("~/tonsa_genomics/figures/all_upSet.pdf",
    height = 4, width = 5)

upset(fromExpression(all_input),
        #order.by = "degree", 
        #group.by = "sets",
        keep.order = TRUE, empty.intersections = "on",
        sets = c("OWA","Warm", "Acidic","Ambient"),
        mainbar.y.label = "Loci intersecting",
        sets.x.label = "Number of significant loci",
        point.size = 3.4, line.size = 1.2, 
          sets.bar.color=rev(c('#6699CC',"#F2AD00", "#00A08A", "#CC3333")),
           text.scale = c(1.1, 1.1, 1, 1, 1.5, 1)
)

dev.off()




# lab adaptation signal

ah_lab <- rep(FALSE, nrow(df.out))
ha_lab <- rep(FALSE, nrow(df.out))
hh_lab <- rep(FALSE, nrow(df.out))

ah_lab[which(df.out$ah_fdr < 0.01 & df.out$aa_sig >= 0.01)] <- TRUE
ha_lab[which(df.out$ha_fdr < 0.01 & df.out$aa_sig >= 0.01)] <- TRUE
hh_lab[which(df.out$hh_fdr < 0.01 & df.out$aa_sig >= 0.01)] <- TRUE

ah <- sum(ah_lab == TRUE & ha_lab == FALSE & hh_lab == FALSE)
ha <- sum(ah_lab == FALSE & ha_lab == TRUE & hh_lab == FALSE)
hh <- sum(ah_lab == FALSE & ha_lab == FALSE & hh_lab == TRUE)

ah_ha <- sum(ah_lab == TRUE & ha_lab == TRUE & hh_lab == FALSE)
ah_hh <- sum(ah_lab == TRUE & ha_lab == FALSE & hh_lab == TRUE)
hh_ha <- sum(ah_lab == FALSE & ha_lab == TRUE & hh_lab == TRUE)

all <- sum(ah_lab == TRUE & ha_lab == TRUE & hh_lab == TRUE)

lab_overlaps <- euler(c("Acidic" = ah,
                        "Warm" = ha,
                        "OWA" = hh,
                "Acidic&Warm" = ah_ha,
                "Acidic&OWA" = ah_hh,
                "OWA&Warm" = hh_ha,
                "Acidic&Warm&OWA" = all))


pdf("~/tonsa_genomics/figures/lab_adaptation_venn.pdf",
    height = 3, width = 3)
plot(lab_overlaps,
     fills = list(fill = c("#F2AD00", "#00A08A", "#CC3333")),
     #edges = list(lty = 1:3),
     labels = list(font = 2),
     quantities = TRUE)

dev.off()


##########################################################
##### Remove lab adaptation- only adaptive loci venn
##########################################################

ah_lab <- rep(FALSE, nrow(df.out))
ha_lab <- rep(FALSE, nrow(df.out))
hh_lab <- rep(FALSE, nrow(df.out))

ah_lab[which(df.out$ah_fdr < 0.01 & df.out$aa_fdr >= 0.01)] <- TRUE
ha_lab[which(df.out$ha_fdr < 0.01 & df.out$aa_fdr >= 0.01)] <- TRUE
hh_lab[which(df.out$hh_fdr < 0.01 & df.out$aa_fdr >= 0.01)] <- TRUE

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



##############################################################
##############################################################
# stats:
##############################################################
##############################################################


sum(aa_sig)
sum(ah_sig)
sum(ha_sig)
sum(hh_sig)


sum(aa_sig)/nrow(df.out)
sum(ah_sig)/nrow(df.out)
sum(ha_sig)/nrow(df.out)
sum(hh_sig)/nrow(df.out)


aa <- df.out$SNP[(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == FALSE & hh_sig == FALSE)]
ah <- df.out$SNP[(aa_sig == FALSE  &ah_sig == TRUE  & ha_sig == FALSE & hh_sig == FALSE)]
ha <- df.out$SNP[(aa_sig == FALSE  &ah_sig == FALSE & ha_sig == TRUE  & hh_sig == FALSE)]
hh <- df.out$SNP[(aa_sig == FALSE  &ah_sig == FALSE & ha_sig == FALSE & hh_sig == TRUE)]

length(aa)
length(ah)
length(ha)
length(hh)

length(aa)/sum(aa_sig)
length(ah)/sum(ah_sig)
length(ha)/sum(ha_sig)
length(hh)/sum(hh_sig)


# find unique significant
aa <- df.out$SNP[(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == FALSE & hh_sig == FALSE)]
ah <- df.out$SNP[(aa_sig == FALSE  &ah_sig == TRUE  & ha_sig == FALSE & hh_sig == FALSE)]
ha <- df.out$SNP[(aa_sig == FALSE  &ah_sig == FALSE & ha_sig == TRUE  & hh_sig == FALSE)]
hh <- df.out$SNP[(aa_sig == FALSE  &ah_sig == FALSE & ha_sig == FALSE & hh_sig == TRUE)]

# pairwise overlaps

aa_ah <- df.out$SNP[(aa_sig == TRUE  & ah_sig == TRUE  & ha_sig == FALSE & hh_sig == FALSE)]
aa_ha <- df.out$SNP[(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == TRUE & hh_sig == FALSE)]
aa_hh <- df.out$SNP[(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == FALSE & hh_sig == TRUE)]
ah_ha <- df.out$SNP[(aa_sig == FALSE  & ah_sig == TRUE  & ha_sig == TRUE & hh_sig == FALSE)]
ah_hh <- df.out$SNP[(aa_sig == FALSE  & ah_sig == TRUE  & ha_sig == FALSE & hh_sig == TRUE)]
ha_hh <- df.out$SNP[(aa_sig == FALSE  & ah_sig == FALSE  & ha_sig == TRUE & hh_sig == TRUE)]

# 3-way overlaps
aa_ah_hh <- df.out$SNP[(aa_sig == TRUE  & ah_sig == TRUE  & ha_sig == FALSE & hh_sig == TRUE)]
aa_ha_hh <- df.out$SNP[(aa_sig == TRUE  & ah_sig == FALSE  & ha_sig == TRUE & hh_sig == TRUE)]
aa_ah_ha <- df.out$SNP[(aa_sig == TRUE  & ah_sig == TRUE  & ha_sig == TRUE & hh_sig == FALSE)]
ah_ha_hh <- df.out$SNP[(aa_sig == FALSE  & ah_sig == TRUE  & ha_sig == TRUE & hh_sig == TRUE)]

all <- length(which(df.out$aa_fdr < 0.01 & df.out$ah_fdr < 0.01 & df.out$ha_fdr < 0.01 & df.out$hh_fdr < 0.01))

library("SuperExactTest")
## make list with snp id for each set.

total=nrow(df.out)
#  calculate the expected overlap size:
length.gene.sets=c(sum(aa_sig), sum(ah_sig))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

#common.genes=length(aa_ah)
num.observed.overlap/(sum(aa_sig) + sum(ah_sig))
(num.observed.overlap=length(aa_ah))
(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)


### OWA and warming
(length.gene.sets=c(sum(hh_sig), sum(ha_sig)))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(num.observed.overlap=(length(ha_hh) + length(ah_ha_hh) + length(aa_ha_hh) + all))

num.observed.overlap/(sum(hh_sig) + sum(ha_sig))
(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

### owa acidification
(length.gene.sets=c(sum(hh_sig), sum(ah_sig)))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(num.observed.overlap=(length(ah_hh) + length(ah_ha_hh) + length(aa_ah_hh) + all))
num.observed.overlap/(sum(hh_sig) + sum(ah_sig))
(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

### warming acidification
(length.gene.sets=c(sum(ha_sig), sum(ah_sig)))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(num.observed.overlap=(length(ah_ha) + length(ah_ha_hh) + length(aa_ah_ha) + all))
num.observed.overlap/(sum(ha_sig) + sum(ah_sig))
(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

#################
# remove lab
#################

ah_lab <- rep(FALSE, nrow(df.out))
ha_lab <- rep(FALSE, nrow(df.out))
hh_lab <- rep(FALSE, nrow(df.out))

ah_lab[which(df.out$ah_fdr < 0.01 & df.out$aa_fdr >= 0.01)] <- TRUE
ha_lab[which(df.out$ha_fdr < 0.01 & df.out$aa_fdr >= 0.01)] <- TRUE
hh_lab[which(df.out$hh_fdr < 0.01 & df.out$aa_fdr >= 0.01)] <- TRUE


## OWA warming
(length.gene.sets=c(sum(hh_lab), sum(ha_lab)))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(num.observed.overlap=(length(ha_hh) + length(ah_ha_hh)))

num.observed.overlap/(sum(hh_lab) + sum(ha_lab))
(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

### owa acidification
(length.gene.sets=c(sum(hh_lab), sum(ah_lab)))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(num.observed.overlap=(length(ah_hh) + length(ah_ha_hh)))
num.observed.overlap/(sum(hh_lab) + sum(ah_lab))
(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

### warming acidification
(length.gene.sets=c(sum(ha_lab), sum(ah_lab)))
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(num.observed.overlap=(length(ah_ha) + length(ah_ha_hh)))
num.observed.overlap/(sum(ha_lab) + sum(ah_lab))

(FE=num.observed.overlap/num.expcted.overlap)
# The probability density of the observed intersection size is therefore:
dpsets(num.observed.overlap, length.gene.sets, n=total)
# The probability of observing intersection genes can be calculated using the cumulative probability function cpsets:
cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)



######################################################################################################
#
# af change in bins
#
######################################################################################################


#### for fdr va;s first:
# set up cut-off values
breaks <- c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,150)
#breaks <- seq(0, 2.75, 0.25)
# specify interval/bin labels
tags <- c("[0-2)","[2-4)", "[4-6)", "[6-8)", "[8-10)", "[10-12)","[12-14)", "[14-16)","[16-18)", "[18-20)", "[20-22)", "[22-24)", "[24-26)", "[26-28)", "[28-30)", "[30-150)")
plotlabs <- c("0-2","2-4", "4-6", "6-8", "8-10", "10-12","12-14", "14-16","16-18", "18-20", "20-22", "22-24", "24-26", "26-28", "28-30", "30-150")
#tags <- c("[0,0.25)","[0.25,0.5)","[0.5,0.75)","[0.75,1)",
  "[1,1.25)","[1.25,1.5)","[1.5,1.75)","[1.75,2)", "[2,2.25)","[2.25,2.5)","[2.5,2.75]")

#plotlabs <- c("0-0.25","0.25-0.5", "0.5-0.75", "0.75-1", "1-1.25", "1.25-1.5","1.5-1.75", "1.75-2",
#                "2-2.25", "2.25-2.5", "2.5-2.75")

# bucketing values into bins
df.af$hh.f25.bin <- cut(-log10(df.af$hh_fdr),
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE)
# inspect bins

n_fun <- function(x){
    return(data.frame(y = 0.02+ max(x),
                      label = round(100*(length(x)/nrow(dfplot)),2)))
}


dfplot <- df.af[!is.na(df.af$hh.f25.bin),]

library(ggbeeswarm)

a <- ggplot(data = dfplot, mapping = aes(x=hh.f25.bin,y=abs(delta_hh))) +
  geom_jitter(color='lightblue3',alpha=0.2, width=0.3) +
  geom_boxplot(fill="lemonchiffon3",color="black",alpha=0.3) +
  guides(color=FALSE) +
  theme_classic() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size=3) +
  ylab("Mean change in allele frequency") +
  xlab("Binned significance (-log10(P))") +
  scale_x_discrete(labels=plotlabs) +
  ggtitle("OWA") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


df.af$ah.f25.bin <- cut(-log10(df.af$ah_fdr),
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE)

dfplot <- df.af[!is.na(df.af$ah.f25.bin),]

b <- ggplot(data = dfplot, mapping = aes(x=ah.f25.bin,y=abs(delta_ah))) +
  geom_jitter(color='lightblue3',alpha=0.2, width=0.3) +
  geom_boxplot(fill="lemonchiffon3",color="black",alpha=0.3) +
  guides(color=FALSE) +
  theme_classic() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size=3) +
  ylab("Mean change in allele frequency") +
  xlab("Binned significance (-log10(P))") +
  scale_x_discrete(labels=plotlabs) +
  ggtitle("Acidification") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# set up cut-off values
#breaks <- c(0,2,4,6,8,10,12,14,16,18,20,22,24,150)
# specify interval/bin labels
#tags <- c("[0-2)","[2-4)", "[4-6)", "[6-8)", "[8-10)", "[10-12)","[12-14)", "[14-16)","[16-18)", "[18-20)", "[20-22)", "[22-24)", "[24-150)")
#plotlabs <- c("0-2","2-4", "4-6", "6-8", "8-10", "10-12","12-14", "14-16","16-18", "18-20", "20-22", "22-24", "24-150")


df.af$ha.f25.bin <- cut(-log10(df.af$ha_fdr),
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE)

dfplot <- df.af[!is.na(df.af$ha.f25.bin),]

c <- ggplot(data = dfplot, mapping = aes(x=ha.f25.bin,y=abs(delta_ha))) +
        geom_jitter(color='lightblue3',alpha=0.2, width=0.3) +
        geom_boxplot(fill="lemonchiffon3",color="black",alpha=0.3) +
        guides(color=FALSE) +
        theme_classic() +
        stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size=3) +
        ylab("Mean change in allele frequency") +
        xlab("Binned significance (-log10(P))") +
        scale_x_discrete(labels=plotlabs) +
        ggtitle("Warming") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))



df.af$aa.f25.bin <- cut(-log10(df.af$aa_fdr),
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE)

dfplot <- df.af[!is.na(df.af$aa.f25.bin),]

d <- ggplot(data = dfplot, mapping = aes(x=aa.f25.bin,y=abs(delta_aa))) +
        geom_jitter(color='lightblue3',alpha=0.2, width=0.3) +
        geom_boxplot(fill="lemonchiffon3",color="black",alpha=0.3) +
        guides(color=FALSE) +
        theme_classic() +
        stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size=3) +
        ylab("Mean change in allele frequency") +
        xlab("Binned significance (-log10(P))") +
        scale_x_discrete(labels=plotlabs) +
        ggtitle("Ambient") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("~/tonsa_genomics/figures/binned_afchange.png", 
            ggarrange(b,c,a,d), height=8 ,width=11, units="in")


############################################################
##############################
# binned prop af change
##############################
############################################################


########
# what about prop of total in each bin.

breaks <- seq(0,0.5,0.02)

# specify interval/bin labels
tags <- paste("[",breaks[1:length(breaks)-1],",",breaks[2:length(breaks)], ")", sep="")

plotlabs <- paste(breaks[1:length(breaks)-1]," - ",breaks[2:length(breaks)], sep="")

df.af$aa.f00.bin <- cut((df.af$aa.f0.af),
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE)

# make df of all groups, with sig cats. in format for ggplot.

pltdf <- as.data.frame(matrix(ncol=3, nrow=0))
colnames(pltdf) <- c("group", "bin", "proportion")

for(i in tags){
  tot_len <- nrow(df.af)
  aa_len <- length(which(df.af$aa_sig == TRUE))
  ah_len <- length(which(df.af$ah_sig == TRUE))
  ha_len <- length(which(df.af$ha_sig == TRUE))
  hh_len <- length(which(df.af$hh_sig == TRUE))

  bin_len <- length(which(df.af$aa.f00.bin == i))
  aa_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$aa_sig == TRUE))
  ah_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$ah_sig == TRUE))
  ha_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$ha_sig == TRUE))
  hh_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$hh_sig == TRUE))

  dftmp <- data.frame(
              group=c("all_loci","aa","ah", "ha", "hh"),
              bin=rep(i, 5),
              proportion = c(bin_len/tot_len,
                            aa_sig_len/aa_len,
                            ah_sig_len/ah_len,
                            ha_sig_len/ha_len,
                            hh_sig_len/hh_len)
              )
  pltdf <- rbind(pltdf, dftmp)

}

pltdf$group <- factor(pltdf$group, levels=c("all_loci","aa", "ah", "ha", "hh"))

# permutation

nrep <- 1000
nsamp <- sum(df.af$ah_sig == TRUE)

bs <- matrix(nrow=nsamp, ncol=nrep)
avg_rep <- c()
for(i in 1:nrep){
    bs[,i] <- sample(df.af$aa.f0.af, size=nsamp, replace=FALSE)
    avg_rep <- c(avg_rep, bs[,i])
}

out <- list()
proplist <- list()
for (j in 1:ncol(bs)){

    out[[j]] <- data.frame(x=bs[,j],
        bin = cut(bs[,j], breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE))
    tmpplt <- as.data.frame(matrix(ncol=3, nrow=0))
    # figure out prop for the sims
    for(i in tags){
        tot_len <- nrow(df.af)
        perm_len <- nsamp

        bin_len <- length(which(out[[j]]$bin == i))

        dftmp <- data.frame(
                    group=c("perm_loci"),
                    bin=i,
                    proportion = c(bin_len/perm_len)
                    )
        tmpplt <- rbind(tmpplt, dftmp)

    }

    proplist[[j]] <- tmpplt

}

out.new <- bind_rows(proplist)

densities.qtiles <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(proportion, 0.025),
            q50 = quantile(proportion, 0.5),
            q95 = quantile(proportion, 0.975))


pltperm <- rbind(data.frame(group=rep("permutations", 50) , bin=densities.qtiles$bin,
                  proportion=densities.qtiles$q50),
        pltdf)

pltdf <- filter(pltdf, group != "all_loci")


pltperm <- rbind(data.frame(group=rep("all_loci", 50) , bin=densities.qtiles$bin, proportion=densities.qtiles$q50),
        pltdf)

pltperm$group <- factor(pltperm$group, levels=c("all_loci","aa", "ah", "ha", "hh"))
pltperm$lower <- c(densities.qtiles$q05, rep(NA, nrow(pltperm)-nrow(densities.qtiles)))
pltperm$upper <- c(densities.qtiles$q95, rep(NA, nrow(pltperm)-nrow(densities.qtiles)))



a <- ggplot(data = pltperm, aes(x=bin,y=proportion, group=group, fill=group, shape=group)) +
  geom_col(color="black", position="dodge", size=0.3) +
  theme_classic() +
  #stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size=3) +
  ylab("Proportion of loci in bin") +
  xlab("Folded allele frequency bin") +
  scale_x_discrete(labels=plotlabs) +
  ggtitle("All significant, lab adaptation included") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=c("#D3DDDC","#6699CC","#F2AD00","#00A08A", "#CC3333"),
                    labels = c("All loci","Ambient","Acidic", "Warming", "Greenhouse"))+
   scale_shape_manual(values=c( 21,21,22, 23, 24)) +
#guides(fill=guide_legend(override.aes=list(
#        shape=c(19,19,19, 19, 19),
#        size=3,
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
#    fill=c('#D3DDDC',"#6699CC","#F2AD00","#00A08A", "#CC3333")),order = 2),
#    shape= FALSE)  +
theme(legend.title = element_blank(), legend.text=element_text(size=11),
        legend.position = c(0.8, 0.8))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.25, position=position_dodge(0.9))

ggsave("~/tonsa_genomics/figures/binned_proportions.pdf", a, h=4, w=10, units="in")

###########################
# remove lab effects
############################
pltdf <- as.data.frame(matrix(ncol=3, nrow=0))
colnames(pltdf) <- c("group", "bin", "proportion")

for(i in tags){
  tot_len <- nrow(df.af)
  aa_len <- length(which(df.af$aa_sig == TRUE))
  ah_len <- length(which(df.af$ah_sig == TRUE & df.af$aa_sig == FALSE))
  ha_len <- length(which(df.af$ha_sig == TRUE & df.af$aa_sig == FALSE))
  hh_len <- length(which(df.af$hh_sig == TRUE & df.af$aa_sig == FALSE))

  bin_len <- length(which(df.af$aa.f00.bin == i))
  aa_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$aa_sig == TRUE))
  ah_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$ah_sig == TRUE & df.af$aa_sig == FALSE))
  ha_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$ha_sig == TRUE & df.af$aa_sig == FALSE))
  hh_sig_len <- length(which(df.af$aa.f00.bin == i & df.af$hh_sig == TRUE & df.af$aa_sig == FALSE))

  dftmp <- data.frame(
              group=c("all_loci","aa","ah", "ha", "hh"),
              bin=rep(i, 5),
              proportion = c(bin_len/tot_len,
                            aa_sig_len/aa_len,
                            ah_sig_len/ah_len,
                            ha_sig_len/ha_len,
                            hh_sig_len/hh_len)
              )
  pltdf <- rbind(pltdf, dftmp)

}

pltdf$group <- factor(pltdf$group, levels=c("all_loci","aa", "ah", "ha", "hh"))


nrep <- 1000
nsamp <- sum(df.af$ah_sig == TRUE & df.af$aa_sig == FALSE)

bs <- matrix(nrow=nsamp, ncol=nrep)
avg_rep <- c()
for(i in 1:nrep){
    bs[,i] <- sample(df.af$aa.f0.af, size=nsamp, replace=FALSE)
    avg_rep <- c(avg_rep, bs[,i])
}

out <- list()
proplist <- list()
for (j in 1:ncol(bs)){

    out[[j]] <- data.frame(x=bs[,j],
        bin = cut(bs[,j], breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE))
    tmpplt <- as.data.frame(matrix(ncol=3, nrow=0))
    # figure out prop for the sims
    for(i in tags){
        tot_len <- nrow(df.af)
        perm_len <- nsamp

        bin_len <- length(which(out[[j]]$bin == i))


        dftmp <- data.frame(
                    group=c("perm_loci"),
                    bin=i,
                    proportion = c(bin_len/perm_len)
                    )
        tmpplt <- rbind(tmpplt, dftmp)

    }

    proplist[[j]] <- tmpplt

}

out.new <- bind_rows(proplist)

densities.qtiles <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(proportion, 0.025),
            q50 = quantile(proportion, 0.5),
            q95 = quantile(proportion, 0.975))


pltperm <- rbind(data.frame(group=rep("permutations", 50) , bin=densities.qtiles$bin, proportion=densities.qtiles$q50),
        pltdf)

pltdf <- filter(pltdf, group != "all_loci")


pltperm <- rbind(data.frame(group=rep("all_loci", 50) , bin=densities.qtiles$bin, proportion=densities.qtiles$q50),
        pltdf)

pltperm$group <- factor(pltperm$group, levels=c("all_loci","aa", "ah", "ha", "hh"))
pltperm$lower <- c(densities.qtiles$q05, rep(NA, nrow(pltperm)-nrow(densities.qtiles)))
pltperm$upper <- c(densities.qtiles$q95, rep(NA, nrow(pltperm)-nrow(densities.qtiles)))


b <- ggplot(data = pltperm, aes(x=bin,y=proportion, group=group, fill=group, shape=group)) +
  geom_col(color="black", position="dodge", size=0.3) +
  theme_classic() +
  #stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size=3) +
  ylab("Proportion of loci in bin") +
  xlab("Folded allele frequency bin") +
  scale_x_discrete(labels=plotlabs) +
  ggtitle("Lab adaptation loci removed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=c("#D3DDDC","#6699CC","#F2AD00","#00A08A", "#CC3333"),
                    labels = c("All loci","Ambient","Acidic", "Warming", "Greenhouse"))+
   scale_shape_manual(values=c(21,21,22, 23, 24)) +
#guides(fill=guide_legend(override.aes=list(
#        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
#    fill=c('#D3DDDC',"#6699CC","#F2AD00","#00A08A", "#CC3333")),order = 2),
#    shape= FALSE)  +
theme(legend.title = element_blank(), legend.text=element_text(size=11),
        legend.position = c(0.8, 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.25, position=position_dodge(0.9))


ggsave("~/tonsa_genomics/figures/binned_proportions_nolab.pdf", b, h=4, w=10, units="in")



########################################################################################
########################################################################################
########################################################################################
# find loci with large af changes.
########################################################################################
########################################################################################


annot2 <-data.frame(SNP=annot$SNP, class = annot$class, annotation = annot$annotation, distance = annot$distance, description = annot$description)


df.all <- full_join(df.af, annot2, by="SNP")

## find hh only large af change
df <- df.all[df.all$aa_fdr > 0.1,]
df <- df.all[df.all$aa_fdr > 0.05,]

# write function to calc df for large and sig af changes for whatever group of interest and output df

large_af_df <- function(group, threshold, dataf){
  tmp_fdr <- dataf[,paste0(group,"_fdr")]
  tmp_af <- dataf[,paste0("delta_",group)]
  sig_thresh <- quantile(-log10(tmp_fdr), threshold)
  af_thresh <- quantile(abs(tmp_af), threshold)
  #tmp_large <- dataf[-log10(tmp_fdr) > sig_thresh & abs(tmp_af) > af_thresh,]
  tmp_large <- dataf[-log10(tmp_fdr) > sig_thresh ,]
  #tmp_large <- (df[abs(tmp_af) > af_thresh,])


  return(tmp_large)

}



library(stringr)

hhlarge <- large_af_df(group="hh", threshold=0.99, dataf=df)
halarge <- large_af_df(group="ha", threshold=0.99, dataf=df)
ahlarge <- large_af_df(group="ah", threshold=0.99, dataf=df)
aalarge <- large_af_df(group="aa", threshold=0.99, dataf=df.all)

genes <- unique(unlist(strsplit(as.character(df.all$annotation), ",")))

####### AH ##########

pdf("~/tonsa_genomics/figures/candidates_ah.pdf", h=7, w=7)

num_plt <- 0
for(i in 2:length(genes)){
  tmpdf <- df.all[grep(genes[i], df.all$annotation),]
  if(i%%500 ==0){print(i)}
  if(sum(tmpdf$ah_sig) > 1 & sum(tmpdf$aa_sig) < 2){
        plot(  x=tmpdf$POS, y=-log10(tmpdf$ah_fdr), col="#F2AD00",  type = "l",
          ylab=c("-log10(p)"), xlab=c("Position"), ylim=c(0, 7),
          main=paste(unique(tmpdf$annotation),":",unique(tmpdf$description)))
        lines( x=tmpdf$POS, y=-log10(tmpdf$hh_fdr), col="#CC3333")
        points( x=tmpdf$POS, y=-log10(tmpdf$hh_fdr), bg="#CC3333", pch=24, cex=1.5)
        lines( x=tmpdf$POS, y=-log10(tmpdf$ha_fdr), col="#00A08A")
        points(x=tmpdf$POS, y=-log10(tmpdf$ha_fdr), bg="#00A08A", pch=23, cex=1.5)
        points(x=tmpdf$POS, y=-log10(tmpdf$ah_fdr), bg="#F2AD00", pch=22, cex=1.5)
        lines( x=tmpdf$POS, y=-log10(tmpdf$aa_fdr), col="#6699CC")
        points(x=tmpdf$POS, y=-log10(tmpdf$aa_fdr), bg="#6699CC", pch=21, cex=1.5)
        num_plt <- num_plt +1
        print(paste("plot num:",num_plt))
        print(paste("gene is:", genes[i]))
      }
  }
}

dev.off()


###### OWA ##########

pdf("~/tonsa_genomics/figures/candidates.pdf", h=7, w=7)

num_plt <- 0
for(i in 2:length(genes)){
  tmpdf <- df.all[grep(genes[i], df.all$annotation),]
  if(sum(tmpdf$hh_sig_nolab) > 3 | sum(tmpdf$ha_sig_nolab) > 3){
        if(sum(-log10(tmpdf$hh_fdr) > 3) > 3){


        plot(  x=tmpdf$POS, y=-log10(tmpdf$hh_fdr), col="#CC3333",  type = "l",
          ylab=c("-log10(p)"), xlab=c("Position"), ylim=c(0, 7),
          main=paste(unique(tmpdf$annotation),":",unique(tmpdf$description)))
        points( x=tmpdf$POS, y=-log10(tmpdf$hh_fdr), bg="#CC3333", pch=24, cex=1.5)
        lines( x=tmpdf$POS, y=-log10(tmpdf$ha_fdr), col="#00A08A")
        points(x=tmpdf$POS, y=-log10(tmpdf$ha_fdr), bg="#00A08A", pch=23, cex=1.5)
        lines( x=tmpdf$POS, y=-log10(tmpdf$ah_fdr), col="#F2AD00")
        points(x=tmpdf$POS, y=-log10(tmpdf$ah_fdr), bg="#F2AD00", pch=22, cex=1.5)
        lines( x=tmpdf$POS, y=-log10(tmpdf$aa_fdr), col="#6699CC")
        points(x=tmpdf$POS, y=-log10(tmpdf$aa_fdr), bg="#6699CC", pch=21, cex=1.5)
        num_plt <- num_plt +1
        print(paste("plot num:",num_plt))
        print(paste("gene is:", i))
      }
  }
}

dev.off()



genein <- c(

                            "TRINITY_DN121198_c0_g1",
                            "TRINITY_DN142673_c0_g3",
                            "TRINITY_DN134578_c0_g1",
                            "TRINITY_DN128043_c0_g1",
                            "TRINITY_DN145316_c0_g7",
                            "TRINITY_DN138308_c0_g6" 

                            )


# make new df with these candidate genes.

dfpl <- df.all[grep(paste(genein,collapse="|"),df.all$annotation), ]
nrow(dfpl)

# polarize by rising allele in hh


for(i in 1:nrow(dfpl)){
  if(dfpl$delta_hh[i] < 0){
    dfpl$aa.f0.af[i] <- 1-dfpl$aa.f0.af[i]
    dfpl$aa.f25.af[i] <- 1-dfpl$aa.f25.af[i]
    dfpl$hh.f25.af[i] <- 1-dfpl$hh.f25.af[i]
    dfpl$hh.f3.af[i] <- 1-dfpl$hh.f3.af[i]
    dfpl$ah.f25.af[i] <- 1-dfpl$ah.f25.af[i]
    dfpl$ha.f25.af[i] <- 1-dfpl$ha.f25.af[i]
  }
}
# recalc af change
dfpl$delta_aa    <- dfpl$aa.f25.af - dfpl$aa.f0.af
dfpl$delta_ah    <- dfpl$ah.f25.af - dfpl$aa.f0.af
dfpl$delta_ha    <- dfpl$ha.f25.af - dfpl$aa.f0.af
dfpl$delta_hh    <- dfpl$hh.f25.af - dfpl$aa.f0.af
dfpl$delta_hh_f3 <- dfpl$hh.f3.af  - dfpl$aa.f0.af

dfpl2 <- data.frame(POS=rep(dfpl$POS,4),
                    group = c(rep("Ambient", nrow(dfpl)),
                              rep("Acidification", nrow(dfpl)),
                              rep("Warming", nrow(dfpl)),
                              rep("Greenhouse", nrow(dfpl))
                              ),
                    pvalue = c(-log10(dfpl$aa_fdr),
                               -log10(dfpl$ah_fdr),
                               -log10(dfpl$ha_fdr),
                               -log10(dfpl$hh_fdr)
                               ),
                    F0_af = rep(dfpl$aa.f0.af, 4),
                    delta_af = c(dfpl$delta_aa, dfpl$delta_ah, dfpl$delta_ha, dfpl$delta_hh),
                    Gene = rep(dfpl$annotation,4)
                    )

dfpl2$F0_af <- ifelse(dfpl2$F0_af > 0.5, 1-dfpl2$F0_af, dfpl2$F0_af)


nmdf <- data.frame(gene = c(
                            "TRINITY_DN121198_c0_g1",
                            "TRINITY_DN142673_c0_g3",
                            "TRINITY_DN134578_c0_g1",
                            "TRINITY_DN128043_c0_g1",
                            "TRINITY_DN145316_c0_g7",
                            "TRINITY_DN138308_c0_g6" 
), # tatabox),

           name = c(
                    "Na+/K+ ATPase subunit alpha-4",
                    "DnaJ Hsp40 family member A3",
                    "NADH dehydrogenase 1 alpha",
                    "NADH dehydrogenase 1 beta",
                    "Focal adhesion kinase 1",
                    "ATP-dependent RNA helicase DDX39")
            )

dfpl2$Gene <- as.character(dfpl2$Gene)
dfpl2$Name <- NA
nmdf$name <- as.character(nmdf$name)

for (i in 1:length(nmdf$gene)){
  dfpl2$Name[grep(nmdf$gene[i],dfpl2$Gene)] <- nmdf$name[i]
}

dfpl2$group <- factor(dfpl2$group, levels = c("Ambient", "Acidification", "Warming", "Greenhouse"))

dfpl2 <- dfpl2[!(dfpl2$Name == "DnaJ Hsp40 family member A3" & dfpl2$POS > 1500),]

library(gtable)
library(cowplot)
library(grid)

shift_legend <- function(p){

  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }

  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")

  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")

  return(gp)
}


dfpl3 <- dfpl2
# remove low maf loci

dfpl2 <- dfpl2[dfpl2$F0_af > 0.1,]

p1 <- ggplot(dfpl2, aes(x=POS, y=delta_af, fill=group, shape=group, color=group)) +
        geom_line(show.legend=F, size=0.5, alpha=0.6) +
        geom_point(size=2.5, color="black") +
        theme_bw() +
        scale_shape_manual(values=c( 21,22, 23, 24))+
  scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"))+
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "Greenhouse"))+
          facet_wrap(~Name, scales="free_x")+ #labeller = label_wrap_gen(width=12))+
        theme(legend.title = element_blank())+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
      xlab("Position") +
      ylab("Allele frequency change from F0") +
      theme(strip.text.x = element_text(size = 8),
        axis.text.x = element_text(size =c(8)))+
       guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24), size=rep(5),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)

ggsave(grid.draw(shift_legend(p1)), file="~/tonsa_genomics/figures/candidate_genes_lines2.pdf",
        h=5, w=9)

ggsave((p1), file="~/tonsa_genomics/figures/candidate_genes_lines2.pdf",
        h=5, w=9)

p1 <- ggplot(dfpl2, aes(x=POS, y=delta_af, fill=group, shape=group, color=group)) +
        #geom_line(show.legend=F, size=0.5) +
        geom_point(size=2.5, color="black") +
        theme_bw() +
        scale_shape_manual(values=c( 21,22, 23, 24))+
  scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"))+
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "Greenhouse"))+
          facet_wrap(~Name, scales="free_x")+ #labeller = label_wrap_gen(width=12))+
        theme(legend.title = element_blank())+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
      xlab("Position") +
      ylab("Allele frequency change from F0") +
      theme(strip.text.x = element_text(size = 8),
        axis.text.x = element_text(size =c(8)))+
       guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24), size=rep(5),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)

ggsave(grid.draw(shift_legend(p1)), file="~/tonsa_genomics/figures/candidate_genes_nolines.pdf",
  h=5, w=9)


p1 <- ggplot(dfpl2, aes(x=POS, y=delta_af, fill=group, shape=group, color=group)) +
        geom_line(show.legend=F, size=0.5, alpha=0.5) +
        geom_point(size=2.5, color="black",aes(alpha=F0_af)) +
        theme_bw() +
        scale_shape_manual(values=c( 21,22, 23, 24))+
  scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"))+
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "Greenhouse"))+
          facet_wrap(~Name, scales="free_x")+ #labeller = label_wrap_gen(width=12))+
        theme(legend.title = element_blank())+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
      xlab("Position") +
      ylab("Allele frequency change from F0") +
      theme(strip.text.x = element_text(size = 8),
        axis.text.x = element_text(size =c(8)))+
       guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24), size=rep(5),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)


ggsave(grid.draw(shift_legend(p1)), file="~/tonsa_genomics/figures/candidate_genes_lines_alpha.pdf",
  h=5, w=9)


p1 <- ggplot(dfpl2, aes(x=POS, y=delta_af, fill=group, shape=group, color=group)) +
        #geom_line(show.legend=F, size=0.5, alpha=0.5) +
        geom_point(size=2.5, color="black",aes(alpha=F0_af)) +
        theme_bw() +
        scale_shape_manual(values=c( 21,22, 23, 24))+
  scale_color_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"))+
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "Greenhouse"))+
          facet_wrap(~Name, scales="free_x")+ #labeller = label_wrap_gen(width=12))+
        theme(legend.title = element_blank())+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
      xlab("Position") +
      ylab("Allele frequency change from F0") +
      theme(strip.text.x = element_text(size = 8),
        axis.text.x = element_text(size =c(8)))+
       guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24), size=rep(5),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)


ggsave(grid.draw(shift_legend(p1)), file="~/tonsa_genomics/figures/candidate_genes_nolines_alpha.pdf",
  h=5, w=7.20472)


