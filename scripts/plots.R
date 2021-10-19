library(ggplot2)
library(tidyr)

###########################
# PCA

##### snps


af <- read.table("~/tonsa_genomics/analysis/filtered_maf_freqs.txt", header=TRUE)

pops <- c(
            "AA_F00_Rep1","AA_F00_Rep2","AA_F00_Rep3","AA_F00_Rep4",
            "AA_F25_Rep1","AA_F25_Rep2","AA_F25_Rep3","AA_F25_Rep4",
            "AH_F25_Rep1","AH_F25_Rep2","AH_F25_Rep3","AH_F25_Rep4",
            "HA_F25_Rep1","HA_F25_Rep2","HA_F25_Rep3","HA_F25_Rep4",
            "HH_F00_Rep1","HH_F00_Rep2","HH_F00_Rep3","HH_F00_Rep4",
            "HH_F03_Rep1","HH_F03_Rep2","HH_F03_Rep3","HH_F03_Rep4",
            "HH_F25_Rep1","HH_F25_Rep2","HH_F25_Rep3","HH_F25_Rep4")

# remove F3 pop

af2 <- af[,grep("F03", colnames(af), invert=T)]
af2 <- af2[,grep("HH_F00", colnames(af2), invert=T)]
pops2 <- pops[grep("F03", pops, invert=T)]
pops2 <- pops2[grep("HH_F00", pops2, invert=T)]
pops <- pops2

varout <- apply(af2[,2:ncol(af2)], 1, var)

freqs <- t(af2[varout != 0,2:ncol(af2)])
colnames(freqs) <-  af$SNP[varout != 0]

nrow(freqs)


##
## plot pca
##

pcaResult <- prcomp(freqs, scale=T)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(id=pops, Line=substr(pops, 1,2),
        gen=substr(pops, 4,6),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

data$Line <- c(rep("Founding population", 4),
                rep("Ambient", 4),
                rep("Acidification", 4),
                rep("Warming", 4),
                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("Founding population","Ambient", "Acidification", "Warming", "OWA"))

data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=Line)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,21,22, 23, 24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Founding population","Ambient", "Acidification",
                               "Warming", "OWA"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank())+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23, 24),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)

ggsave("~/tonsa_genomics/figures/pca_afs_noF3.pdf",d, w=5.5, h=3.7)

#########################################################################
#########################################################################
# include F3
#########################################################################
#########################################################################
af <- read.table("~/tonsa_genomics/analysis/filtered_maf_freqs.txt", header=TRUE)

af2 <- af[,grep("HH_F00", colnames(af), invert=T)]
af <- af2


pops <- c(
            "AA_F00_Rep1","AA_F00_Rep2","AA_F00_Rep3","AA_F00_Rep4",
            "AA_F25_Rep1","AA_F25_Rep2","AA_F25_Rep3","AA_F25_Rep4",
            "AH_F25_Rep1","AH_F25_Rep2","AH_F25_Rep3","AH_F25_Rep4",
            "HA_F25_Rep1","HA_F25_Rep2","HA_F25_Rep3","HA_F25_Rep4",
            "HH_F03_Rep1","HH_F03_Rep2","HH_F03_Rep3","HH_F03_Rep4",
            "HH_F25_Rep1","HH_F25_Rep2","HH_F25_Rep3","HH_F25_Rep4")


freqs <- t(af[,2:ncol(af)])
colnames(freqs) <-  af$SNP

nrow(freqs)



####
##
## plot pca
##
####

pcaResult <- prcomp(freqs, scale=T)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(id=pops, Line=substr(pops, 1,2),
        gen=substr(pops, 4,6),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

data$Line <- c(rep("Founding population", 4),
                rep("Ambient", 4),
                rep("Acidification", 4),
                rep("Warming", 4),
                rep("OWA F3", 4),
                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("Founding population","Ambient", "Acidification",
                            "Warming","OWA F3", "OWA"))

#data$PC1[18] <- data$PC1[18] - 10



d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=Line)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,21,22, 23,24, 24))+
        scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00",
                            "#00A08A", "lightcoral", "#CC3333"),
                    labels = c("Founding ambient","Ambient", "Acidification",
                               "Warming","OWA F3", "OWA"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank())+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23,24, 24),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00",
                            "#00A08A", "lightcoral", "#CC3333")),order = 2),
    shape= FALSE)

ggsave("~/tonsa_genomics/figures/pca_afs.pdf",d, w=5.5, h=3.7)




################################################################################################################
################################################################################################################
# cvtk plots
################################################################################################################
################################################################################################################


########################################################
# shared variance:
########################################################


pdat_labs <- c("Total", "Lab\nadaptation\nremoved")

##### snps
dat <- read.csv("~/tonsa_genomics/analysis/total_variance.csv", header=TRUE)

# remove rows I don't need for plotting:

dat <- dat[grep("total", dat$var_group, invert=T),]
dat <- dat[grep("lab", dat$var_group, invert=T),]
# replace names:
dat$category <- dat$var_group

dat$category <- rep(c("total", "lab removed"), 3)

dat$category <- factor(dat$category, levels = c("total", "lab removed"))

dat$treatment <- c("OWA", "OWA", "Warming", "Warming", "Acidification", "Acidification")
dat$treatment <- factor(dat$treatment, levels = c("Acidification", "Warming", "OWA"))
#dat$cat_share <- factor(dat$cat_share, levels = c("shared_Acidic", "shared_Warm", "shared_Greenhoue","no_lab_Acidic", "no_lab_Warm", "no_lab_Greenhoue"))
d <- ggplot(dat, aes(x=category, y = estimate, fill=treatment, shape=treatment)) +
        stat_summary(fun.data=mean_cl_boot, size=0.5, geom="line", position = position_dodge(width = 0), aes(group=treatment)) +
        geom_point(size=5, position = position_dodge(width=0)) +
  geom_point(data = subset(dat, treatment == 'OWA'),
                aes(x = category, y = estimate, fill=treatment, shape=treatment),
                size=5,position = position_dodge(width=0)) +
        geom_errorbar(aes(ymin=lower_err, ymax=upper_err), width=0.1, size=0.5, position = position_dodge(width = 0)) +
        #geom_point(size=3,  position = position_jitterdodge(), alpha=0.2) +

        xlab(" ") +
        ylab("Proportion of total variance\ndue to selection") +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 22,23, 24))+
        scale_colour_manual(values=c('black')) +
        scale_fill_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                labels = c("Acidification",
                            "Warming", "OWA"))+
        #facet_wrap(~treatment) +
         guides(fill = guide_legend(override.aes = list(size=5))) +
        theme(legend.title = element_blank())+
        #geom_hline(yintercept=0,
        #        color = "black", size=0.5)+
        #theme(legend.position = "none") +
        ylim(0,70) + scale_x_discrete(labels= pdat_labs)+
        theme(axis.title=element_text(size=12),
            axis.text=element_text(size=10),
            strip.text = element_text(size=12))

ggsave("~/tonsa_genomics/figures/total_vs_lab_variance_cvtk.pdf",d, w=4, h=3.4)


### total variance

dat <- read.csv("~/tonsa_genomics/analysis/total_variance.csv", header=TRUE)
dat <- dat[grep("total", dat$var_group, invert=F),]
# replace names:
dat$category <- dat$var_group

dat$treatment <- c("OWA", "Warming", "Acidification")
dat$treatment <- factor(dat$treatment, levels = c("Acidification", "Warming", "OWA"))
#dat$cat_share <- factor(dat$cat_share, levels = c("shared_Acidic", "shared_Warm", "shared_Greenhoue","no_lab_Acidic", "no_lab_Warm", "no_lab_Greenhoue"))
d <- ggplot(dat, aes(x=treatment, y = estimate, fill=treatment, shape=treatment)) +
        geom_point(size=5, position = position_dodge(width=0)) +
        geom_errorbar(aes(ymin=lower_err, ymax=upper_err), width=0.1, size=0.5, position = position_dodge(width = 0)) +
        #geom_point(size=3,  position = position_jitterdodge(), alpha=0.2) +

        xlab("") +
        ylab("Total variance in\nallele frequency change") +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 22,23, 24))+
        scale_colour_manual(values=c('black')) +
        scale_fill_manual(values=c("#F2AD00","#00A08A", "#CC3333"),
                labels = c("Acidification",
                            "Warming", "OWA"))+
        #facet_wrap(~treatment) +
         guides(fill = guide_legend(override.aes = list(size=5))) +
        theme(legend.title = element_blank())+
        #geom_hline(yintercept=0,
        #        color = "black", size=0.5)+
        #theme(legend.position = "none") +
        theme(axis.title=element_text(size=12),
            axis.text=element_text(size=10),
            strip.text = element_text(size=12))

ggsave("~/tonsa_genomics/figures/total_variance_cvtk.pdf",d, w=4, h=3.4)





#################
#
# shared variance between treatments
#
#################

dat <- read.csv("~/tonsa_genomics/analysis/shared_variance.csv", header=TRUE)

dat <- dat[grep("total", dat$var_group, invert=T),]
dat <- dat[grep("lab", dat$var_group, invert=T),]

dat$treatment <- c("OWA vs.\nWarming","OWA vs.\nWarming",
                    "Acidification vs.\nWarming", "Acidification vs.\nWarming",
                    "OWA vs.\nAcidification", "OWA vs.\nAcidification")

dat$category <- rep(c("total", "lab removed"), 3)
dat$category <- factor(dat$category, levels = c("total", "lab removed"))

d <- ggplot(dat, aes(x=category, y = estimate, color=treatment,fill=treatment, shape=treatment)) +
        geom_hline(yintercept=0,
                color = "gray40", size=0.5, linetype="dashed")+
     stat_summary(fun.data=mean_cl_boot, size=0.5, geom="line", position = position_dodge(width = 0), aes(group=treatment)) +
        geom_point(size=4, position = position_dodge(width=0)) +
        geom_errorbar(aes(ymin=lower_err, ymax=upper_err), width=0.1, size=0.5, position = position_dodge(width = 0), color="black") +

       #         stat_summary(fun.data=mean_cl_boot, size=1, geom="line", position = position_dodge(width = 0.9), aes(group=treatment)) +
        #stat_summary(size=1.3, color="black", position = position_dodge(width = 0.9), stroke=0.5,
                        #fun.data=mean_cl_boot) +
        xlab(" ") +
        ylab("Proportion of total variance shared") +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,21,21))+
        scale_colour_manual(values=c('gray40','gray40','gray40')) +
        scale_fill_manual(values=c('gray70',"gray70", "gray70"),
                labels = c("Acidic",
                            "Warming", "OWA"))+
        facet_wrap(~treatment) +

        theme(legend.title = element_blank())+

        theme(legend.position = "none")+
        ylim(-2,50) + scale_x_discrete(labels= pdat_labs)+
        theme(axis.title=element_text(size=12),
            axis.text=element_text(size=10),
            strip.text = element_text(size=12))


ggsave("~/tonsa_genomics/figures/between_trt_vs_lab_variance_cvtk.pdf",d, w=5.5, h=3.4)



######################################################
# stacked bar plot
######################################################

dat <- read.csv("~/tonsa_genomics/analysis/total_variance.csv", header=TRUE)

# remove rows I don't need for plotting:

dat <- dat[grep("total", dat$var_group, invert=T),]
dat <- dat[grep("shared", dat$var_group, invert=T),]

#dat$var_group <- c("")

driftdf <- data.frame(treatment = c("OWA", "warm", "acidic"),
            var_group = c("drift", "drift", "drift"),
            lower_err = c(0,0,0),
            estimate= c(100-66.888, 100-64.298, 100-32.176),
            upper_err = c(0,0,0))


df <- rbind(dat, driftdf)
dat <- df

dat$category <- dat$var_group

dat$treatment <- c("OWA","OWA", "Warming", "Warming",
                    "Acidification", "Acidification", "OWA", "Warming", "Acidification")
dat$treatment <- factor(dat$treatment, levels = c("Acidification", "Warming", "OWA"))
#dat$cat_share <- factor(dat$cat_share, levels = c("shared_Acidic", "shared_Warm", "shared_Greenhoue","no_lab_Acidic", "no_lab_Warm", "no_lab_Greenhoue"))

library(RColorBrewer)

d <- ggplot(dat, aes(x=treatment, y = estimate, fill=var_group, shape=treatment)) +
        geom_col() +
        scale_fill_manual(values = (brewer.pal(n = 3, name = "Set2"))) +
        theme_bw() +
        #geom_point(size=3,  position = position_jitterdodge(), alpha=0.2) +
        xlab("") +
        ylab("Percent of total") +
       # ylim(-30, 23) + xlim(-50, 65)+
        #facet_wrap(~treatment) +
         guides(fill = guide_legend(override.aes = list(size=5))) +
        theme(legend.title = element_blank())+
        #geom_hline(yintercept=0,
        #        color = "black", size=0.5)+
        #theme(legend.position = "none") +
        theme(axis.title=element_text(size=12),
            axis.text=element_text(size=10),
            strip.text = element_text(size=12))


ggsave("~/tonsa_genomics/figures/all_variance.pdf",d, w=4, h=3.4)


################################
# covariance plots
################################

dat <- read.csv("~/tonsa_genomics/analysis/covariance_pairwise.csv", header=T)

dat$labels <- gsub("AA", "Ambient",dat$labels)
dat$labels <- gsub("AH", "Acidification",dat$labels)
dat$labels <- gsub("HA", "Warming",dat$labels)
dat$labels <- gsub("HH", "OWA",dat$labels)

dat1 <- dat %>% separate(labels, sep=" ", c("id_1", "id_2"),remove = FALSE)

dat1 <- dat1 %>% separate(id_1, sep="-", c("pop1"),remove = FALSE) %>%
                    separate(id_2, sep="-", c("pop2"),remove = FALSE)


# want to sort so AA's come first,

odat <- as.data.frame(matrix(nrow = 0, ncol=ncol(dat1)))
colnames(odat) <- colnames(dat1)

tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Ambient"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Acidification"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "Acidification"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "Warming" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Warming" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "OWA" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)


# Plot
odat$plt_order <- seq(1, nrow(odat), 1)

odat$comparison <- "inter-line"
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "Ambient")] <- "Ambient"
odat$comparison[which(odat$pop1 == "Acidification" & odat$pop2 == "Acidification")] <- "Acidification"
odat$comparison[which(odat$pop1 == "Warming" & odat$pop2 == "Warming")] <- "Warming"
odat$comparison[which(odat$pop1 == "OWA" & odat$pop2 == "OWA")] <- "OWA"

odat$comparison <- factor(odat$comparison, levels=c("Ambient", "Acidification","Warming", "OWA", "inter-line"))


p <- ggplot(odat, aes(x=plt_order, y=mean, fill=comparison)) +
  geom_hline(yintercept=0, color="grey50")+
  geom_segment( aes(x=plt_order, xend=plt_order, y=0, yend=mean), color="grey") +
  geom_errorbar(aes(ymin=mean-lower, ymax=mean+lower), width=.1) +
  geom_point( size=4, pch=21) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_x_continuous(breaks=odat$plt_order, labels=odat$labels)+
    theme(axis.text.x=element_text(angle=90, vjust=.5)) +
  xlab("") +
  ylab("Covariance") +
        scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333","gray78"),
                labels = c("Ambient", "Acidic",
                            "Warming", "OWA", "inter-line"))

ggsave("~/tonsa_genomics/figures/covariance_pairwise.pdf",p, h=5, w=14)


########################################################################################
############################################
# convergence correlation
############################################
########################################################################################

dat <- read.csv("~/tonsa_genomics/scripts/cvtk/combined_conv_corrs.csv", header=F)

colnames(dat) <- c("lower", "mean", "upper", "labels")
dat$labels <- gsub("AA", "Ambient",dat$labels)
dat$labels <- gsub("AH", "Acidification",dat$labels)
dat$labels <- gsub("HA", "Warming",dat$labels)
dat$labels <- gsub("HH", "OWA",dat$labels)

dat1 <- dat %>% separate(labels, sep=" ", c("id_1", "id_2"),remove = FALSE)

dat1 <- dat1 %>% separate(id_1, sep="-", c("pop1", "rep1"),remove = FALSE) %>%
                    separate(id_2, sep="-", c("pop2", "rep2"),remove = FALSE)
# want to sort so AA's come first,

odat <- as.data.frame(matrix(nrow = 0, ncol=ncol(dat1)))
colnames(odat) <- colnames(dat1)

tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Ambient"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Acidification"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "Acidification"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "Warming" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Warming" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "OWA" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)


# Plot
odat$plt_order <- seq(1, nrow(odat), 1)

odat$comparison <- "inter-line"
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "Ambient")] <- "Ambient"
odat$comparison[which(odat$pop1 == "Acidification" & odat$pop2 == "Acidification")] <- "Acidification"
odat$comparison[which(odat$pop1 == "Warming" & odat$pop2 == "Warming")] <- "Warming"
odat$comparison[which(odat$pop1 == "OWA" & odat$pop2 == "OWA")] <- "OWA"

odat$comparison <- factor(odat$comparison, levels=c("Ambient", "Acidification","Warming", "OWA", "inter-line"))



p <- ggplot(odat, aes(x=plt_order, y=mean, fill=comparison)) +
  geom_hline(yintercept=0, color="grey50")+
  geom_segment( aes(x=plt_order, xend=plt_order, y=0, yend=mean), color="grey") +
  geom_point( size=4, pch=21) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_x_continuous(breaks=odat$plt_order, labels=odat$labels)+
    theme(axis.text.x=element_text(angle=90, vjust=.5)) +
  xlab("") +
  ylab("Convergent correlation") +
        scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333","gray78"),
                labels = c("Ambient", "Acidic",
                            "Warming", "OWA", "inter-line"))

ggsave("~/tonsa_genomics/figures/convergence_correlation_pairwise.pdf",p, h=5, w=14)


###############
############### drop the shared replicates out:

dat1 <- dat1[dat1$rep1 != dat1$rep2,]

odat <- as.data.frame(matrix(nrow = 0, ncol=ncol(dat1)), stringsAsFactors=F)
colnames(odat) <- colnames(dat1)

tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Ambient"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Acidification"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Ambient" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "Acidification"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Acidification" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "Warming" & dat1$pop2 == "Warming"),]
odat <- rbind(odat, tmpdf)
tmpdf <- dat1[which(dat1$pop1 == "Warming" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)

tmpdf <- dat1[which(dat1$pop1 == "OWA" & dat1$pop2 == "OWA"),]
odat <- rbind(odat, tmpdf)


# Plot
odat$plt_order <- seq(1, nrow(odat), 1)

odat$comparison <- "inter-line"
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "Ambient")] <- "Ambient"
odat$comparison[which(odat$pop1 == "Acidification" & odat$pop2 == "Acidification")] <- "Acidification"
odat$comparison[which(odat$pop1 == "Warming" & odat$pop2 == "Warming")] <- "Warming"
odat$comparison[which(odat$pop1 == "OWA" & odat$pop2 == "OWA")] <- "OWA"

odat$comparison <- factor(odat$comparison, levels=c("Ambient", "Acidification","Warming", "OWA", "inter-line"))


p <- ggplot(odat, aes(x=plt_order, y=mean, fill=comparison)) +
  geom_hline(yintercept=0, color="grey50")+
  geom_segment( aes(x=plt_order, xend=plt_order, y=0, yend=mean), color="grey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  geom_point( size=4, pch=21) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_x_continuous(breaks=odat$plt_order, labels=odat$labels)+
    theme(axis.text.x=element_text(angle=90, vjust=.5)) +
  xlab("") +
  ylab("Convergence correlation") +
        scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333","gray78"),
                labels = c("Ambient", "Acidification",
                            "Warming", "OWA", "Inter-line"))

ggsave("~/tonsa_genomics/figures/convergence_correlation_pairwise_shared_rm.pdf",p, h=5, w=14)

#######
# merged plot
#######
odat$comparison <- c("inter-line")
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "Ambient")] <- "Ambient"
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "Acidification")] <- "Ambient vs. Acidification"
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "Warming")] <- "Ambient vs. Warming"
odat$comparison[which(odat$pop1 == "Ambient" & odat$pop2 == "OWA")] <- "Ambient vs. OWA"
odat$comparison[which(odat$pop1 == "Acidification" & odat$pop2 == "Acidification")] <- "Acidification"
odat$comparison[which(odat$pop1 == "Acidification" & odat$pop2 == "Warming")] <- "Acidification vs. Warming"
odat$comparison[which(odat$pop1 == "Acidification" & odat$pop2 == "OWA")] <- "Acidification vs. OWA"
odat$comparison[which(odat$pop1 == "Warming" & odat$pop2 == "Warming")] <- "Warming"
odat$comparison[which(odat$pop1 == "Warming" & odat$pop2 == "OWA")] <- "Warming vs. OWA"
odat$comparison[which(odat$pop1 == "OWA" & odat$pop2 == "OWA")] <- "OWA"


odat$comparison <- factor(odat$comparison, levels = c("Ambient", "Ambient vs. Acidification",
                          "Ambient vs. Warming", "Ambient vs. OWA",
                                                        "Acidification", "Acidification vs. Warming", "Acidification vs. OWA",
                                                        "Warming","Warming vs. OWA", "OWA"))

labs <- c("Ambient", "Ambient\nvs.\nAcidic", "Ambient\nvs.\nWarming", "Ambient\nvs.\nOWA",
                                                        "Acidic", "Acidic\nvs.\nWarming", "Acidic\nvs.\nOWA",
                                                        "Warming","Warming\nvs.\nOWA", "OWA")


odat$comp_cols <- "inter-line"
odat$comp_cols[which(odat$pop1 == "Ambient" & odat$pop2 == "Ambient")] <- "Ambient"
odat$comp_cols[which(odat$pop1 == "Acidification" & odat$pop2 == "Acidification")] <- "Acidification"
odat$comp_cols[which(odat$pop1 == "Warming" & odat$pop2 == "Warming")] <- "Warming"
odat$comp_cols[which(odat$pop1 == "OWA" & odat$pop2 == "OWA")] <- "OWA"

odat$comp_cols <- as.factor(odat$comp_cols)

odat$comp_cols <- factor(odat$comp_cols, levels=c("Ambient", "Acidification","Warming", "OWA", "inter-line"))

library(ggbeeswarm)
p <- ggplot(odat, aes(x=comparison, y=mean, fill=comp_cols, shape=comp_cols)) +
  #geom_hline(yintercept=0, color="grey50")+
  #geom_segment( aes(x=plt_order, xend=plt_order, y=0, yend=mean), color="grey") +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, alpha=0.3) +
  geom_violin(alpha=0.3, show.legend=F) +
    stat_summary(size=6, fun = mean,
                 geom = "point")  +
  geom_beeswarm( size=2, alpha=0.3, width=0.1, show.legend=F) +

    theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
  ) +
        scale_shape_manual(values=c(21,22,23,24, 21))+
   # theme(axis.text.x=element_text(angle=90, vjust=.5)) +
  xlab("") +
  ylab("Convergent correlation") +
        scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333","gray78"),
                labels = c("Ambient", "Acidic",
                            "Warming", "OWA", "Inter-line")) +
  guides(fill=guide_legend(override.aes=list(
        shape=c(21,22, 23, 24, 21), size=rep(5.5),
    fill=c('#6699CC',"#F2AD00","#00A08A", "#CC3333","gray78")),order = 2),
    shape= FALSE) +
  theme(legend.title=element_blank()) +
  scale_x_discrete(labels = labs)


ggsave("~/tonsa_genomics/figures/convergence_correlation_pairwise_grouped.pdf",p, h=4, w=8)

