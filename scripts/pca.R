
library(ggplot2)


##### snps


af <- read.table("~/tonsa_genomics/analysis/filtered_allele_freqs.txt", header=TRUE)

pops <- c(
            "AA_F00_Rep1","AA_F00_Rep2","AA_F00_Rep3","AA_F00_Rep4",
            "AA_F25_Rep1","AA_F25_Rep2","AA_F25_Rep3","AA_F25_Rep4",
            "AH_F25_Rep1","AH_F25_Rep2","AH_F25_Rep3","AH_F25_Rep4",
            "HA_F25_Rep1","HA_F25_Rep2","HA_F25_Rep3","HA_F25_Rep4",
            "HH_F00_Rep1","HH_F00_Rep2","HH_F00_Rep3","HH_F00_Rep4",
            "HH_F03_Rep1","HH_F03_Rep2","HH_F03_Rep3","HH_F03_Rep4",
            "HH_F25_Rep1","HH_F25_Rep2","HH_F25_Rep3","HH_F25_Rep4")


freqs <- t(af[,2:ncol(af)])
colnames(freqs) <-  af$SNP

####
##
## plot pca
##
####

# f1 alone

pcaResult <- prcomp(freqs, scale=TRUE)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(id=pops, Line=substr(pops, 1,2),
        gen=substr(pops, 4,6),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

#data$Treatment <- gsub("AA", "AM", data$Treatment)

d <- ggplot(data, aes(PC1, PC2, fill=gen, shape=Line)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22, 23, 24))+
        scale_color_manual(values=c('black')) +
        scale_fill_manual(values=c('brown3','cornflowerblue', "darkorchid2"))+
        #theme(legend.position = c(0.88,0.17))+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
        guides(fill=guide_legend(override.aes=list(shape=21, fill=c('brown3', 'cornflowerblue','darkorchid2' )),order = 2))
                #shape=guide_legend(override.aes=list(shape=c(16, 17), size=c(7,5)),order = 1),
          #      size=FALSE)+
       # scale_size_manual(values=c(7,5)) + theme(plot.title = element_text(hjust = 0.5))

ggsave("~/tonsa_genomics/figures/pca.png",d)
