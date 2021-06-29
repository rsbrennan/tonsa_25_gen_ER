library(ggplot2)
library(ggpubr)

l1 <- read.csv("~/tonsa_genomics/analysis/count.aligned.txt", header=F)
colnames(l1) <- c("indiv","total", "mapped", "q20", "dups", "avg","tot2", "cov", "map_rate")
l2 <- read.csv("~/tonsa_genomics/analysis/count.aligned_lane2.txt", header=F)
colnames(l2) <- c("indiv","total", "mapped", "q20", "dups", "avg","tot2", "cov", "map_rate")

dat <- read.csv("~/tonsa_genomics/analysis/count.aligned_merged.txt", header=F)
colnames(dat) <- c("indiv","total", "mapped", "q20", "dups", "avg","tot2", "cov", "map_rate")

l1$total + l2$total == dat$total
# great.

all <- data.frame(indiv=l1$indiv, group=substr(l1$indiv, 1,6),
                  total_l1 = l1$total, total_l2 = l2$total,
                  mapped_l1 = l1$q20,mapped_l2 = l2$q20,
                  total_all = dat$total, mapped_all = dat$q20,
                  map_rate = dat$map_rate)


p1 <- ggplot(data=all,aes(x= indiv, y=mapped_l1, fill=group, group=indiv)) +
   geom_col(color="black") + theme_classic() +
   ylab("mapped reads")+
   ggtitle("lane1")+
   ylim(0, 30000000) +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=10),
          axis.title.y = element_text( size=16))

p2 <- ggplot(data=all,aes(x= indiv, y=mapped_l2, fill=group, group=indiv)) +
   geom_col(color="black") + theme_classic() +
   ylab("mapped reads")+
   ggtitle("lane2")+
   ylim(0, 30000000) +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=10),
          axis.title.y = element_text( size=16))

p3 <- ggplot(data=all,aes(x= indiv, y=mapped_all, fill=group, group=indiv)) +
   geom_col(color="black") + theme_classic() +
   ylab("mapped reads")+
   ggtitle("combined")+
   ylim(0, 30000000) +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=10),
          axis.title.y = element_text( size=16))
ggarrange(p1, p2, p3, nrow=1, common.legend=T)

ggsave("~/tonsa_genomics/figures/reads.png", ggarrange(p1, p2, p3, nrow=1, common.legend=T), h=6, w=12)
