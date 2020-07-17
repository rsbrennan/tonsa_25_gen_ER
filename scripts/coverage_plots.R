library(ggplot2)

# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'
setwd("~/tonsa_genomics/analysis/coverage/all_probes")
# Get a list of the bedtools output files you'd like to read in
print(files <- list.files(pattern="all.txt$"))

print(labs <- paste(gsub("\\.bam\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))

#dat <- read.table("~/tonsa_genomics/analysis/coverage/all_probes/MeanCoverageBED.bedgraph", header=F)

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i], sep="\t", header=FALSE)
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
    print(paste(i, "done"))

}


cov <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i], sep="\t", header=FALSE)[,c(2,5)]
    cov_cumul=1-cumsum(cov[[i]][,2])
    cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")
cov_df$gp <- substr(cov_df$sample, 1,6)

p<-ggplot(data=cov_df,aes(x= depth, y=cov_cumul, color=gp, group=sample)) + 
    geom_vline(xintercept= 50, color="gray60") +
    geom_vline(xintercept= 75, color="gray60") +
    geom_hline(yintercept= 0.5, color="gray60") +
    geom_hline(yintercept= 0.9, color="gray60") +
    geom_line() + #xlim(0,200) +
    scale_color_brewer(palette="Set1") + theme_classic() +
    ylab("Fraction of probe bases \u2265 depth") +
    ggtitle("Coverage across all probes") +
    scale_x_continuous( breaks=c(0,50, 75, 100, 200), limits=c(0,200))+
    scale_y_continuous( breaks=c(0,.25,.5, .75, 0.9, 1))


ggsave("~/tonsa_genomics/figures/coverage_all_probes.png", p, h=4, w=5)


#### correlate with total reads:

total <- read.csv("~/tonsa_genomics/analysis/count.aligned.txt",header=F)

dp50 <- cov_df[which(cov_df$depth == 50),]

tot_all <- merge(dp50, total, by.x="sample", by.y="V1" )
colnames(tot_all) <- c("sample", "depth","fraction", "cov_cumul", "gp", "Total_reads", colnames(total)[3:ncol(total)])

plot(x=tot_all$Total_reads, y=tot_all$cov_cumul)

p <- ggplot(data=tot_all,aes(x= Total_reads, y=cov_cumul, fill=gp, group=sample)) +
     geom_point(size=4, pch=21, color="black") + theme_classic() +
     ylab("Proportion of probes with coverage >= 50x")


ggsave("~/tonsa_genomics/figures/cov_vs_depth.png", p, h=4, w=5)


#### plot total reads:

total <- read.csv("~/tonsa_genomics/analysis/count.aligned.txt",header=F)

dp50 <- cov_df[which(cov_df$depth == 50),]

tot_all <- merge(dp50, total, by.x="sample", by.y="V1" )
colnames(tot_all) <- c("sample", "depth","fraction", "cov_cumul", "gp", "Total_reads", colnames(total)[3:ncol(total)])

plot(x=tot_all$Total_reads, y=tot_all$cov_cumul)

p <- ggplot(data=tot_all,aes(x= sample, y=Total_reads, fill=gp, group=sample)) +
     geom_col(color="black") + theme_classic() +
     ylab("total number of reads")+
          theme(axis.text.x = element_text(angle=45, hjust=1, size=10),
          axis.title.y = element_text( size=16))
p

ggsave("~/tonsa_genomics/figures/reads.png", p, h=4, w=8)


####################
# tonsa
####################


setwd("~/tonsa_genomics/analysis/coverage/atonsa")

print(files <- list.files(pattern="all.txt$"))

print(labs <- paste(gsub("\\.bam\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
    print(paste(i, "done"))

}

cov <- list()
for (i in 1:length(files)) {
cov[[i]] <- read.table(files[i])[,c(2,5)]
cov_cumul=1-cumsum(cov[[i]][,2])
cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
cov[[i]]$sample=labs[i]
}
cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")
cov_df$gp <- substr(cov_df$sample, 1,6)

p<-ggplot(data=cov_df,aes(x= depth, y=cov_cumul, color=gp, group=sample)) + 
    geom_vline(xintercept= 50, color="gray60") +
    geom_vline(xintercept= 75, color="gray60") +
    geom_hline(yintercept= 0.5, color="gray60") +
    geom_hline(yintercept= 0.9, color="gray60") +
    geom_line() + #xlim(0,200) +
    scale_color_brewer(palette="Set1") + theme_classic() +
    ylab("Fraction of probe bases \u2265 depth") +
    ggtitle("Coverage across tonsa probes") +
    scale_x_continuous( breaks=c(0,50, 75, 100, 200), limits=c(0,200))+
    scale_y_continuous( breaks=c(0,.25,.5, .75, 0.9, 1))


ggsave("~/tonsa_genomics/figures/coverage_tonsa_probes.png", p, h=4, w=5)




####################
# promoter
####################


setwd("~/tonsa_genomics/analysis/coverage/promoter")

print(files <- list.files(pattern="all.txt$"))

print(labs <- paste(gsub("\\.bam\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
    print(paste(i, "done"))

}

cov <- list()
for (i in 1:length(files)) {
cov[[i]] <- read.table(files[i])[,c(2,5)]
cov_cumul=1-cumsum(cov[[i]][,2])
cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
cov[[i]]$sample=labs[i]
}
cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")
cov_df$gp <- substr(cov_df$sample, 1,6)

p<-ggplot(data=cov_df,aes(x= depth, y=cov_cumul, color=gp, group=sample)) + 
    geom_vline(xintercept= 50, color="gray60") +
    geom_vline(xintercept= 75, color="gray60") +
    geom_hline(yintercept= 0.5, color="gray60") +
    geom_hline(yintercept= 0.9, color="gray60") +
    geom_line() + #xlim(0,200) +
    scale_color_brewer(palette="Set1") + theme_classic() +
    ylab("Fraction of probe bases \u2265 depth") +
    ggtitle("Coverage across promoter probes") +
    scale_x_continuous( breaks=c(0,50, 75, 100, 200), limits=c(0,200))+
    scale_y_continuous( breaks=c(0,.25,.5, .75, 0.9, 1))


ggsave("~/tonsa_genomics/figures/coverage_promoter_probes.png", p, h=4, w=5)

####################
# both
####################

setwd("~/tonsa_genomics/analysis/coverage/both")

print(files <- list.files(pattern="all.txt$"))

print(labs <- paste(gsub("\\.bam\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
    print(paste(i, "done"))

}

cov <- list()
for (i in 1:length(files)) {
cov[[i]] <- read.table(files[i])[,c(2,5)]
cov_cumul=1-cumsum(cov[[i]][,2])
cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
cov[[i]]$sample=labs[i]
}
cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")
cov_df$gp <- substr(cov_df$sample, 1,6)

p<-ggplot(data=cov_df,aes(x= depth, y=cov_cumul, color=gp, group=sample)) + 
    geom_vline(xintercept= 50, color="gray60") +
    geom_vline(xintercept= 75, color="gray60") +
    geom_hline(yintercept= 0.5, color="gray60") +
    geom_hline(yintercept= 0.9, color="gray60") +
    geom_line() + #xlim(0,200) +
    scale_color_brewer(palette="Set1") + theme_classic() +
    ylab("Fraction of probe bases \u2265 depth") +
    ggtitle("Coverage across \"both\" probes") +
    scale_x_continuous( breaks=c(0,50, 75, 100, 200), limits=c(0,200))+
    scale_y_continuous( breaks=c(0,.25,.5, .75, 0.9, 1))


ggsave("~/tonsa_genomics/figures/coverage_both_probes.png", p, h=4, w=5)


####################
# hudsonica
####################

setwd("~/tonsa_genomics/analysis/coverage/hudsonica")

print(files <- list.files(pattern="all.txt$"))

print(labs <- paste(gsub("\\.bam\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
    print(paste(i, "done"))

}

cov <- list()
for (i in 1:length(files)) {
cov[[i]] <- read.table(files[i])[,c(2,5)]
cov_cumul=1-cumsum(cov[[i]][,2])
cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
cov[[i]]$sample=labs[i]
}
cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")
cov_df$gp <- substr(cov_df$sample, 1,6)

p<-ggplot(data=cov_df,aes(x= depth, y=cov_cumul, color=gp, group=sample)) + 
    geom_vline(xintercept= 50, color="gray60") +
    geom_vline(xintercept= 75, color="gray60") +
    geom_hline(yintercept= 0.5, color="gray60") +
    geom_hline(yintercept= 0.9, color="gray60") +
    geom_line() + #xlim(0,200) +
    scale_color_brewer(palette="Set1") + theme_classic() +
    ylab("Fraction of probe bases \u2265 depth") +
    ggtitle("Coverage across hudsonica probes") +
    scale_x_continuous( breaks=c(0,50, 75, 100, 200), limits=c(0,200))+
    scale_y_continuous( breaks=c(0,.25,.5, .75, 0.9, 1))


ggsave("~/tonsa_genomics/figures/coverage_hudsonica_probes.png", p, h=4, w=5)

