# ld analysis, tonsa multigen genomics

library(scales)
library(plyr)
library(dplyr)

setwd("~/tonsa_genomics/analysis/ldx")
filelist = list.files(pattern = "*.ld.out")
data_list = lapply(filelist, read.table, header=FALSE)

file_name <- gsub(".ld.out", "", filelist)


names(data_list)<-file_name
colnames <- c("chr", "snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
dat <- lapply(data_list, setNames, colnames)


dat <- lapply(dat, function(x) transform(x, id1 = paste(chr, snp1, sep=":")))
dat <- lapply(dat, function(x) transform(x, id2 = paste(chr, snp2, sep=":")))

##############
###
### plot decay in ld
###
##############
    
# calc distance between snp1 and snp2
df <- lapply(dat, function(x) transform(x, distance = abs(x$snp1-x$snp2)))

df.1 <- df

df <- list()
for (i in 1:length(df.1)){

    df[[i]] <- df.1[[i]][which(df.1[[i]]$distance <=400),] 

}

names(df) <- names(dat)

# add read depth filter
df.1 <- df

df <- list()
for (i in 1:length(df.1)){

    df[[i]] <- df.1[[i]][which(df.1[[i]]$dp_intersect >5),] 

}

names(df) <- names(dat)





# difference in slopes
library(nlme)

model_out <- list()
for(i in 1:length(df)){
    model_out[[i]] <- lm(mle_est ~ log10(distance),data=df[[i]])
}



png("~/tonsa_genomics/figures/ld.png", height=120, width=150, units="mm", res=300)
par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(0,type='n', xlim=c(1,400), ylim=c(0,0.3),
    main="",
    ylab="",
    xlab="",
    cex.lab=0.9, cex.axis=0.7,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7, tcl=-0.2, at=c(1, 100, 200, 300, 400)) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Distance between SNPs in base pairs", line=1.5, cex.lab=0.9)
title(ylab="Estimated Linkage Disequilibrium", line=1.5, cex.lab=0.9)

for(i in 1:4){
    lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='#a8bcba')
}

for(i in 5:8){
    lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
    lwd=1.5, lty=1, col='#6699CC')
}

for(i in 9:12){
    lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
    lwd=1.5, lty=1, col='#F2AD00')
}

for(i in 13:16){
    lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
    lwd=1.5, lty=1, col='#00A08A')
}

#for(i in 17:20){
#   lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
#    lwd=1.5, lty=2, col='#CC3333')
#}

for(i in 21:24){
    lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='#CC3333')
}

for(i in 26:28){
    lines(sort(df[[i]]$distance, decreasing=FALSE), sort(fitted(model_out[[i]]), decreasing=TRUE),
    lwd=1.5, lty=1, col='#CC3333')
}


legend("topright",
    legend=c("Founding population",
            "Ambient",
            "Acidic",
            "Warming",
            "Greenhouse F3",
            "Greenhouse"),
       col=c("#a8bcba",'#6699CC',"#F2AD00","#00A08A", "#CC3333", "#CC3333"),
       lty=c(2,1,1,1,2,1), lwd=2.2, cex=0.75, bty = "n")
dev.off()



##############################################################################
#### stats
##############################################################################

bin_list <- list()

for(i in 1:length(df)){
    bin_list[[i]] <- transform(df[[i]], bin = cut(df[[i]]$distance, breaks=seq(from=0,to=500, by=10)))
    bin_list[[i]] <- transform(bin_list[[i]], bin1 = gsub( '\\(', "", sapply(strsplit(as.character(bin_list[[i]]$bin),","), '[', 1)))

}

for(i in 1:length(bin_list)){
    bin_list[[i]]$id <- names(df)[i]
}

library(dplyr)
library(ggplot2)

gdat <- bind_rows(bin_list)
gdat$group <- substr(gdat$id,1,6)


m1 <- lme(mle_est ~ log10(distance)*group,random=~1|id,data=gdat, method = "ML")
summary(m1)


library(emmeans)
emm = emmeans(m1, ~ log10(distance)*group)

pairs(emm)
# or for simple comparisons
pairs(emm, simple = "each")

 contrast         estimate      SE df t.ratio p.value
 AA_F00 - AA_F25 -0.019647 0.00154 21 -12.753 <.0001
 AA_F00 - AH_F25 -0.012294 0.00154 21  -7.981 <.0001
 AA_F00 - HA_F25 -0.017934 0.00154 21 -11.642 <.0001
 AA_F00 - HH_F00  0.000287 0.00154 21   0.187 1.0000
 AA_F00 - HH_F03 -0.000598 0.00154 21  -0.389 0.9997
 AA_F00 - HH_F25 -0.022036 0.00154 21 -14.301 <.0001
 AA_F25 - AH_F25  0.007353 0.00154 21   4.771 0.0017
 AA_F25 - HA_F25  0.001713 0.00154 21   1.112 0.9175
 AA_F25 - HH_F00  0.019935 0.00154 21  12.941 <.0001
 AA_F25 - HH_F03  0.019049 0.00154 21  12.365 <.0001
 AA_F25 - HH_F25 -0.002389 0.00154 21  -1.550 0.7135
 AH_F25 - HA_F25 -0.005640 0.00154 21  -3.660 0.0209
 AH_F25 - HH_F00  0.012581 0.00154 21   8.169 <.0001
 AH_F25 - HH_F03  0.011696 0.00154 21   7.593 <.0001
 AH_F25 - HH_F25 -0.009742 0.00154 21  -6.320 0.0001
 HA_F25 - HH_F00  0.018221 0.00154 21  11.830 <.0001
 HA_F25 - HH_F03  0.017335 0.00154 21  11.254 <.0001
 HA_F25 - HH_F25 -0.004102 0.00154 21  -2.661 0.1578
 HH_F00 - HH_F03 -0.000886 0.00154 21  -0.575 0.9969
 HH_F00 - HH_F25 -0.022323 0.00154 21 -14.489 <.0001
 HH_F03 - HH_F25 -0.021437 0.00154 21 -13.913 <.0001

# https://stats.stackexchange.com/questions/18391/how-to-calculate-the-difference-of-two-slopes

library(nlme)

din_aa_f25 <- gdat[which(gdat$group == "AA_F25"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_aa_f25)
anova(m1)
summary(m1)
# intercept: .24; distance: -0.05621456

din_aa_f00 <- gdat[which(gdat$group == "AA_F00"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_aa_f00)
summary(m1)
# intercept: .21; distance: -0.04972532

din_ah_f25 <- gdat[which(gdat$group == "AH_F25"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_ah_f25)
summary(m1)
# intercept: .23; distance: -0.05490757

din_ha_f25 <- gdat[which(gdat$group == "HA_F25"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_ha_f25)
summary(m1)
# intercept: 0.24; distance: -0.05561034

din_hh_f00 <- gdat[which(gdat$group == "HH_F00"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_hh_f00)
summary(m1)
# intercept: 0.215 ; distance: -0.05145203

din_hh_f03 <- gdat[which(gdat$group == "HH_F03"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_hh_f03)
summary(m1)
# intercept: 0.22; distance: -0.05444736

din_hh_f25 <- gdat[which(gdat$group == "HH_F25"),]
m1 <- lme(mle_est ~ log10(distance),random=~1|id,data=din_hh_f25)
summary(m1)
# intercept: .25; distance: -0.05802259


##############################################################################
#### plot by bin avg:
##############################################################################


bin_list <- list()

for(i in 1:length(df)){
    bin_list[[i]] <- transform(df[[i]], bin = cut(df[[i]]$distance, breaks=seq(from=0,to=500, by=10)))
    bin_list[[i]] <- transform(bin_list[[i]], bin1 = gsub( '\\(', "", sapply(strsplit(as.character(bin_list[[i]]$bin),","), '[', 1)))

}

for(i in 1:length(bin_list)){
    bin_list[[i]]$id <- names(df)[i]
}

library(dplyr)
library(ggplot2)

gdat <- bind_rows(bin_list)
gdat$group <- substr(gdat$id,1,6)

p <- ggplot(gdat, aes(x=distance, y=mle_est, shape=group, fill=group)) +
        stat_summary(fun.data="mean_se") +
        theme_bw() +
       scale_shape_manual(values=c( 21,21,22, 23, 24, 24, 24),
                    labels = c("Founding population-Ambient","Ambient", "Acidic", 
                                "Warming", "Founding population-Warm", "Warm- F3",
                               "OWA"))+
        scale_color_manual(values=c('black')) +
        xlim(0,300) +
          #scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A", "#CC3333"),
          scale_fill_manual(values=c(NA,'#6699CC',"#F2AD00",
                                     "#CC3333", NA, '#D3DDDC',
                                    "#00A08A"),
                    labels = c("Founding population-Ambient","Ambient", "Acidic", 
                                "Warming", "Founding population-Warm", "Warm- F3",
                               "OWA"))
ggsave(file="~/tonsa_genomics/figures/ld_bins.png", p, height=8, width=10)

