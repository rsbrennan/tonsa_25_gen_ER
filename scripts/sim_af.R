library(poolSeq)
library(scales)

source("~/tonsa_genomics/scripts/sim_function.R")

#dat <- read.table("~/tonsa_genomics/analysis/variants.sync", header=F)

#synchronized (sync) files contain allele frequencies at specific genomic loci in multiple populations. Suppose you want to load a sync-file containing allele frequency trajectories of 2 populations (F0.R1, F10.R1, F20.R1, F0.R2, F10.R2, F20.R2). The following command allows you to read such file with poolSeq:

allSync <- read.sync(file="~/tonsa_genomics/analysis/variants.sync", polarization="reference" , keepOnlyBiallelic = FALSE,
              gen=c(rep(0,28)),
              repl=c(seq(1,28,1))) # note that the labels are incorrect, just getting it read in

aaSync <- read.sync(file="~/tonsa_genomics/analysis/AA.sync", polarization="reference" , keepOnlyBiallelic = FALSE,
              gen=c(rep(0,4), rep(25,4)),
              repl=c(seq(1,4,1),seq(1,4,1)))

ahSync <- read.sync(file="~/tonsa_genomics/analysis/AH.sync", polarization="reference",
              gen=c(rep(0,4), rep(25,4)),
              repl=c(seq(1,4,1),seq(1,4,1)))

haSync <- read.sync(file="~/tonsa_genomics/analysis/HA.sync", polarization="reference",
              gen=c(rep(0,4), rep(25,4)),
              repl=c(seq(1,4,1),seq(1,4,1)))

hhSync <- read.sync(file="~/tonsa_genomics/analysis/HH.sync", polarization="reference",
              gen=c(rep(0,4), rep(25,4)),
              repl=c(seq(1,4,1),seq(1,4,1)))


listSync <- list(aaSync, ahSync, haSync, hhSync)
names(listSync) <- c("aaSync", "ahSync", "haSync", "hhSync")
#   af(sync, chr, pos, repl, gen)
# af(mySync, "LS387016.1", 58, 1, 0)
# coverage(mySync, "LS387016.1", 58, 1, 0)
# af.traj(mySync, "LS387016.1", 58, 1)

## estimate effective population size:

NeOut <- as.data.frame(matrix(nrow=16, ncol=3))
colnames(NeOut) <- c("treatment", "replicate", "Ne")

for(i in 1:length(listSync)){
  for(j in 1:4){
        tmpNe <- substr(names(listSync)[i], 1,2)
        tmpNe[2] <- j
        tmpNe[3] <- estimateNe(p0=af(listSync[[i]],,, j, 0), pt=af(listSync[[i]],,, j, 25),
              cov0=coverage(listSync[[i]],,, j, 0), covt=coverage(listSync[[i]],,, j, 25),
              t=25, Ncensus=3000, poolSize=c(50, 50))
       NeOut[min(which(is.na(NeOut$treatment))),] <- tmpNe
  }

}

ne_in <- median(as.numeric(NeOut$Ne))
# [1] 413.97

write.table(NeOut, "~/tonsa_genomics/analysis/Ne_est.txt", quote=FALSE, sep="\t")

# run simulation based on the mean starting allele frequency of AA and HH (check this correlation)
## Ne is based on above

# mean of AA af

AA_t0 <- (rowMeans(af(listSync[['aaSync']],,, 1:4, 0)))
HH_t0 <- (rowMeans(af(listSync[['hhSync']],,, 1:4, 0)))
aa.0 <- data.frame(snp= names(AA_t0), af = AA_t0)
hh.0 <- data.frame(snp= names(HH_t0), af = HH_t0)

startaf <- merge(aa.0, hh.0, by="snp")
nrow(startaf)

summary(lm(startaf$af.x~startaf$af.y))
#R-squared = 0.9854

png("~/tonsa_genomics/figures/F0_af_corr.png",
      res=300, height=6, width=6, units="in")
plot(startaf$af.x,startaf$af.y, pch=19, col=alpha("black", 0.1),
  xlab=c("AA F0 mean allele frequency"),
  ylab=c("HH F0 mean allele frequency"),
  main=c("R-squared = 0.9868"))
abline(lm(startaf$af.x~startaf$af.y), col="firebrick3", lwd=3, lty=2)
#points(x=startaf$af.x[which(startaf$af.y > 0.65)],
#       y=startaf$af.y[which(startaf$af.y > 0.65)], pch=19, col="red")

dev.off()

# whats the mean depth for t0 and t25?

mCov <- read.table("~/tonsa_genomics/analysis/mean_coverage.txt", header=T)

f0_cov <- mean(mCov$coverage[c(1:4,17:20)])
f25_cov <- mean(mCov$coverage[c(5:16,25:28)])


#############
# simulate data
#############

# get starting allele freq:

mean_start_af <- rowMeans(af(allSync,,, c(1:4,17:20), 0))


# each run takes about 50 seconds

n=500
for (i in 1:n){
  outt <- simulate_allele_freq(startingAF= mean_start_af, effective_pop= round(ne_in), T0=0, T_N=25, samp_size=50, census_size=3000,
                          F0_depth = round(f0_cov), FN_depth = round(f25_cov))
  write.table(outt, file=paste('~/tonsa_genomics/analysis/cmh_simulation/simulated.',i, '.sync',sep=""), sep = '\t',
                            col.names = FALSE, row.names = FALSE, quote=FALSE)

if(i%%25 ==0){print(i)}
}


#### Simulartion for F25 comparisons:


n=500
for (i in 1:n){
  outt <- simulate_allele_f25(startingAF= mean_start_af, effective_pop= round(ne_in), T0=0, T_N=25, samp_size=50, census_size=3000,
                          F0_depth = round(f0_cov), FN_depth = round(f25_cov))
 write.table(outt, file=paste('~/tonsa_genomics/analysis/cmh_simulation/simulated.f25.',i, '.sync',sep=""), sep = '\t',
                            col.names = FALSE, row.names = FALSE, quote=FALSE)

if(i%%25 ==0){print(i)}
}



# plot distribution of differences in allele frequency before and after sampling
png("~/tonsa_genomics/figures/simulate_cov_error.png",
      res=300, height=7, width=11, units="in")
par(mfrow=c(2,2))

hist(simTraj_1b$F0-simTraj_1b$F0.ind.adj, main="F0 error: sampling individuals", xlab="Error in allele frequency (%)", ylab="Occurrences",
        xlim=c(-.3,.3), col="grey70")
hist(simTraj_1b$F0-simTraj_1b$F0.cov.adj, main="F0 error: cov of 173x and sampling indivs", xlab="Error in allele frequency (%)", ylab="Occurrences",
        xlim=c(-.3,.3), col="grey70")

hist(simTraj_1b$F25-simTraj_1b$F25.ind.adj, main="F25 error: sampling individuals", xlab="Error in allele frequency (%)", ylab="Occurrences",
        xlim=c(-.3,.3), col="grey70")

hist(simTraj_1b$F25-simTraj_1b$F25.cov.adj, main="F25 error: cov of 108x and sampling indivs", xlab="Error in allele frequency (%)", ylab="Occurrences",
        xlim=c(-.3,.3), col="grey70")

dev.off()
