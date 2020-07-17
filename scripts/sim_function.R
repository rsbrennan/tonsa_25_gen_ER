simulate_allele_freq <- function(startingAF, effective_pop, T0, T_N, samp_size, census_size, F0_depth, FN_depth){

  # simulate the trajectory
  simTraj_1 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_2 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_3 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_4 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)

# then need to add noise due to sampling only a subset of individuals as well as subset available DNA (ie, coverage)
## first, due to sampling subset of indivs:

  simTraj_1a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_1), dimnames=dimnames( simTraj_1))
  simTraj_2a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_2), dimnames=dimnames( simTraj_2))
  simTraj_3a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_3), dimnames=dimnames( simTraj_3))
  simTraj_4a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_4), dimnames=dimnames( simTraj_4))

  # introduce sampling variance to mimic Pool-seq of the entire population at coverage for each timepoint
  simTraj_1b <- cbind(simTraj_1,cbind((simTraj_1a),(sample.alleles(simTraj_1a[,1], size=F0_depth, mode="coverage"))))
  simTraj_1b <- cbind(simTraj_1b,(sample.alleles(simTraj_1a[,2], size=FN_depth, mode="coverage")))

  simTraj_2b <- cbind(simTraj_2,cbind((simTraj_2a),(sample.alleles(simTraj_2a[,1], size=F0_depth, mode="coverage"))))
  simTraj_2b <- cbind(simTraj_2b,(sample.alleles(simTraj_2a[,2], size=FN_depth, mode="coverage")))

  simTraj_3b <- cbind(simTraj_3,cbind((simTraj_3a),(sample.alleles(simTraj_3a[,1], size=F0_depth, mode="coverage"))))
  simTraj_3b <- cbind(simTraj_3b,(sample.alleles(simTraj_3a[,2], size=FN_depth, mode="coverage")))

  simTraj_4b <- cbind(simTraj_4,cbind((simTraj_4a),(sample.alleles(simTraj_4a[,1], size=F0_depth, mode="coverage"))))
  simTraj_4b <- cbind(simTraj_4b,(sample.alleles(simTraj_4a[,2], size=FN_depth, mode="coverage")))

  colnames(simTraj_1b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_2b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_3b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_4b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )

  # convert to sync format
  #first need to convert to read counts

  cnt_1_f0 <- round(simTraj_1b$F0.cov.adj * simTraj_1b$F0.cov.size, digits=0)
  cnt_1_f25 <- round(simTraj_1b$F25.cov.adj * simTraj_1b$F25.cov.size, digits=0)
  cov_1_f0 <- simTraj_1b$F0.cov.size
  cov_1_f25 <- simTraj_1b$F25.cov.size

  cnt_2_f0 <- round(simTraj_2b$F0.cov.adj * simTraj_2b$F0.cov.size, digits=0)
  cnt_2_f25 <- round(simTraj_2b$F25.cov.adj * simTraj_2b$F25.cov.size, digits=0)
  cov_2_f0 <- simTraj_2b$F0.cov.size
  cov_2_f25 <- simTraj_2b$F25.cov.size

  cnt_3_f0 <-  round(simTraj_3b$F0.cov.adj *  simTraj_3b$F0.cov.size, digits=0)
  cnt_3_f25 <- round(simTraj_3b$F25.cov.adj * simTraj_3b$F25.cov.size, digits=0)
  cov_3_f0 <-  simTraj_3b$F0.cov.size
  cov_3_f25 <- simTraj_3b$F25.cov.size

  cnt_4_f0 <-  round(simTraj_4b$F0.cov.adj *  simTraj_4b$F0.cov.size, digits=0)
  cnt_4_f25 <- round(simTraj_4b$F25.cov.adj * simTraj_4b$F25.cov.size, digits=0)
  cov_4_f0 <-  simTraj_4b$F0.cov.size
  cov_4_f25 <- simTraj_4b$F25.cov.size

  #concatenate counts and coverages for all replicates
  syncTable <- c()
  ListCnt <- list(cnt_1_f0,cnt_1_f25,
                  cnt_2_f0,cnt_2_f25,
                  cnt_3_f0,cnt_3_f25,
                  cnt_4_f0,cnt_4_f25)

  ListCov <- list(cov_1_f0,cov_1_f25,
                 cov_2_f0,cov_2_f25,
                 cov_3_f0,cov_3_f25,
                cov_4_f0,cov_4_f25)


  snpNum <- length(cnt_1_f0)

  f <- function(x){paste(ListCnt[[i]][x], 0,ListCov[[i]][x]-ListCnt[[i]][x] , 0, 0, 0, sep=":")}

  for (i in 1:length(ListCnt)) {
    syncTable <- cbind(syncTable, sapply(1:snpNum,f))
  }


  final_sync <- cbind(allSync@alleles$chr,  allSync@alleles$pos, c(rep('A', snpNum)),syncTable)


  return(final_sync)

}



############
## sim allele freqs for f25 comparisions
#############
simulate_allele_f25 <- function(startingAF, effective_pop, T0, T_N, samp_size, census_size, F0_depth, FN_depth){

  # simulate the trajectory
  simTraj_1 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_2 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_3 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_4 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_5 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_6 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_7 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)
  simTraj_8 <- wf.traj(p0=startingAF, Ne=effective_pop, t=c(0, 25), s=0, h=0.5, haploid=FALSE, approximate=FALSE)

# then need to add noise due to sampling only a subset of individuals as well as subset available DNA (ie, coverage)
## first, due to sampling subset of indivs:

  simTraj_1a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_1), dimnames=dimnames( simTraj_1))
  simTraj_2a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_2), dimnames=dimnames( simTraj_2))
  simTraj_3a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_3), dimnames=dimnames( simTraj_3))
  simTraj_4a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_4), dimnames=dimnames( simTraj_4))
  simTraj_5a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_5), dimnames=dimnames( simTraj_5))
  simTraj_6a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_6), dimnames=dimnames( simTraj_6))
  simTraj_7a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_7), dimnames=dimnames( simTraj_7))
  simTraj_8a <- matrix(sample.alleles(simTraj_1, size=samp_size, mode="individuals", Ncensus=census_size), nrow=nrow(simTraj_8), dimnames=dimnames( simTraj_8))


  # introduce sampling variance to mimic Pool-seq of the entire population at coverage for each timepoint
  simTraj_1b <- cbind(simTraj_1,cbind((simTraj_1a),(sample.alleles(simTraj_1a[,1], size=F0_depth, mode="coverage"))))
  simTraj_1b <- cbind(simTraj_1b,(sample.alleles(simTraj_1a[,2], size=FN_depth, mode="coverage")))

  simTraj_2b <- cbind(simTraj_2,cbind((simTraj_2a),(sample.alleles(simTraj_2a[,1], size=F0_depth, mode="coverage"))))
  simTraj_2b <- cbind(simTraj_2b,(sample.alleles(simTraj_2a[,2], size=FN_depth, mode="coverage")))

  simTraj_3b <- cbind(simTraj_3,cbind((simTraj_3a),(sample.alleles(simTraj_3a[,1], size=F0_depth, mode="coverage"))))
  simTraj_3b <- cbind(simTraj_3b,(sample.alleles(simTraj_3a[,2], size=FN_depth, mode="coverage")))

  simTraj_4b <- cbind(simTraj_4,cbind((simTraj_4a),(sample.alleles(simTraj_4a[,1], size=F0_depth, mode="coverage"))))
  simTraj_4b <- cbind(simTraj_4b,(sample.alleles(simTraj_4a[,2], size=FN_depth, mode="coverage")))

  simTraj_5b <- cbind(simTraj_5,cbind((simTraj_5a),(sample.alleles(simTraj_5a[,1], size=F0_depth, mode="coverage"))))
  simTraj_5b <- cbind(simTraj_5b,(sample.alleles(simTraj_5a[,2], size=FN_depth, mode="coverage")))

  simTraj_6b <- cbind(simTraj_6,cbind((simTraj_6a),(sample.alleles(simTraj_6a[,1], size=F0_depth, mode="coverage"))))
  simTraj_6b <- cbind(simTraj_6b,(sample.alleles(simTraj_6a[,2], size=FN_depth, mode="coverage")))

  simTraj_7b <- cbind(simTraj_7,cbind((simTraj_7a),(sample.alleles(simTraj_7a[,1], size=F0_depth, mode="coverage"))))
  simTraj_7b <- cbind(simTraj_7b,(sample.alleles(simTraj_7a[,2], size=FN_depth, mode="coverage")))

  simTraj_8b <- cbind(simTraj_8,cbind((simTraj_8a),(sample.alleles(simTraj_8a[,1], size=F0_depth, mode="coverage"))))
  simTraj_8b <- cbind(simTraj_8b,(sample.alleles(simTraj_8a[,2], size=FN_depth, mode="coverage")))

  colnames(simTraj_1b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_2b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_3b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_4b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_5b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_6b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_7b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )
  colnames(simTraj_8b) <- c("F0", "F25", "F0.ind.adj", "F25.ind.adj","F0.cov.adj", "F0.cov.size", "F25.cov.adj", "F25.cov.size" )

  # convert to sync format
  #first need to convert to read counts

  cnt_1_f0 <- round(simTraj_1b$F0.cov.adj * simTraj_1b$F0.cov.size, digits=0)
  cnt_1_f25 <- round(simTraj_1b$F25.cov.adj * simTraj_1b$F25.cov.size, digits=0)
  cov_1_f0 <- simTraj_1b$F0.cov.size
  cov_1_f25 <- simTraj_1b$F25.cov.size

  cnt_2_f0 <- round(simTraj_2b$F0.cov.adj * simTraj_2b$F0.cov.size, digits=0)
  cnt_2_f25 <- round(simTraj_2b$F25.cov.adj * simTraj_2b$F25.cov.size, digits=0)
  cov_2_f0 <- simTraj_2b$F0.cov.size
  cov_2_f25 <- simTraj_2b$F25.cov.size

  cnt_3_f0 <-  round(simTraj_3b$F0.cov.adj *  simTraj_3b$F0.cov.size, digits=0)
  cnt_3_f25 <- round(simTraj_3b$F25.cov.adj * simTraj_3b$F25.cov.size, digits=0)
  cov_3_f0 <-  simTraj_3b$F0.cov.size
  cov_3_f25 <- simTraj_3b$F25.cov.size

  cnt_4_f0 <-  round(simTraj_4b$F0.cov.adj *  simTraj_4b$F0.cov.size, digits=0)
  cnt_4_f25 <- round(simTraj_4b$F25.cov.adj * simTraj_4b$F25.cov.size, digits=0)
  cov_4_f0 <-  simTraj_4b$F0.cov.size
  cov_4_f25 <- simTraj_4b$F25.cov.size

  cnt_5_f0 <-  round(simTraj_5b$F0.cov.adj *  simTraj_5b$F0.cov.size, digits=0)
  cnt_5_f25 <- round(simTraj_5b$F25.cov.adj * simTraj_5b$F25.cov.size, digits=0)
  cov_5_f0 <-  simTraj_5b$F0.cov.size
  cov_5_f25 <- simTraj_5b$F25.cov.size

  cnt_6_f0 <-  round(simTraj_6b$F0.cov.adj *  simTraj_6b$F0.cov.size, digits=0)
  cnt_6_f25 <- round(simTraj_6b$F25.cov.adj * simTraj_6b$F25.cov.size, digits=0)
  cov_6_f0 <-  simTraj_6b$F0.cov.size
  cov_6_f25 <- simTraj_6b$F25.cov.size

  cnt_7_f0 <-  round(simTraj_7b$F0.cov.adj *  simTraj_7b$F0.cov.size, digits=0)
  cnt_7_f25 <- round(simTraj_7b$F25.cov.adj * simTraj_7b$F25.cov.size, digits=0)
  cov_7_f0 <-  simTraj_7b$F0.cov.size
  cov_7_f25 <- simTraj_7b$F25.cov.size

  cnt_8_f0 <-  round(simTraj_8b$F0.cov.adj *  simTraj_8b$F0.cov.size, digits=0)
  cnt_8_f25 <- round(simTraj_8b$F25.cov.adj * simTraj_8b$F25.cov.size, digits=0)
  cov_8_f0 <-  simTraj_8b$F0.cov.size
  cov_8_f25 <- simTraj_8b$F25.cov.size

  #concatenate counts and coverages for all replicates
  syncTable <- c()
  ListCnt <- list(cnt_1_f25,
                  cnt_2_f25,
                  cnt_3_f25,
                  cnt_4_f25,
                  cnt_5_f25,
                  cnt_6_f25,
                  cnt_7_f25,
                  cnt_8_f25)

  ListCov <- list(cov_1_f25,
                  cov_2_f25,
                  cov_3_f25,
                  cov_4_f25,
                  cov_5_f25,
                  cov_6_f25,
                  cov_7_f25,
                  cov_8_f25)


  snpNum <- length(cnt_1_f0)

  f <- function(x){paste(ListCnt[[i]][x], 0,ListCov[[i]][x]-ListCnt[[i]][x] , 0, 0, 0, sep=":")}

  for (i in 1:length(ListCnt)) {
    syncTable <- cbind(syncTable, sapply(1:snpNum,f))
  }


  final_sync <- cbind(allSync@alleles$chr,  allSync@alleles$pos, c(rep('A', snpNum)),syncTable)


  return(final_sync)
}
