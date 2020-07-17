#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
# Contig_0000008    582 .   A   G   .   PASS    .   GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR0/0:99:118:118:118:0:0%:1E0:36:0:112:6:0:0 0/1:15:157:157:152:5:3.18%:3.0256E-2:35:35:148:4:5:0    0/0:99:125:125:125:0:0%:1E0:35:0:121:4:0:0   0/0:99:87:87:86:1:1.15%:5E-1:34:35:81:5:1:0 0/0:99:179:179:179:0:0%:1E0:34:0:174:5:0:0 0/0:99:128:128:128:0:0%:1E0:35:0:128:0:0:0

#Sync file:

#2R  2302    N   0:7:0:0:0:0 0:7:0:0:0:0
#2R  2303    N   0:8:0:0:0:0 0:8:0:0:0:0
#2R  2304    N   0:0:9:0:0:0 0:0:9:0:0:0
#2R  2305    N   1:0:9:0:0:0 0:0:9:1:0:0

#col1: reference contig
#col2: position within the refernce contig
#col3: reference base
#col4: allele frequencies of population number 1
#col5: allele frequencies of population number 2
#coln: allele frequencies of population number n

#The allele frequencies are in the format A:T:C:G:N:del, i.e: count of bases 'A', count of bases 'T',... and deletion count in the end (character '*' in the mpileup)


# for varscan file:
#       Each sample has six values separated by colons:
#           Cons - consensus genotype in IUPAC format
#           Cov - total depth of coverage
#           Reads1 - number of reads supporting reference
#           Reads2 - number of reads supporting variant
#           Freq - the variant allele frequency by read count
#           P-value - FET p-value of observed reads vs expected non-variant


library(stringr)

dat <- read.table("~/tonsa_genomics/analysis/filtered_variants.txt", header=TRUE, stringsAsFactors=FALSE)

pops <- colnames(dat)[11:ncol(dat)]


syncOut <- as.data.frame(matrix(nrow=nrow(dat), ncol=length(pops) + 3))
colnames(syncOut) <- c("Chr", "Pos", "Ref", pops)

# remember A:T:C:G:N:del

for (i in 1:nrow(dat)){

    syncOut[i, 1:3] <-  (dat[i, 1:3])

    for ( pop in 1:length(pops)){

        pop_in <- pops[pop]

        if(dat$Ref[i] == "A"){

            if(dat$Var[i] == "T"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste(pop_split[3], pop_split[4], "0", "0", "0", "0", sep=":")
            }

            if(dat$Var[i] == "C"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste(pop_split[3], "0",pop_split[4], "0", "0", "0", sep=":")
            }

            if(dat$Var[i] == "G"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste(pop_split[3], "0", "0", pop_split[4], "0", "0", sep=":")
            }

        }

        if(dat$Ref[i] == "T"){

            if(dat$Var[i] == "A"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste(pop_split[4], pop_split[3], "0", "0", "0", "0", sep=":")
            }

            if(dat$Var[i] == "C"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste("0", pop_split[3], pop_split[4], "0", "0", "0", sep=":")
            }

            if(dat$Var[i] == "G"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste("0", pop_split[3] , "0", pop_split[4], "0", "0", sep=":")
            }

            }

        if(dat$Ref[i] == "C"){

            if(dat$Var[i] == "A"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste(pop_split[4], "0", pop_split[3] , "0", "0", "0", sep=":")
            }

            if(dat$Var[i] == "T"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste("0", pop_split[4], pop_split[3],  "0", "0", "0", sep=":")
            }

            if(dat$Var[i] == "G"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste("0", "0", pop_split[3] , pop_split[4], "0", "0", sep=":")
            }

            }

        if(dat$Ref[i] == "G"){

            if(dat$Var[i] == "A"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste(pop_split[4], "0","0", pop_split[3],  "0", "0", sep=":")
            }

            if(dat$Var[i] == "T"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste("0", pop_split[4], "0", pop_split[3] , "0", "0", sep=":")
            }

            if(dat$Var[i] == "C"){

                pop_split <- str_split_fixed(dat[i,pop_in], ":", n=6)
                syncOut[i,pop_in] <- paste("0", "0", pop_split[4], pop_split[3] , "0", "0", sep=":")
            }

            }
    }

if(i%%10000 == 0){ print(i)}
}

write.table(syncOut, file="~/tonsa_genomics/analysis/variants.sync", col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
