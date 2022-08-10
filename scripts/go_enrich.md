

# assign go terms to each gene



## SNPs

~/tonsa_genomics/analysis/cmh_results.txt
~/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt

```r

library(dplyr)

snp <- read.csv("~/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt", header=T, sep="\t")
afs <- read.csv("~/tonsa_genomics/analysis/filtered_allele_freqs.txt", header=T, sep="\t")
cmh <- read.csv("~/tonsa_genomics/analysis/cmh.fdr.txt", header=T, sep="\t")

cmh$class <- snp$class
cmh$annotation <- snp$annotation
cmh$distance <- snp$distance

# drop the lab significant loci from the list:
#cmh$ah_fdr[cmh$ah_sig_nolab == TRUE] <- NA
#cmh$ha_fdr[cmh$ha_sig_nolab == TRUE] <- NA
#cmh$hh_fdr[cmh$hh_sig_nolab == TRUE] <- NA

# make empty list to hold dataframes that we'll fill in
snp_out <- vector(mode = "list", length = 4)
names(snp_out) <- c("all", "genic", "exon", "promoter")

class_all <- c("promoter", "exon", "downstream", "intron")
class_ids <- list(class_all, c("exon", "intron"), c("exon"), c("promoter"))
names(class_ids) <- c("all", "genic", "exon", "promoter")

for(class_index in 1){
    # filter to get only the loci in the class of interest.
    tmp_filt <- cmh %>% filter(class %in% as.vector(class_ids[[class_index]]) )
    # get the unique genes from the filtered dataset we just made
    genes_filt <- (unique(unlist(strsplit(as.character(tmp_filt$annotation),split = ";"))))
    # make the dataframe that we need to fill in
    tmp_df <- as.data.frame(matrix(nrow=length(genes_filt), ncol=6))
    colnames(tmp_df) <- c("gene", "n_snp",
                        "aa_pval_min","ah_pval_min","ha_pval_min","hh_pval_min"
                        )
    # loop over all the genes in this df, save mean, number snps, etc.
    for(i in 1:length(genes_filt)){
        # pull down loci that match each gene
        tmp_gene <- tmp_filt[grep(genes_filt[i], tmp_filt$annotation),]
        # add gene to output df
        tmp_df$gene[i] <- genes_filt[i]

        # save number of loci
        tmp_df$n_snp[i] <- nrow(tmp_gene)

        # let's also take the min p value for each set, so I can parse by significant etc, later.
        tmp_df$aa_pval_min[i] <- mean(tmp_gene$aa_fdr, na.rm=TRUE)
        tmp_df$ah_pval_min[i] <- mean(tmp_gene$ah_fdr, na.rm=TRUE)
        tmp_df$ha_pval_min[i] <- mean(tmp_gene$ha_fdr, na.rm=TRUE)
        tmp_df$hh_pval_min[i] <- mean(tmp_gene$hh_fdr, na.rm=TRUE)

        if(i %% 500 == 0){print(i)}
    }

    # add to the list of dataframes
    snp_out[[names(class_ids)[class_index]]] <- tmp_df
    print("done with")
    print(names(class_ids)[class_index])
}


#load(file = "~/tonsa_genomics/analysis/gene_level_analysis/gene_level_summary.RData")


# save genes to pull down all possible gene identities (make gene "world")
genes_snp <- (unique(unlist(strsplit(as.character(snp$annotation),split = ";"))))

snp_rm <- genes_snp[which(genes_snp != "-")]
write.table(file="~/tonsa_genomics/analysis/go_enrich/snp_genes.txt", snp_rm,sep="\t", row.names=F, quote=FALSE, col.names=F)

# log transform, then save each to csv

for(i in 1:length(snp_out)){

    aa_sv <- as.data.frame(cbind(snp_out[[i]]$gene, -log10(snp_out[[i]]$aa_pval_min)))
    ah_sv <- as.data.frame(cbind(snp_out[[i]]$gene, -log10(snp_out[[i]]$ah_pval_min)))
    ha_sv <- as.data.frame(cbind(snp_out[[i]]$gene, -log10(snp_out[[i]]$ha_pval_min)))
    hh_sv <- as.data.frame(cbind(snp_out[[i]]$gene, -log10(snp_out[[i]]$hh_pval_min)))


    aa_sv$V2[which(hh_sv$V2 == -Inf)] <- 0
    ah_sv$V2[which(ah_sv$V2 == -Inf)] <- 0
    ha_sv$V2[which(ha_sv$V2 == -Inf)] <- 0
    hh_sv$V2[which(hh_sv$V2 == -Inf)] <- 0

    write.table(file=paste("~/tonsa_genomics/analysis/go_enrich/snp_min_aa_",
                            names(snp_out)[i],
                            "_goenrich.txt", sep = ""),
                aa_sv, col.names=TRUE,
                row.names=FALSE, quote=FALSE,sep=",")

    write.table(file=paste("~/tonsa_genomics/analysis/go_enrich/snp_min_ah_",
                            names(snp_out)[i],
                            "_goenrich.txt", sep = ""),
                ah_sv, col.names=TRUE,
                row.names=FALSE, quote=FALSE,sep=",")

    write.table(file=paste("~/tonsa_genomics/analysis/go_enrich/snp_min_ha_",
                            names(snp_out)[i],
                            "_goenrich.txt", sep = ""),
                ha_sv, col.names=TRUE,
                row.names=FALSE, quote=FALSE,sep=",")

    write.table(file=paste("~/tonsa_genomics/analysis/go_enrich/snp_min_hh_",
                            names(snp_out)[i],
                            "_goenrich.txt", sep = ""),
                hh_sv, col.names=TRUE,
                row.names=FALSE, quote=FALSE,sep=",")
}



```


assign go terms to final gene sets

```python

import pandas as pd
import numpy as np
import gzip
import csv
import re
from collections import OrderedDict

# assign GO terms
filepath = '/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz'


print("Starting go assignments")

# gene assignment is in col 17

#make empty array
num_lines = sum(1 for line in open('/users/r/b/rbrennan/tonsa_genomics/analysis/go_enrich/snp_genes.txt'))

go_df = np.empty(shape=(num_lines,2), dtype = object)

ict=0 # start counter

with open('/users/r/b/rbrennan/tonsa_genomics/analysis/go_enrich/snp_genes.txt') as master_file:
        for idx, line in enumerate(master_file):
            tmp_gene = line.split("\n")[0]
            tmp_gene2 = tmp_gene.split(";")
            out_go = []
            # match GO terms
            match = []
            for entry in tmp_gene2:
                with gzip.open(filepath, mode="rt") as file:
                    for gene_line in file:
                        if len(re.findall(entry, gene_line)) > 0:
                            gene_spl = gene_line.split("\n")[0].split("\t")
                            go_spl = [i for i in gene_spl[8].split(";") if i.startswith('Ontology_id')]
                            if len(go_spl) > 0:
                                if len(go_spl[0].split()) > 1:
                                    go_tmp = go_spl[0].split()[1]
                                if len(out_go) == 0:
                                    out_go = go_tmp
                                if len(out_go) > 0:
                                    out_go = out_go + "," + go_tmp
            go_df[idx,0] = tmp_gene
            if len(out_go) > 0:
                out_go = ";".join(set(out_go.split(","))) # remove duplicates, join together with string.
                go_df[idx,1] = out_go
                #go_df[idx,1] = out_go.replace('"',"").replace(',',";")
            if len(out_go) == 0:
                go_df[idx,1] = "unknown"
            ict=ict+1
            if ict % 1000 == 0: print(ict)

np.savetxt('/users/r/b/rbrennan/tonsa_genomics/analysis/go_enrich/cmh_GOterms.out', go_df,fmt='%s', delimiter='\t')

```



```bash
# there are quotes in these files, remove them. Could do it earlier, but.. I haven't.


cat ~/tonsa_genomics/analysis/go_enrich/cmh_GOterms.out | sed 's/"//g' | grep -v 'unknown' > ~/tonsa_genomics/analysis/go_enrich/cmh_GOterms.corrected.out


```


## topGO

need tab delim table with gene id in col1, and comma delim go terms in col 2.

convert to commas:

```bash

cat ~/tonsa_genomics/analysis/go_enrich/cmh_GOterms.corrected.out | sed 's/;/,/g' > ~/tonsa_genomics/analysis/go_enrich/cmh_GOterms.topgo.out
```



```r
library(topGO)
library(dplyr)


# find af change:
df_fdr <- read.csv("~/tonsa_genomics/analysis/cmh.fdr.txt", sep="\t", header=T)

geneID2GO <- readMappings(file = "~/tonsa_genomics/analysis/go_enrich/cmh_GOterms.topgo.out")

dat <- read.csv("~/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt", header=TRUE, sep="\t")
annot <- data.frame(SNP = dat$SNP, class=dat$class, annotation=dat$annotation, distance = dat$distance)

goterms <- read.csv("~/tonsa_genomics/analysis/go_enrich/cmh_GOterms.topgo.out", header=FALSE, sep="\t")
colnames(goterms) <- c("annotation", "goterm")


# need to split out genes with more than one annotation. They need to be two rows.
unlist_genes <- strsplit(as.character(annot$annotation),split = ";")

length(which(lengths(unlist_genes) > 1))
length(which(lengths(unlist_genes) == 1))


dup_genes <- which(lengths(unlist_genes) > 1)

#num rows?
num_row <- sum(lengths(unlist_genes)[dup_genes])

dedup_df <- as.data.frame(matrix(ncol=ncol(annot), nrow=num_row))
colnames(dedup_df) <- colnames(annot)
dedup_df$class <- as.character(dedup_df$class)
dedup_df$SNP <- as.character(dedup_df$SNP)

ict <- 0
for(i in dup_genes){
    tmp_df <- as.data.frame(matrix(nrow=length(unlist_genes[[i]]), ncol=ncol(annot)))

    # add vals to this:
    tmp_df[] <- annot[i,]
    tmp_df$V3 <- unlist_genes[[i]]
    tmp_df[] <- lapply(tmp_df, as.character)
    #bind to new df
    #starting row:
    start_row <- min(which(is.na(dedup_df$annot)))
    dedup_df[start_row:(start_row+nrow(tmp_df)-1),] <- tmp_df
    ict <- ict +1
    if(ict%%5000 == 0){print(ict)}
}

nondup_genes <- which(lengths(unlist_genes) == 1)

all_dedup <- rbind(dedup_df,annot[nondup_genes,])

# need to parse down dat to only those with go terms.
all <- inner_join(all_dedup, goterms, "annotation")
all2 <- inner_join(all, df_fdr, "SNP")

all <- all2
nrow(all)
# 92926

all$aa_all <- FALSE
all$ah_all <- FALSE
all$ha_all <- FALSE
all$hh_all <- FALSE
all$aa_all[which(all$aa_fdr < 0.01) ] <- TRUE
all$ah_all[which(all$ah_fdr < 0.01) ] <- TRUE
all$ha_all[which(all$ha_fdr < 0.01) ] <- TRUE
all$hh_all[which(all$hh_fdr < 0.01) ] <- TRUE

all$aa_sig <- FALSE
all$ah_sig <- FALSE
all$ha_sig <- FALSE
all$hh_sig <- FALSE
all$aa_sig[which(all$aa_fdr < 0.01 & all$ah_fdr >= 0.01 & all$ha_fdr >= 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$ah_sig[which(all$ah_fdr < 0.01 & all$aa_fdr >= 0.01 & all$ha_fdr >= 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$ha_sig[which(all$ha_fdr < 0.01 & all$ah_fdr >= 0.01 & all$aa_fdr >= 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$hh_sig[which(all$hh_fdr < 0.01 & all$ah_fdr >= 0.01 & all$ha_fdr >= 0.01 & all$aa_fdr >= 0.01)] <- TRUE

# two-way overlaps
all$aa_ah <- FALSE
all$aa_ha <- FALSE
all$aa_hh <- FALSE
all$ah_ha <- FALSE
all$ah_hh <- FALSE
all$ha_hh <- FALSE

all$aa_ah[which(all$aa_fdr < 0.01 & all$ah_fdr < 0.01 & all$ha_fdr >= 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$aa_ha[which(all$aa_fdr < 0.01 & all$ha_fdr < 0.01 & all$ah_fdr >= 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$aa_hh[which(all$aa_fdr < 0.01 & all$hh_fdr < 0.01 & all$ha_fdr >= 0.01 & all$ah_fdr >= 0.01)] <- TRUE
all$ah_ha[which(all$ah_fdr < 0.01 & all$ha_fdr < 0.01 & all$aa_fdr >= 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$ah_hh[which(all$ah_fdr < 0.01 & all$hh_fdr < 0.01 & all$ha_fdr >= 0.01 & all$aa_fdr >= 0.01)] <- TRUE
all$ha_hh[which(all$ha_fdr < 0.01 & all$hh_fdr < 0.01 & all$aa_fdr >= 0.01 & all$ah_fdr >= 0.01)] <- TRUE

# 3-way overlaps
all$aa_ah_ha <- FALSE
all$aa_ah_hh <- FALSE
all$aa_ha_hh <- FALSE
all$ah_ha_hh <- FALSE

all$aa_ah_ha[which(all$aa_fdr < 0.01 & all$ah_fdr < 0.01 & all$ha_fdr < 0.01 & all$hh_fdr >= 0.01)] <- TRUE
all$aa_ah_hh[which(all$aa_fdr < 0.01 & all$ah_fdr < 0.01 & all$hh_fdr < 0.01 & all$ha_fdr >= 0.01)] <- TRUE
all$aa_ha_hh[which(all$aa_fdr < 0.01 & all$ha_fdr < 0.01 & all$hh_fdr < 0.01 & all$ah_fdr >= 0.01)] <- TRUE
all$ah_ha_hh[which(all$ha_fdr < 0.01 & all$ah_fdr < 0.01 & all$hh_fdr < 0.01 & all$aa_fdr >= 0.01)] <- TRUE

# all overlaps
all$aa_ah_ha_hh   <- FALSE
all$aa_ah_ha_hh[which(all$aa_fdr < 0.01 & all$ah_fdr < 0.01 & all$ha_fdr < 0.01 & all$hh_fdr < 0.01)] <- TRUE


table(all$aa_all)
table(all$ah_all)
table(all$ha_all)
table(all$hh_all)
table(all$aa_sig)
table(all$ah_sig)
table(all$ha_sig)
table(all$hh_sig)
table(all$aa_hh)
table(all$aa_ah_ha)
table(all$aa_ah_hh)
table(all$aa_ah)
table(all$aa_ha)
table(all$aa_hh)
table(all$ah_ha)
table(all$ah_hh)
table(all$ha_hh)
table(all$aa_ah_ha_hh)


#all <- dat %>% filter(annotation != "-")
genes <- unique(all$annotation)
toTest <- c("aa_all","ah_all", "ha_all", "hh_all",
            "aa_sig","ah_sig", "ha_sig", "hh_sig",
            "aa_ah","aa_ha","aa_hh","ah_ha","ah_hh","ha_hh",
            "aa_ah_ha","aa_ah_hh","aa_ha_hh","ah_ha_hh",
            "aa_ah_ha_hh")


for(test_set in toTest){
    out <- as.data.frame(matrix(nrow=length(genes), ncol=2))
    ofInterest <- c()
    out.save <- as.data.frame(matrix(ncol=8, nrow=0))
    colnames(out.save) <- c("GO.ID", "Term","Annotated","Significant","Expected","classicFisher","weight", "ontology")
    for(i in 1:length(genes)){
        a <- all[which(all$annotation == genes[i]),]
        out[i,1] <- as.character(a$annotation[1])
        out[i,2] <- as.character(a$goterm[1])
        ofInterest[i] <- ifelse(sum(a[,test_set]) > 0, TRUE, FALSE)
    }

    print(test_set)
    print(sum(ofInterest))
    print(length(ofInterest))

   write.table(file= "~/tonsa_genomics/analysis/go_enrich/test.go", out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

    system('cat ~/tonsa_genomics/analysis/go_enrich/test.go |  sed \'s/;/,/g\' | sed \'s/,$//g\' > ~/tonsa_genomics/analysis/go_enrich/test.annotation')


geneID2GO <- readMappings(file = "~/tonsa_genomics/analysis/go_enrich/test.annotation")
geneID2GO <- geneID2GO[2:length(geneID2GO)]
    # set gene background
    geneUniverse <- names(geneID2GO)

    genesOfInterest <- names(geneID2GO[ofInterest])

    #show genes of interest in universe vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse

 for(j in c("BP", "CC", "MF")){

        myGOdata <- (new("topGOdata", description="My project",
            ontology=j, allGenes=geneList,
            annot = annFUN.gene2GO, gene2GO = geneID2GO,
            nodeSize = 5 ))
        resultWeight <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
        resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
        allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight"  ,
            topNodes = length(resultWeight@score))
        allRes$ontology <- j
        all_filt <- allRes[which(allRes$weight < 0.05 & allRes$Significant > 1 ),]
        if(nrow(all_filt) > 0){


            out.save <- rbind(out.save,all_filt)

            #### get genes in go terms

            myterms =all_filt$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
           mygenes = genesInTerm(myGOdata, myterms)

          genes_out <- as.data.frame(matrix(nrow=0, ncol=2))
          colnames(genes_out) <- c("goterm", "genes")
          for (i in 1:length(myterms)){
             myterm=myterms[i]
             mygenesforterm= mygenes[myterm][[1]]
             mygenesforterm=paste(mygenesforterm, collapse=',')
             tmp_term <- data.frame(goterm = myterm, genes = mygenesforterm)
             genes_out <- rbind(genes_out, tmp_term)
           }
        }
    write.table(genes_out,
        paste("~/tonsa_genomics/analysis/go_enrich/",test_set, "_", j, "_genetoGOmapping.txt", sep=""),
            sep="\t",quote=F, row.names=F)

    }

    write.table(file=paste("~/tonsa_genomics/analysis/go_enrich/", test_set, "_CMH_GO.txt", sep=""),
             out.save, col.names=TRUE,
            row.names=FALSE, quote=FALSE,sep="\t")

}



```


### plot results

```r

library(GOSemSim)
library(ggplot2)
library(ggdendro)
library(dplyr)
library(org.Dm.eg.db)
library(ggpubr)
library(UpSetR)


hh_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/hh_all_CMH_GO.txt", header=TRUE, sep="\t")
ha_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/ha_all_CMH_GO.txt", header=TRUE, sep="\t")
ah_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_all_CMH_GO.txt", header=TRUE, sep="\t")
aa_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_all_CMH_GO.txt", header=TRUE, sep="\t")

hh_all_bp <- hh_all$GO.ID[hh_all$ontology == "BP"]
ha_all_bp <- ha_all$GO.ID[ha_all$ontology == "BP"]
ah_all_bp <- ah_all$GO.ID[ah_all$ontology == "BP"]
aa_all_bp <- aa_all$GO.ID[aa_all$ontology == "BP"]

aa_all[(aa_all_bp %in% ah_all_bp),]
aa_all[(aa_all_bp %in% ha_all_bp),]
aa_all[(aa_all_bp %in% hh_all_bp),]
ah_all[(ah_all_bp %in% ha_all_bp),]
ah_all[(ah_all_bp %in% hh_all_bp),]
ha_all[(ha_all_bp %in% hh_all_bp),]
Reduce(intersect, list(aa_all_bp,ah_all_bp,ha_all_bp))
Reduce(intersect, list(aa_all_bp,ah_all_bp,hh_all_bp))
Reduce(intersect, list(aa_all_bp,ha_all_bp,hh_all_bp))
Reduce(intersect, list(ah_all_bp,ha_all_bp,hh_all_bp))
Reduce(intersect, list(aa_all_bp,ah_all_bp,ha_all_bp,hh_all_bp))

all <- length(Reduce(intersect, list(aa_all_bp,ah_all_bp,ha_all_bp,hh_all_bp)))


listInput <- list(Ambient = aa_all_bp,
                  Acidification = ah_all_bp,
                  Warming = ha_all_bp,
                  OWA = hh_all_bp)


### note that this figure in not in the supplemental. 

pdf("~/tonsa_genomics/figures/all_upSet_GOTerms_overlaps.pdf",
   height = 4, width = 5)
upset(fromList(listInput),keep.order = TRUE, empty.intersections = "on",
        sets=c("OWA","Warming","Acidification","Ambient"),
        mainbar.y.label = "Significant GO enrichment in intersecting SNP set",
        sets.x.label = "Total number of significan GO terms",
        point.size = 3.4, line.size = 1.2, 
          sets.bar.color=rev(c('#6699CC',"#F2AD00", "#00A08A", "#CC3333")),
           text.scale = c(1.1, 1.1, 1, 1, 1.5, 1))
dev.off()


#######
# go enrichment for overlapping snp sets.


hh_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/hh_all_CMH_GO.txt", header=TRUE, sep="\t")
ha_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/ha_all_CMH_GO.txt", header=TRUE, sep="\t")
ah_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_all_CMH_GO.txt", header=TRUE, sep="\t")
aa_all <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_all_CMH_GO.txt", header=TRUE, sep="\t")

hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/hh_sig_CMH_GO.txt", header=TRUE, sep="\t")
ha <- read.csv("~/tonsa_genomics/analysis/go_enrich/ha_sig_CMH_GO.txt", header=TRUE, sep="\t")
ah <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_sig_CMH_GO.txt", header=TRUE, sep="\t")
aa <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_sig_CMH_GO.txt", header=TRUE, sep="\t")

aa_ah <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_ah_CMH_GO.txt", header=TRUE, sep="\t")
aa_ha <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_ha_CMH_GO.txt", header=TRUE, sep="\t")
aa_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_hh_CMH_GO.txt", header=TRUE, sep="\t")
ah_ha <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_ha_CMH_GO.txt", header=TRUE, sep="\t")
ah_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_hh_CMH_GO.txt", header=TRUE, sep="\t")
ha_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/ha_hh_CMH_GO.txt", header=TRUE, sep="\t")

aa_ah_ha <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_ah_ha_CMH_GO.txt", header=TRUE, sep="\t")
aa_ah_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_ah_hh_CMH_GO.txt", header=TRUE, sep="\t")
aa_ha_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_ha_hh_CMH_GO.txt", header=TRUE, sep="\t")
ah_ha_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_ha_hh_CMH_GO.txt", header=TRUE, sep="\t")

aa_ah_ha_hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/aa_ah_ha_hh_CMH_GO.txt", header=TRUE, sep="\t")


tstNm <-  c("aa_all","ah_all", "ha_all", "hh_all",
            "aa","ah", "ha", "hh",
            "aa_ah","aa_ha","aa_hh","ah_ha","ah_hh","ha_hh",
            "aa_ah_ha","aa_ah_hh","aa_ha_hh","ah_ha_hh",
            "aa_ah_ha_hh")

toTest <- list(aa_all,ah_all, ha_all, hh_all,
            aa,ah, ha, hh,
            aa_ah,aa_ha,aa_hh,ah_ha,ah_hh,ha_hh,
            aa_ah_ha,aa_ah_hh,aa_ha_hh,ah_ha_hh,
            aa_ah_ha_hh)

names(toTest) <- tstNm
bp.out <- list()

for (i in 1:length(toTest)){
    bp.out[[i]] <- toTest[[i]]$GO.ID[toTest[[i]]$ontology == "BP"]
    names(bp.out)[i] <- names(toTest)[i]
}

aa_sig <- length(bp.out[['aa']])
ah_sig <- length(bp.out[['ah']])
ha_sig <- length(bp.out[['ha']])
hh_sig <- length(bp.out[['hh']])

aa_ah_sig <- length(bp.out[['aa_ah']])
aa_ha_sig <- length(bp.out[['aa_ha']])
aa_hh_sig <- length(bp.out[['aa_hh']])
ah_ha_sig <- length(bp.out[['ah_ha']])
ah_hh_sig <- length(bp.out[['ah_hh']])
ha_hh_sig <- length(bp.out[['ha_hh']])

aa_ah_hh_sig <- length(bp.out[['aa_ah_hh']])
aa_ha_hh_sig <- length(bp.out[['aa_ha_hh']])
aa_ah_ha_sig <- length(bp.out[['aa_ah_ha']])
ah_ha_hh_sig <- length(bp.out[['ah_ha_hh']])

aa_ah_ha_hh_sig <- length(bp.out[['aa_ah_ha_hh']])


# find overlap counts
library(UpSetR)
all_input <- all_overlaps <- c(
                        "Ambient" = aa_sig,
                        "Acidic" = ah_sig,
                        "Warm" = ha_sig,
                        "OWA" = hh_sig,
                "Ambient&Acidic" = aa_ah_sig,
                "Ambient&Warm" = aa_ha_sig,
                "Ambient&OWA" = aa_hh_sig,
                "Acidic&Warm" = ah_ha_sig,
                "Acidic&OWA" = ah_hh_sig,
                "Warm&OWA" = ha_hh_sig,
                "Ambient&Acidic&OWA" = aa_ah_hh_sig,
                "Ambient&Warm&OWA" = aa_ha_hh_sig,
                "Ambient&Acidic&Warm" = aa_ah_ha_sig,
                "Acidic&Warm&OWA" = ah_ha_hh_sig,
                "Ambient&Acidic&Warm&OWA" = aa_ah_ha_hh_sig)




pdf("~/tonsa_genomics/figures/all_upSet_GOTerms.pdf",
    height = 4, width = 5)

upset(fromExpression(all_input),
        #order.by = "degree", 
        #group.by = "sets",
        keep.order = TRUE, empty.intersections = "on",
        sets = c("OWA","Warm", "Acidic","Ambient"),
        mainbar.y.label = "Significant GO enrichment in intersecting SNP set",
        sets.x.label = "Total number of significan GO terms",
        point.size = 3.4, line.size = 1.2, 
          sets.bar.color=rev(c('#6699CC',"#F2AD00", "#00A08A", "#CC3333")),
           text.scale = c(1.1, 1.1, 1, 1, 1.5, 1)
)

dev.off()


```
