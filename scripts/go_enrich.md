# assign go terms to each gene



## SNPs 

~/tonsa_genomics/analysis/cmh_results.txt

~/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt'

```r

library(dplyr)

snp <- read.csv("~/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt", header=T, sep="\t")
afs <- read.csv("~/tonsa_genomics/analysis/filtered_allele_freqs.txt", header=T, sep="\t")
cmh <- read.csv("~/tonsa_genomics/analysis/cmh.fdr.txt", header=T, sep="\t")

cmh$class <- snp$class
cmh$annotation <- snp$annotation
cmh$distance <- snp$distance

# drop the lab significant loci from the list:
cmh$ah_fdr[cmh$ah_sig_nolab == TRUE] <- NA
cmh$ha_fdr[cmh$ha_sig_nolab == TRUE] <- NA
cmh$hh_fdr[cmh$hh_sig_nolab == TRUE] <- NA

# make empty list to hold dataframes that we'll fill in
snp_out <- vector(mode = "list", length = 4)
names(snp_out) <- c("all", "genic", "exon", "promoter")

class_all <- c("promoter", "exon", "downstream", "intron")
class_ids <- list(class_all, c("exon", "intron"), c("exon"), c("promoter"))
names(class_ids) <- c("all", "genic", "exon", "promoter")

for(class_index in 1:length(class_ids)){
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



use the simulated af changes



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

#dedup_df <- as.data.frame(matrix(ncol=ncol(annot), nrow=0))

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

# no lab sig:
#all <- dat %>% filter(annotation != "-")
genes <- unique(all$annotation)
toTest <- c("aa_sig","ah_sig_nolab", "ha_sig_nolab", "hh_sig_nolab" )

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
        all_filt <- allRes[which(allRes$weight < 0.05 ),]
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

hh <- read.csv("~/tonsa_genomics/analysis/go_enrich/hh_sig_nolab_CMH_GO.txt", header=TRUE, sep="\t")
ha <- read.csv("~/tonsa_genomics/analysis/go_enrich/ha_sig_nolab_CMH_GO.txt", header=TRUE, sep="\t")
ah <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_sig_nolab_CMH_GO.txt", header=TRUE, sep="\t")
aa <- read.csv("~/tonsa_genomics/analysis/go_enrich/ah_sig_CMH_GO.txt", header=TRUE, sep="\t")

hh$group <- "Greenhouse"
ha$group <- "Warm"
ah$group <- "Acidic"
aa$group <- "Ambient"

dat <- rbind(rbind(hh, ha), ah)

# change p-value threshold for plotting

dat <- dat[(dat$weight < 0.05),]
nrow(dat)

#BiocManager::install("org.Dm.eg.db")
### IMPORTANT!
# trace(godata, edit=TRUE)
## need to edit this to work with topgo
### on line 12: columns = c("GO", "ONTOLOGY")))
#### change GO to GOALL
## note that this makes the go mapping data much slower to compute.
## this needs to be done every time youopen a new R session

hsGObp <- godata('org.Dm.eg.db', ont="BP", computeIC=FALSE)
hsGOmf <- godata('org.Dm.eg.db', ont="MF", computeIC=FALSE)
hsGOcc <- godata('org.Dm.eg.db', ont="CC",computeIC=FALSE)

# pull out the go terms from the df
gobp <- unique(as.character(dat$GO.ID[which(dat$ontology == "BP")]))
gocc <- unique(as.character(dat$GO.ID[which(dat$ontology == "CC")]))
gomf <- unique(as.character(dat$GO.ID[which(dat$ontology == "MF")]))

# calculate distance matrix
bpdist <- (mgoSim(gobp, gobp, semData=hsGObp, measure="Wang", combine=NULL))
ccdist <- (mgoSim(gocc, gocc, semData=hsGOcc, measure="Wang", combine=NULL))
mfdist <- (mgoSim(gomf, gomf, semData=hsGOmf, measure="Wang", combine=NULL))
# using lapply, loop over columns and match values to the look up table. store in "new".
# this gives us the actual names of the GO terms
row.names(bpdist) <- paste(row.names(bpdist),(unlist(lapply(row.names(bpdist), function(x) dat$Term[match(x, dat$GO.ID)]))), sep=": ")
colnames(bpdist) <- paste(colnames(bpdist),(unlist(lapply(colnames(bpdist), function(x) dat$Term[match(x, dat$GO.ID)]))), sep=": ")

row.names(ccdist) <- paste(row.names(ccdist),(unlist(lapply(row.names(ccdist), function(x) dat$Term[match(x, dat$GO.ID)]))), sep=": ")
colnames(ccdist) <- paste(colnames(ccdist),(unlist(lapply(colnames(ccdist), function(x) dat$Term[match(x, dat$GO.ID)]))), sep=": ")

row.names(mfdist) <- paste(row.names(mfdist),(unlist(lapply(row.names(mfdist), function(x) dat$Term[match(x, dat$GO.ID)]))), sep=": ")
colnames(mfdist) <- paste(colnames(mfdist),(unlist(lapply(colnames(mfdist), function(x) dat$Term[match(x, dat$GO.ID)]))), sep=": ")

clusterbp <- hclust(1-as.dist(bpdist), method = "ward.D2")
clustercc <- hclust(1-as.dist(ccdist), method = "ward.D2")
clustermf <- hclust(1-as.dist(mfdist), method = "ward.D2")

#convert cluster object to use with ggplot
dendrbp <- dendro_data(clusterbp, type="rectangle") 
dendrcc <- dendro_data(clustercc, type="rectangle") 
dendrmf <- dendro_data(clustermf, type="rectangle") 

#your own labels (now rownames) are supplied in geom_text() and label=label
bpplot <- ggplot() + 
  geom_segment(data=segment(dendrbp), aes(x=x, y=y, xend=xend, yend=yend)) + 
  #geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=4) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
 scale_x_continuous(breaks = label(dendrbp)$x,
                    labels=label(dendrbp)$label,
                     position = "top")+
   theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = rel(1.1), hjust=10),
        panel.grid = element_blank(),
         axis.text.x=element_blank(),
)
ccplot <- ggplot() + 
  geom_segment(data=segment(dendrcc), aes(x=x, y=y, xend=xend, yend=yend)) + 
  #geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=4) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
 scale_x_continuous(breaks = label(dendrcc)$x,
                    labels=label(dendrcc)$label,
                     position = "top")+
   theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = rel(1.1), hjust=10),
        panel.grid = element_blank(),
         axis.text.x=element_blank(),
)

mfplot <- ggplot() + 
  geom_segment(data=segment(dendrmf), aes(x=x, y=y, xend=xend, yend=yend)) + 
  #geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=4) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
 scale_x_continuous(breaks = label(dendrmf)$x,
                    labels=label(dendrmf)$label,
                     position = "top")+
   theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = rel(1.1), hjust=10),
        panel.grid = element_blank(),
         axis.text.x=element_blank(),
)

ggsave("~/tonsa_genomics/figures/topgo_enrich_bp_cluster.pdf",bpplot, h=11, w=7)
ggsave("~/tonsa_genomics/figures/topgo_enrich_cc_cluster.pdf",ccplot, h=3.6, w=7)
ggsave("~/tonsa_genomics/figures/topgo_enrich_mf_cluster.pdf",mfplot, h=4, w=7)


# make plot in same order as dendrogram

dat$group <- factor(dat$group, levels = c("Acidic","Greenhouse", "Warm"))

bp_dat <- dat[dat$ontology == "BP",]
cc_dat <- dat[dat$ontology == "CC",]
mf_dat <- dat[dat$ontology == "MF",]

bp_labs <- substring(label(dendrbp)$label,13)
cc_labs <- substring(label(dendrcc)$label,13)
mf_labs <- substring(label(dendrmf)$label,13)

dat$Term <- factor(dat$Term, levels = c(bp_labs,cc_labs,mf_labs))

bp_dat$Term <- factor(bp_dat$Term, levels = bp_labs)
cc_dat$Term <- factor(cc_dat$Term, levels = cc_labs)
mf_dat$Term <- factor(mf_dat$Term, levels = mf_labs)

p1 <- ggplot(dat, aes(x=group, y=factor(Term), group=Term, fill=-log10(weight), size=Significant)) + 
    geom_line(size=0.5, col="grey45")+
  geom_point(pch=21, color="black") +
  facet_grid(rows=vars(ontology), scales = "free", space = "free") +
  scale_fill_gradient(low = "grey", high = "firebrick3") +
  theme_bw() +
    scale_size(name   = "Number of Significant\nGenes",
             breaks = c(10, 50, 100, 250, 500)) +
    labs(fill="-log10(p)") + ylab("")

ggsave("~/tonsa_genomics/figures/topgo_enrich.pdf",p1, h=15, w=9)


# indiv plots

pbp <- ggplot(bp_dat, aes(x=group, y=factor(Term), group=Term, fill=-log10(weight), size=Significant)) + 
    geom_line(size=0.5, col="grey45")+
  geom_point(pch=21, color="black") +
  #facet_grid(rows=vars(ontology), scales = "free", space = "free") +
  scale_fill_gradient(low = "grey", high = "firebrick3") +
  theme_bw() +
    scale_size(name   = "Number of Significant\nGenes",
             breaks = c(10, 50, 100, 250, 500)) +
    labs(fill="-log10(p)") + ylab("") +
    theme(axis.text=element_text(size=6))

pmf <- ggplot(mf_dat, aes(x=group, y=factor(Term), group=Term, fill=-log10(weight), size=Significant)) + 
    geom_line(size=0.5, col="grey45")+
  geom_point(pch=21, color="black") +
  #facet_grid(rows=vars(ontology), scales = "free", space = "free") +
  scale_fill_gradient(low = "grey", high = "firebrick3") +
  theme_bw() +
    scale_size(name   = "Number of Significant\nGenes",
             breaks = c(10, 50, 100, 250, 500)) +
    labs(fill="-log10(p)") + ylab("")

pcc <- ggplot(cc_dat, aes(x=group, y=factor(Term), group=Term, fill=-log10(weight), size=Significant)) + 
    geom_line(size=0.5, col="grey45")+
  geom_point(pch=21, color="black") +
  #facet_grid(rows=vars(ontology), scales = "free", space = "free") +
  scale_fill_gradient(low = "grey", high = "firebrick3") +
  theme_bw() +
    scale_size(name   = "Number of Significant\nGenes",
             breaks = c(10, 50, 100, 250, 500)) +
    labs(fill="-log10(p)") + ylab("")

ggsave("~/tonsa_genomics/figures/topgo_enrich_bp.pdf",pbp, h=9, w=6)
ggsave("~/tonsa_genomics/figures/topgo_enrich_mf.pdf",pmf, h=6, w=7)
ggsave("~/tonsa_genomics/figures/topgo_enrich_cc.pdf",pcc, h=6, w=7)

```
