
import pandas as pd
import numpy as np
import os

# files:
snp_file = os.path.expanduser("~/tonsa_genomics/analysis/cmh_results.txt")
annotation_file = os.path.expanduser("~/tonsa_annotation/gawn/annotation_merged_2020_05_27.txt")
gene_file = os.path.expanduser("~/tonsa_genomics/analysis/gene_level_analysis/snp.genes.bed")
exon_file = os.path.expanduser("~/tonsa_genomics/analysis/gene_level_analysis/snp.exons.bed")

############
# link each snp to its category
############

# read in snp file
snp_in = pd.read_csv(snp_file, delimiter="\t",header='infer')
# subset to only the first 3 cols
snp_ids = snp_in.iloc[:, range(3)]

# read in gene file
gene_in = pd.read_csv(gene_file, delimiter="\t",header=None)
# read in snp file
exon_in = pd.read_csv(exon_file, delimiter="\t",header=None)

gene_in['SNP'] = gene_in[0].astype(str) + ":" + gene_in[2].astype(str)
exon_in['SNP'] = exon_in[0].astype(str) + ":" + exon_in[2].astype(str)

# split out name:
gene_in['annotation'] = gene_in[11].str.split("Name=",n = 1, expand = True)[1]
exon_in['annotation'] = exon_in[11].str.split("Name=",n = 1, expand = True)[1]
gene_in.annotation.fillna(value="-", inplace=True)
exon_in.annotation.fillna(value="-", inplace=True)

gene_in[5] = gene_in[5].replace(".", "-")
exon_in[5] = exon_in[5].replace(".", "-")

snp_in["class"] = str(np.nan)
snp_in["annotation"] = str(np.nan)
snp_in["distance"] = str(np.nan)

idx_dbl = []
idx_exon_dbl = []
idx_zero_dbl = []
idx_intron_dbl = []
for idx, row in snp_in.iterrows():
  gene_match = gene_in[gene_in['SNP'] == (row['SNP'])].reset_index(drop=True)
  exon_match = exon_in[exon_in['SNP'] == (row['SNP'])].reset_index(drop=True)
  #### check for not matches.
  if (exon_match.empty and gene_match.empty):
    #print("no match!")
    snp_in.at[idx,'class'] = "-"
    snp_in.at[idx,'annotation'] = "-"
    snp_in.at[idx,'distance'] = "-"
  #### check if there's more than 1 match. multiple genes
  if len(gene_match.index) > 1 and len(exon_match.index) > 1:
    idx_dbl.append(idx)
  # if both match and not missing see which is closer. or if value is 0, take exon
    ### if both equal to zero, take exon. dif exon in one, intron in another. call it a exon.
    if any(x == 0 for x in exon_match[12]):
      snp_in.at[idx,'class'] = "exon"
      # find which are 0
      idx_zero = exon_match.index[exon_match[12] == 0].tolist()
      for ex_row in idx_zero:
        snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + exon_match['annotation'][ex_row].split(";")[0]
      snp_in.at[idx,'distance'] = 0
      # remove the nan
      snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
      idx_exon_dbl.append(idx)
    #### check for 0 in gene, and no 0 in exon these are introns. assing to both genes
    if all(x != 0 for x in exon_match[12]) and any(x == 0 for x in gene_match[12]):
      idx_intron_dbl.append(idx)
      snp_in.at[idx,'class'] = "intron"
      # find the annotation that is 0
      idx_zero = gene_match.index[gene_match[12] == 0].tolist()
      for ex_row in idx_zero:
        snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + gene_match['annotation'][ex_row].split(";")[0]
      snp_in.at[idx,'distance'] = 0
      # remove the nan
      snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
    # check if no zeros, then find lowest.
    # In this case, it can't be an exon. must be either up or downstream
    if all(x  != 0 for x in exon_match[12]) and all(x != 0 for x in gene_match[12]):
      #print("none are zero")
      idx_zero_dbl.append(idx)
      snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[gene_match[12].abs().idxmin()].split(";")[0]
      snp_in.at[idx,'distance'] = gene_match[12].iloc[gene_match[12].abs().idxmin()]
      if snp_in.at[idx,'distance'] < 0:
        snp_in.at[idx,'class'] = "promoter"
      if snp_in.at[idx,'distance'] > 0:
        snp_in.at[idx,'class'] = "downstream"
  #### check if one exon, but two genes:
  if len(gene_match.index) > 1 and len(exon_match.index) == 1:
    # if exon is 0, choose this one.
    if any(x == 0 for x in exon_match[12]):
      snp_in.at[idx,'class'] = "exon"
      idx_zero = exon_match.index[exon_match[12] == 0].tolist()
      snp_in.at[idx,'annotation'] = exon_match['annotation'].iloc[0].split(";")[0]
      snp_in.at[idx,'distance'] = 0
    #  check for 0 in gene, and no 0 in exon these are introns. assing to both genes
    if all(x != 0 for x in exon_match[12]) and any(x == 0 for x in gene_match[12]):
      #print("exon not 0, gene is 0")
      idx_intron_dbl.append(idx)
      snp_in.at[idx,'class'] = "intron"
      # find the annotation that is 0
      idx_zero = gene_match.index[gene_match[12] == 0].tolist()
      for ex_row in idx_zero:
        snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + gene_match['annotation'][ex_row].split(";")[0]
      snp_in.at[idx,'distance'] = 0
      # remove the nan
      snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
    # check if no zeros, then find lowest.
    # In this case, it can't be an exon. must be either up or downstream
    if all(x  != 0 for x in exon_match[12]) and all(x != 0 for x in gene_match[12]):
      idx_zero_dbl.append(idx)
      snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[gene_match[12].abs().idxmin()].split(";")[0]
      snp_in.at[idx,'distance'] = gene_match[12].iloc[gene_match[12].abs().idxmin()]
      if snp_in.at[idx,'distance'] < 0:
        snp_in.at[idx,'class'] = "promoter"
      if snp_in.at[idx,'distance'] > 0:
        snp_in.at[idx,'class'] = "downstream"
  # check if two exons, but one genes. only a few hundred of these:
  if len(gene_match.index) == 1 and len(exon_match.index) > 1:
    # if both equal to zero, take exon. if exon in one, intron in another. call it a exon.
    if any(x == 0 for x in exon_match[12]):
      snp_in.at[idx,'class'] = "exon"
      idx_zero = exon_match.index[exon_match[12] == 0].tolist()
      for ex_row in idx_zero:
        snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + exon_match['annotation'][ex_row].split(";")[0]
      snp_in.at[idx,'annotation'] = ";".join(pd.unique(snp_in.at[idx,'annotation'].split(";")))
      snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
      snp_in.at[idx,'distance'] = 0
    # first check for 0 in gene, and no 0 in exon these are introns. assing to both genes
    if all(x != 0 for x in exon_match[12]) and any(x == 0 for x in gene_match[12]):
      #print("exon not 0, gene is 0")
      idx_intron_dbl.append(idx)
      snp_in.at[idx,'class'] = "intron"
      # find the annotation that is 0
      idx_zero = gene_match.index[gene_match[12] == 0].tolist()
      snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
      snp_in.at[idx,'distance'] = 0
    # check if no zeros, then find lowest.
    # In this case, it can't be an exon. must be either up or downstream
    if all(x  != 0 for x in exon_match[12]) and all(x != 0 for x in gene_match[12]):
      snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0]
      snp_in.at[idx,'distance'] = gene_match[12].iloc[0]
      if snp_in.at[idx,'distance'] < 0:
        snp_in.at[idx,'class'] = "promoter"
      if snp_in.at[idx,'distance'] > 0:
        snp_in.at[idx,'class'] = "downstream"
  #if both match, do they have any actual annotation?
  if len(gene_match.index) == 1 and len(exon_match.index) == 1:
    if (exon_match['annotation'].iloc[0] == "-" and gene_match['annotation'].iloc[0] == "-"):
      snp_in.at[idx,'class'] = "-"
      snp_in.at[idx,'annotation'] = "-"
      snp_in.at[idx,'distance'] = "-"
  # if both match and not missing take exon if its 0.
  # otherwise, take intron, up, or down stream
    if (exon_match['annotation'].iloc[0] != "-" and gene_match['annotation'].iloc[0] != "-"):
      # if both equal to zero, take exon
      if (exon_match[12].iloc[0] == 0 and gene_match[12].iloc[0] == 0 ):
        snp_in.at[idx,'class'] = "exon"
        snp_in.at[idx,'annotation'] = exon_match['annotation'].iloc[0].split(";")[0]
        snp_in.at[idx,'distance'] = 0
      # if one is zero, take that one.
      else:
        # first check for 0 in one
        if (exon_match[12].iloc[0] == 0):
          snp_in.at[idx,'class'] = "exon"
          snp_in.at[idx,'annotation'] = exon_match['annotation'].iloc[0].split(";")[0]
          snp_in.at[idx,'distance'] = 0
        if (exon_match[12].iloc[0] != 0 and gene_match[12].iloc[0] == 0 ):
          snp_in.at[idx,'class'] = "intron"
          snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
          snp_in.at[idx,'distance'] = 0
        # check if upstream
        # if negative,
        if (gene_match[12].iloc[0] < 0 ):
          snp_in.at[idx,'class'] = "promoter"
          snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
          snp_in.at[idx,'distance'] = gene_match[12].iloc[0]
        # check if downstream
        if (gene_match[12].iloc[0] > 0 ):
          snp_in.at[idx,'class'] = "downstream"
          snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
          snp_in.at[idx,'distance'] = gene_match[12].iloc[0]
  if idx % 10000 == 0:
    print(idx)

# running at 11:15 6/28
# screen -r 34591

snp_in.to_csv('/users/r/b/rbrennan/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt', index=False, sep='\t')


snp_in['class'].value_counts()

