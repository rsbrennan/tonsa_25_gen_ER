
import pandas as pd
import numpy as np
import os

# files:
snp_file = os.path.expanduser("~/tonsa_genomics/analysis/gene_level_analysis/annotation_table1.txt")
annotation = "/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz"

############
# link each snp to its category
############

# read in snp file
snp_in = pd.read_csv(snp_file, delimiter="\t",header='infer')
# subset to only the first 3 cols

# read in gene file
gene_in = pd.read_csv(annotation, delimiter="\t",header=None, skiprows=2)


gene_in['gene'] = [x[:-3] for x in gene_in[0]]

snp_in['description'] = 'NA'

for idx, row in snp_in.iterrows():
  # split out gene from the cmh table
  tmp_annot = row['annotation'].split(";")
  if len(tmp_annot) == 1 and tmp_annot[0] == "-":
    snp_in.at[idx,'description'] = "NA"
  if len(tmp_annot) > 1:
    out_desc = ""
    for mtchrow in tmp_annot:
      long_id = gene_in[gene_in['gene'].str.contains(mtchrow)][8].str.cat(sep=";  ")
      tmp_desc = [i for i in long_id.split(";") if i.startswith('Description')]
      tmp_uniq = ";".join(set(tmp_desc[0].split(","))) # remove duplicates, join together with string.
      tmp_list = [out_desc, tmp_uniq]
      out_desc = ";".join(filter(None,tmp_list))
      snp_in.at[idx,'description'] = out_desc
  if len(tmp_annot) == 1 and tmp_annot[0] != "-":
    long_id = gene_in[gene_in['gene'].str.contains(tmp_annot[0])][8].str.cat(sep=";  ")
    tmp_desc = [i for i in long_id.split(";") if i.startswith('Description')]
    tmp_uniq = ";".join(set(tmp_desc[0].split(","))) # remove duplicates, join together with string.
    snp_in.at[idx,'description'] = tmp_uniq
  if idx % 5000 == 0:
    print(idx)


snp_in.to_csv('/users/r/b/rbrennan/tonsa_genomics/analysis/results_with_gene_ids.txt', index=False, sep='\t')

