# genome annotation

## GAWN genome annotation

found here: https://github.com/enormandeau/gawn

./gawn 02_infos/gawn_config.sh

Need to download the blast database. `/data/blastDB`

To avoid redundancy, should use collapsed isoforms. In this case I can use supertranscripts, which includes all exons in the gene. See trinity github page: https://github.com/trinityrnaseq/trinityrnaseq/tree/18bb3182005a61d236525185d8dda65867ca4cc0/Analysis/SuperTranscripts

I have done this previously:

/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/atonsa_super_transcript.fasta

remember:

Exons = gene - introns

CDS = gene - introns - UTRs

CDS = Exons - UTRs

remember that gmap does try to identify UTR regions and leave CDS


#### Run GAWN

My GAWN config file looks like this:

```
#!/bin/bash

# Modify the following parameter values according to your experiment
# Do not modify the parameter names or remove parameters
# Do not add spaces around the equal (=) sign

# Global parameters
NCPUS=10                    # Number of CPUs to use for analyses (int, 1+)

# Genome indexing
SKIP_GENOME_INDEXING=1      # 1 to skip genome indexing, 0 to index it

# Genome annotation with transcriptome
# NOTE: do not use compressed fasta files
GENOME_NAME="GCA_900241095.1_Aton1.0_genomic.fa"                  # Name of genome fasta file found in 03_data
TRANSCRIPTOME_NAME="atonsa_super_transcript.fasta"    # Name of transcriptome fasta file found in 03_data

# Swissprot
SWISSPROT_DB="/data/programs/blastDB/protein_db/swissprot"

```


```bash


#update blast databases
update_blastdb.pl --decompress swissprot

md5sum -c *.md5

## run GAWN

# remember to turn off genome indexing if already done.
cd ~/tonsa_annotation/gawn

nohup ~/tonsa_annotation/gawn/gawn 02_infos/gawn_config.sh 2> ~/tonsa_annotation/gawn/log_out/gawn.stderr_$(date +"%F_%R").txt 1> ~/tonsa_annotation/gawn/log_out/gawn.stdout_$(date +"%F_%R").txt &
echo $! > ~/tonsa_annotation/gawn/log_out/gawn.pid


```

46,045 aligned genes out of 72511 in the FASTA

Remember that each gene can be split into exons, etc., so total lines in alignment is 69243

how many have a swiss prot annotation?

```bash
cat ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.tsv| tail -n +2 | awk '$7 != "-"' | cut -f 5 |  sort | uniq | wc -l
# 15249

# no annotation
cat ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.tsv | tail -n +2 | awk '$7 == "-"' | cut -f 5 |  sort | uniq | wc -l
#30796
```

15,249 have annotation, 30,796 with no annotation


## Use transcriptome annotations to update the gff, etc.

These are from the following paper:

```

Jørgensen, T.S., Petersen, B., Petersen, H.C.B., Browne, P.D., Prost, S., Stillman, J.H., Hansen, L.H. and Hansen, B.W., 2019. The genome and mRNA transcriptome of the cosmopolitan calanoid copepod acartia tonsa dana improve the understanding of copepod genome size evolution. Genome biology and evolution, 11(5), pp.1440-1450.

```

`/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz`

`~/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.ALLTRANSCRIPTS.gff3`

Update the GAWN gff with improved annotations.  



```python

import gzip
import pandas as pd
import numpy as np
import re

gff_loc = "/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.gff.gz"

gff_in = pd.read_csv(gff_loc, delimiter="\t", skiprows= 2,header=None)
gff_in.columns = ['Transcript_id','source','feature','start','end','score','strand','frame','attributes']
# add gene id
gff_in['gene_id'] = gff_in['Transcript_id'].str.rsplit("_",1, expand=True)[0]

genes_loc = "/data/copepods/tonsa_transcriptome/Atonsa_gen_trans_agp_gff/blast2go_gff_export_20170613_0815.GENES_ONLY.txt"
num_gff_lines = sum(1 for line in open(genes_loc))

dict_out = np.empty(shape=(num_gff_lines,9), dtype = object)

with open(genes_loc, "rt") as genes_loc:
    for genes_idx, genes_line in enumerate(genes_loc):
        ts_id = genes_line.split("\n")[0]
        ln_match = gff_in[gff_in['gene_id'].str.match(ts_id)]
        gff_spl = ln_match.head(1)['attributes'].str.split(";").tolist()[0]
        # separate out the attributes:
        # find gene name"
        indices = [i for i, s in enumerate(gff_spl) if 'Description' in s]
        if len(indices) > 0:
            gene_name = gff_spl[1].replace('Description ','').replace('"','')
        else:
            gene_name = '-'
        # gene alt names
        indices = [i for i, s in enumerate(gff_spl) if 'Gene' in s]
        if len(indices) > 0:
            gene_alt_names = '; '.join([gff_spl[i] for i in indices]).replace('Gene ','').replace('"','')
        else:
            gene_alt_names = '-'
        #Gene go
        indices = [i for i, s in enumerate(gff_spl) if 'Ontology_id' in s]
        if len(indices) > 0:
            go_term = gff_spl[indices[0]].replace('Ontology_id ','').replace('"','')
        else:
            go_term = '-'
        # Cellular component
        indices = [i for i, s in enumerate(gff_spl) if 'Ontology_term' in s]
        if len(indices) > 0:
            go_cellular = gff_spl[indices[0]].replace('Ontology_term ','').replace('"','')
        else:
            go_cellular = '-'
        # gene id
        dict_out[genes_idx,0] = ts_id
        # gene accession
        dict_out[genes_idx,1] = '-'
        # gene name
        dict_out[genes_idx,2] = gene_name
        # gene alt names
        dict_out[genes_idx,3] = gene_alt_names
        # pfam
        dict_out[genes_idx,4] = '-'
        #Gene go
        dict_out[genes_idx,5] = go_term
        # Cellular component
        dict_out[genes_idx,6] = go_cellular
        #Mol component
        dict_out[genes_idx,7] = '-'
        # Bio process
        dict_out[genes_idx,8] = '-'
        # replace the NA in the gene name
        na_ind = [i for i, s in enumerate(gff_spl) if '---NA---' in s]
        if len(na_ind) > 0:
            dict_out[genes_idx,2] = '-'
        # if there are multiople matches, need to merge them
        if len(ln_match) > 1:
            for index, row in ln_match.iloc[1:].iterrows():
                gff_spl = row['attributes'].split(";")
                indices = [i for i, s in enumerate(gff_spl) if 'Description' in s]
                if len(indices) > 0:
                    gene_name = gff_spl[1].replace('Description ','').replace('"','')
                else:
                    gene_name = '-'
                # gene alt names
                indices = [i for i, s in enumerate(gff_spl) if 'Gene' in s]
                if len(indices) > 0:
                    gene_alt_names = '; '.join([gff_spl[i] for i in indices]).replace('Gene ','').replace('"','')
                else:
                    gene_alt_names = '-'
                #Gene go
                indices = [i for i, s in enumerate(gff_spl) if 'Ontology_id' in s]
                if len(indices) > 0:
                    go_term = gff_spl[indices[0]].replace('Ontology_id ','').replace('"','')
                else:
                    go_term = '-'
                # Cellular component
                indices = [i for i, s in enumerate(gff_spl) if 'Ontology_term' in s]
                if len(indices) > 0:
                    go_cellular = gff_spl[indices[0]].replace('Ontology_term ','').replace('"','')
                else:
                    go_cellular = '-'
                # gene id
                tmp_nm = '; '.join(set([dict_out[genes_idx,2] + "; " + gene_name][0].split("; ")))
                dict_out[genes_idx,2] = tmp_nm
                # alt names
                tmp_alt = '; '.join(set([dict_out[genes_idx,3]+ "; " + gene_alt_names][0].split("; ")))
                dict_out[genes_idx,3] = tmp_alt
                # go terms
                if go_term != "-":
                    tmp_go = '; '.join(set([dict_out[genes_idx,5]+ "; " + go_term][0].split("; ")))
                    dict_out[genes_idx,5] = tmp_go
                # go cellular
                if go_cellular != "-":
                    tmp_cell = '; '.join(set([dict_out[genes_idx,6]+ "; " + go_cellular][0].split("; ")))
                    dict_out[genes_idx,6] = tmp_cell
        if genes_idx % 1000 == 0:
            print(genes_idx)

############
# use dictionary to fill in annotation table:
############

df1 = pd.DataFrame(dict_out)
df1.columns = ['Transcript_id','GeneAccession','GeneName','GeneAltNames','GenePfam','GeneGo','CellularComponent','MolecularFunction','BiologicalProcess']

gff_dict = df1.set_index('Transcript_id').to_dict(orient="index")

# make empty array to fill:

annotation_file = '/users/r/b/rbrennan/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.tsv'

num_lines = sum(1 for line in open(annotation_file))-1

first_line = next(open(annotation_file)).split("\n")[0]
num_cols = len(first_line.split("\t"))

df_out = np.empty(shape=(num_lines,num_cols), dtype = object)

ct_not_miss = 0
ct_new_annot = 0

with open(annotation_file) as annotation:
    h0 = next(annotation)
    for idx, line in enumerate(annotation):
        line = line.split("\n")[0]
        # check if gene name is missing. if it is, all other columns are missing
        if len(line.split("\t")[7]) > 1:
            df_out[idx,] = line.split("\t")
            ct_not_miss = ct_not_miss + 1
        else:
            df_out[idx,0:6] = line.split("\t")[0:6]
            # look up values to add with dictionary
            tmp_values = gff_dict[line.split("\t")[4]]
            # check if new name is empty:
            #if tmp_values['GeneName'] != "-" :
            df_out[idx,first_line.split("\t").index("GeneName")] = tmp_values['GeneName']
            # check if gene alt name is empty:
            df_out[idx,first_line.split("\t").index("GeneAltNames")] = tmp_values['GeneAltNames']
            df_out[idx,first_line.split("\t").index("GeneGo")] = tmp_values['GeneGo']
            df_out[idx,first_line.split("\t").index("GeneGo")] = "-"
            # pfam is always empty if replacing
            df_out[idx,first_line.split("\t").index("GenePfam")] = "-"
            # the rest just fill in with key
            df_out[idx,first_line.split("\t").index("GeneAccession")] = tmp_values['CellularComponent']
            df_out[idx,first_line.split("\t").index("CellularComponent")] = tmp_values['CellularComponent']
            df_out[idx,first_line.split("\t").index("MolecularFunction")] = tmp_values['MolecularFunction']
            df_out[idx,first_line.split("\t").index("BiologicalProcess")] = tmp_values['BiologicalProcess']
            #df_out[idx,first_line.split("\t").index("GeneAltNames")] = tmp_values['GeneAltNames']
            ct_new_annot = ct_new_annot + 1

ct_not_miss
# 24268
ct_new_annot
# 44974

df_sv = np.vstack((first_line.split("\t"),df_out))

np.savetxt('/users/r/b/rbrennan/tonsa_annotation/gawn/annotation_merged_2020_05_27.txt', df_sv,fmt='%s', delimiter='\t')


```

count up how the annotation looks. A little more sophisticated than the numbers above.

``` python

import pandas as pd
import numpy as np
import os

# files: 
dat_file = os.path.expanduser("~/tonsa_annotation/gawn/annotation_merged_2020_05_27.txt")
############
# link each snp to its category
############

# read in snp file
dfn = pd.read_csv(dat_file, delimiter="\t",header='infer')

# how many of the transcripts have annotations. ie, not counting split mappings, etc.

dfn.columns.values

len(dfn.TranscriptName.unique())
# 46045 unique genes

u_trans = dfn.TranscriptName.unique()

ct_no_annot = 0
ct_annot = 0

for count, trans in enumerate(u_trans,1):
    tmp_df = trans
    t_loc = dfn.index[dfn['TranscriptName'] == trans].tolist()
    t_rows = dfn.loc[[t_loc[0]]]
    if t_rows['GeneName'].values[0] is "-":
        ct_no_annot = ct_no_annot + 1
    if t_rows['GeneName'].values[0] is not "-":
        ct_annot = ct_annot + 1
    if count % 5000 == 0:
        print(count)

ct_no_annot
# 24957
ct_annot
# 21088
count
# 46045

```

somewhere I missed the na... just do it here:

```bash

sed -i 's/---NA---;/-/g' /users/r/b/rbrennan/tonsa_annotation/gawn/annotation_merged_2020_05_27.txt
sed -i 's/- -/-/g' /users/r/b/rbrennan/tonsa_annotation/gawn/annotation_merged_2020_05_27.txt 

```

Adding the annotations from the transcriptome got us another 5k or so annotations.


## What lineage is this reference genome?


```bash

/data/programs/ncbi-blast/makeblastdb -in /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa -parse_seqids -dbtype nucl -out reference_blast_db


cd /data/copepods/tonsa_genome

blastn -db /data/copepods/tonsa_genome/reference_blast_db -query ~/tonsa_genomics/analysis/mitotypes/coI_F_clade.fasta -outfmt '6 sseqid sseq'


```
