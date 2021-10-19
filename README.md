# A. tonsa 25 gen evolve and resequence

This repository hold scripts used to analyse the genomic data from a 25 generation selection experiment of Acartia tonsa. 

## Data availability

Made public once preprint is posted

Note that there were 2 lanes of sequencing. Data processing code below is provided for only one lane, but the methods were identical. Data were merged after aligning.

## plotting details

For my own sanity, here are the colors I'm using for each treatment:

- Founding population: `#D3DDDC`, pch= 21  
- Ambient: `#6699CC`, pch= 21  
- Acidic: `#F2AD00`, pch=22  
- Warm: `#00A08A`, pch= 23  
- Greenhouse: `#CC3333`, pch=24 


## description of scripts

Below here is a description of scripts used in the analyses.


### Genome annotation




### Data processing

- Check quality of sequence data: `fastqc.sh`
- Trim the minor adapter presence from samples: `trim.sh`
- re-check quality to make sure trimming was good: `fastqc_posttrim.sh`

### Align, call snps

- align with bwa mem: `align.sh`
- alignment stats: `align_stats.sh`
- merge resulting bams from each lane/sample: `merge_bams.sh`
- plot resulting stats: `plot_alignment_stats.R`
- Call snps with varscan: `varscan.sh`
- Filter raw SNPs for depth, missingnesss, etc: `filter_raw_snps.R`

### identify candidate SNPs

- convert af file to sync file: `to_sync.py`
- calculate cmh from sync file above: `cmh.sh`
- annotate all snps: `annotate.py`

Simulate drift expectations for significance corrections
- simulation function: `sim_function.R`
- simulate data: `sim_af.R`
- calc cmh: `~/tonsa_genomics/scripts/sim_af.sh`
- run actual sims: `sim_af.sh`
- then take these sims and the cmh results from above, and combine into table `cmh_results.txt`: `cmh_process.R`

and finally, `calc_significance.R` outputs `cmh.fdr.txt`


### cvtk


### Figures

- Fig. 1: 
- Fig. 2: 
- Fig. 3: 
- Fig. 4:
- Fig. 5:
- Fig. 6:

Supplemental figures:
- Fig. S1:
- Fig. S2:
- Fig. S3:
- Fig. S4:
- Fig. S5:
- Fig. S6:
- Fig. S7:
- Fig. S8:


