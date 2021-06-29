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

### Data processing

- Check quality of sequence data: `fastqc.sh`
- Trim the minor adapter presence from samples: `trim.sh`
- re-check quality to make sure trimming was good: `fastqc_posttrim.sh`

### Align, call snps

- align with bwa mem: `align.sh`
- alignment stats: `align_stats.sh`
- merge resulting bams from each lane/sample: `merge_bams.sh`
- plot resulting stats: `plot_alignment_stats.R`

### drift simulations



### cvtk




