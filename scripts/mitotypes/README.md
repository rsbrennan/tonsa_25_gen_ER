# pipeline to analyze mitotype data
 
This is to be run after all of the steps in `mitotypes.md`, in this same directory. See this file for instructions.

This pipeline holds an automated method to analyze mitotype data for a tonsa based on the phylogeny of Figueroa et al., 2020: Phylogeography ofAcartia tonsaDana, 1849 (Calanoida: Copepoda)and phylogenetic reconstruction of the genusAcartiaDana, 1846.

View the help file with: `./mitotype_pipeline.sh --help` or `./mitotype_pipeline.sh -h`

```bash
./mitotype_pipeline.sh your.fasta plot_name

```

where `your.fasta` contains your sequence data and `plot_name` is the name you want for your output plot.

All output will be directed to `output/`, which will be created if it doesn't already exist. 

I assume some R packages are installed: `treeio`, `tidyverse`, `ggtree`, `ape`, `stringr`.

All of this could be run wherever, but paths to mr bayes and muscle are hard coded right now, so you'd need to edit that if exporting this elsewhere.
