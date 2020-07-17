cd ~/tonsa_genomics/data/aligned

ls *.bam > bam.file.in

~/bin/angsd/angsd -bam bam.file.in -minMapQ 30 -minQ 20 -GL 2  -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -minMaf 0.025 -P 5
