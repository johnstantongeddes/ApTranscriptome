*Aphaenogaster* transcriptome assembly
===================================

This directory contains scripts used for assembly of the *Aphaenogaster* transcriptome.

The raw Illumina short-read data is in the local directory `/home/data/Aphaeno_transcriptome/130509_SN1073_0326_BD25DAACXX/Project_Stanton-Geddes_Project_001`.

General workflow is to

1) Cut adapters, quality trim reads and generate QC reports using [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
2) Perform ['digital normalization'](http://ged.msu.edu/papers/2012-diginorm/) of trimmed reads
3) Merge paired-reads using [FLASH](http://bioinformatics.oxfordjournals.org/content/early/2011/09/07/bioinformatics.btr507)
4) Assemble using [velvet-oases](http://www.ebi.ac.uk/~zerbino/oases/))

Intermediate files and assembly are stored locally in /data directory as these files are large and easily recreated from scripts so should not be tracked by git. 

/results directory will include summary statistics and plots of assembly from various methods.

TODO
a) in silico mRNA spike-in to validate assembly
