# Script to run BUSCO analysis for ApTranscriptome
# John Stanton-Geddes
# 7 January 2016 

# Install BUSCO and dependencies (NCBI Blast, HMMER, EMBOSS
# following directions at 
# http://busco.ezlab.org/files/README.html
# Installed BUSCO to the /scripts directory

# Installed Blast and HMMER to ~/software directory

# NCBI: downloaded binaries from http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# HMMER 3.1b2: downloaded binaries from http://hmmer.janelia.org/
# Installed EMBOSS using homebrew 
# https://www.biostars.org/p/107092/
# brew install homebrew/science/emboss

# Add binaries to PATH

PATH=$PATH:/Users/JSG/software/ncbi-blast-2.3.0+/bin/:/Users/JSG/software/hmmer-3.1b2-macosx-intel/binaries/:/usr/local/Cellar/emboss/6.6.0/bin

# change directory
cd ../results/trinity-full

# rename FASTA headers for BUSCO input
awk '/^>/{print ">chromosome" ++i; next}{print}' < Trinity_cap3_uclust_clean.fa > ../../scripts/BUSCO_v1.1b1/short_headers.fa

# run BUSCO
cd ../../scripts/BUSCO_v1.1b1
python3 BUSCO_v1.1b1.py -o Aph20160107 -in short_headers.fa -l arthropoda -m trans
