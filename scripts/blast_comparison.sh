###########################################################################
## Script to run BLAST to compare:
##  a) Trinity transcriptome assembly using all reads
##  b) Trinity transcriptome assembly using digitally-normalized reads
##
## Author: John Stanton-Geddes
## Created: 2013-12-02
###########################################################################


# Make BLAST database for diginorm assembly
makeblastdb -dbtype nucl -in ../results/trinity/Trinity.fasta

# BLAST full assembly against diginorm assembly
# `-outfmt 6` to specify tabular format
# `-evalue 1e-20` to only keep significant hits
# `-max_target_seqs 1` to only keep top hit

blastn -query ../results/trinity-full/Trinity.fasta -db ../results/trinity/Trinity.fasta -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -out blast-full-vs-diginorm.txt

# Evaluate percent of diginorm transcript target length being aligned by best matching full transcript

/opt/software/trinityrnaseq_r2013_08_14/util/analyze_blastPlus_topHit_coverage.pl blast-full-vs-diginorm.txt ../results/trinity/Trinity.fasta ../results/trinity/Trinity.fasta

#### Reciprocal

makeblastdb -dbtype nucl -in ../results/trinity-full/Trinity.fasta

blastn -query ../results/trinity/Trinity.fasta -db ../results/trinity-full/Trinity.fasta -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -out blast-diginorm-vs-full.txt

/opt/software/trinityrnaseq_r2013_08_14/util/analyze_blastPlus_topHit_coverage.pl blast-diginorm-vs-full.txt ../results/trinity-full/Trinity.fasta ../results/trinity-full/Trinity.fasta
