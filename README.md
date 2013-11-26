*Aphaenogaster* transcriptome assembly
=======================================

*Forked* from [github.com/johnstantongeddes/climate-cascade/projects/ApTranscriptome]

Initial development of assembly in that repository, but decided that 
each analysis should be it's own repo, housed within each users account
and forked by the overall `organization` (maybe GotelliLab or Climate-Cascade)

Upon transfering development to this repo, dropped practice of individual
shell scripts for each step (QC, normalization, assembly) in favor of a
single R script because (A) I'm most proficient with R and (B) I really 
like the `cache` feature of `knitr` that allows me to fully implement a 
code chunk and then not re-run it, but continue analysis within the same 
file.

Repository contains scripts for assembly of the *Aphaenogaster* 
transcriptome.

**Data**
  
  - 2 lanes of Illumina HiSeq 2x100bp paired-end reads
  - 1 lane for each of 2 colonies
  - 12 treatments per lane
    * heat shock every 3.5C from 0 to 38C

**Analysis**

  - QC
  - Merge overlapping reads
  - Digital normalization 
  - Assembly with Trinity
  - Reduce transcript redundancy and chimeric transcripts
  - Map reads for each sample (treatment by colony) to transcriptome
  - Functional annotation





