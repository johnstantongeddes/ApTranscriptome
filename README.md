*Aphaenogaster* transcriptome assembly
=======================================

*Forked* from [github.com/johnstantongeddes/climate-cascade/projects/ApTranscriptome]

Initial development of assembly in that repository, but decided that 
each analysis should be it's own repo, housed within each users account
and forked by the overall `organization` (maybe GotelliLab or Climate-Cascade)

Repository contains scripts for assembly of the *Aphaenogaster* 
transcriptome and identification of thermally-responsonsive genes.

**Data**
  
  - 2 lanes of Illumina HiSeq 2x100bp paired-end reads
  - 1 lane for each of 2 colonies
  - 12 treatments per lane
    * heat shock every 3.5C from 0 to 38C

**Analysis**

  - QC
  - Assembly with Trinity
  - Reduce transcript redundancy and chimeric transcripts
  - Map reads for each sample (treatment by colony) to transcriptome
  - Functional annotation
  - Identification of responsive genes





