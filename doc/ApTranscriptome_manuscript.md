Thermal reactionome of a common ant species
============================================

**Authors**: John Stanton-Geddes, Andrew Nguyen, Rob Dunn, Aaron Ellison, Nate Sanders, Nick Gotelli, Sara Helms-Cahan  

## Introduction

Temperature plays a prominent role in structuring biological diversity from molecular to macroecological scales [@Kingsolver2009]. At the population level, thermal tolerance has been assessed in dozens of studies [...] often demonstrating local adaptation [...]. In combination with environmental data, this thermal tolerance data has provided the backbone to dozens of studies examining limits to species geographic ranges [..], and the potential for species to track [...] and adapt to warming climates [...and climate change impacts on population growth [@Deutsch2008, Kingsolver2013]. At the organismal level, studies of thermal tolerance have largely focused on candidate genes such as heat shock proteins [@Dahlgaard1998, ...]. However, estimates of the standing genetic variation [@Krebs1997, Krebs1997, Williams2012], studies of differential gene expression [@Leemans2000], quantitative trait loci (QTL) [@Morgan2006] and knock-outs [@Takahasi2011] in the model system *Drosophila melanogaster* all implicate dozens of additional genes underlying thermal tolerance. Little is known about the diversity of these genes, especially in non-model species. 

Ants are ecologically-important species involved in multiple symbioses and ecosystem-functioning processes. For example, in eastern North America, the ant species *Aphaenogaster rudis* disperses seeds for 30% of herbacious species [@...]. In addition, ants are key contributors to nutrient cycling [...] and .... . As ectothermic species, ants experience a wide range of temperature conditions. Throughout the year, ant colonies persist through freezings and hot summers, and at shorter time scales, individual workers can experience substantial thermal shock as they leave the buffered ant nest to forage for food [...]. Thus, we expect that at the organismal level, ants utilize multiple molecular pathways to deal with thermal stress in addition to colony-level traits such as changes in foraging behavior [@Gordon2013] and nest location [...]. 
 
To explore the molecular toolkit used by ant species to cope with thermal stress, we performed a whole organism RNA sequencing (RNAseq) study using two colonies in the *Aphaenogaster fulva-picea-rudis* complex. One colony was from the interior of the species range in North Carolina, while the other was from a location approaching the northern range edge in Vermont. Our specific goals were to 1) generate the set of all genes (i.e. transcriptome) expressed by two colonies in the *Aphaenogaster fulva-picea-rudis* complex across physiologically-relevant temperatures, 2) identify genes that have responsive expression across this thermal gradient (i.e. thermal reactionome), 3) explore the functional annotations of these genes, and 4) examine differences in the thermal reactionome between the two colonies. Specifically, we hypothesize that the more southern of the two colonies will have higher constitutive expression of many of the genes that are up-regulated in response to heat shock of the more northern colony. Moreover, this study makes a methodological contribution as we analyze expression in a linear model, looking for *responsive* genes, rather than genes that are differentially-expressed.


## Methods

**Samples**

Colonies of *Aphaenogaster* were collected from Molly Bog, Vermont (University of Vermont Natural Areas; lat, lon) and Durham, North Carolina ( 36.037340° N -78.874140° W). The phylogeny of *Aphaenogaster* is currently unresolved, but putative identification based on morphology determined these colonies to be *A. picea* and *A. rudis*, respectively (Field Guide to Ants of New England). These colonies were maintained in the lab at 26C ?% humidity with a 12 hr light for 6 months prior to sample collection. 

From each species, we haphazardly collected 12 ants at the same time on 12 days. Each day, the 12 ants were placed in glass tubes into a water bath (...) at one of 12 randomly-assigned temperatures, every 3.5C between 0 and 38.5C, for one hour. The minimum and maximum temperatures were selected based on preliminary work showing that these temperatures are close to the critical minimum and maximum temperatures for an *Aphaenogaster* from VT (???) . At the end of the hour, the ants were flash frozen in liquid nitrogen and stored at -80C. mRNA was extracted from three pooled ants using a DNeasy kit (Qiagen, ...). Ants were homogenized in Buffer ATL with zirconium silicate beads in a Bullet Blender (...) and the standard Qiagen protocol was followed from this point. 

**Sequencing and transcriptome assembly**

For each species, the 12 samples were barcoded and sequenced in a single lane of 2 x 100bp paired-end reads on an Illumina HiSeq ???? at the Biomedical Genomics Research Center at the University of Minnesota, yielding 160 and 200 million reads for the Ar and Ap samples respectively. Reads were filtered to remove Illumina adapter sequences and low quality bases using the program Trim_galore! [cite], removing ... and ... reads for the Ar and Ap samples respectively.

We assembled the sequenced reads into the full set of mRNA transcript (the `transcriptome`) for each species using the Trinity *de novo* transcriptome assembly program [@Grabherr]. Prior to assembling with Trinity, we performed two steps to increase assembly quality and decrease computational resources. First, as the DNA was size selected to 200 bp for sequencing and we sequenced 2 x 100 bp reads, we expected some fraction of the reads to be overlapping. The program `FLASH` identifies and merges overlapping reads. For our data, xx% and xx% of reads were merged into extended fragments for Ar and Ap, respectively. Second, in mRNAseq studies, high coverage data is necessary to capture lowly expressed genes, but results in highly uneven coverage across the transcriptome. Only 7-20x coverage is adequate for accurate assembly. We applied digital normalization (implemented in the *khmer* python library [@diginorm]) to remove reads with greater than 20x coverage (-k 20). For our data, digital normalization retained ... (x%) and ... (x%) of the total reads. Moreover, Illumina reads are known to have about a 0.1% error rate, and to correct for this we removed rare reads using the `filt-abund.py` script in *khmer*. These filtering steps make our assembly conservative by removing reads containing sequencing errors. In addition, after this stringent filtering, we were able to assemble both transcriptomes in under x hours using less than 50GB of memory, an order of magnitude improvement over using all the data for assembly. 

*De novo* transcriptome assembly, partially due to the absence of a genome reference, is prone to falsely identifying alternative transcripts and identifying inaccurate transcripts that are chimeric (e.g. regions of two separate transcripts with enough overlap to assemble into a false, or chimeric, third transcript). To remove this redundancy, we use the program `CAP3` to cluster sequences with greater than 98% similarity and merge transcripts with overlaps greater than ?? and 98% similar in length. `CAP3` further identifies and merges chimeric transcripts by ... . Specific parameters can be found in the analysis script (link below), which reduced the number of transcripts to YYYY and YYYY for Ar and Ap, respectively. 

** Functional annotation and gene ontology **

To determine the function of the transcripts, we 
	
	* `tblastx` to the NCBI non-redundant protein database.
	* Trinotate
	* blast2GO

**Thermally-responsive genes**

We quantified expression of each transcript using the program `Sailfish` [@Patro]. `Sailfish` reports gene expression as transcripts per million (TPM) which is the preferred measure of gene expression as values are comparable among sequencing runs [@...]. In addition, `Sailfish` corrects for GC-content and ... as documented in ... [@...]. 

To identify genes that have significant changes in expression across the thermal gradient, we fit the linear model

$y = x + x^2$
    
where *y* is the expression level and *x* is the temperature, for every transcript. We retained only those transcripts with an overall significant model hit ... As we performed 80k and 90k unique tests for A22 and Ar, respectively, we applied false discovery rate (FDR) to correct for multiple testing. Specifically, we applied the method of Storey and Tibrashani [@...] where ... 

While we do not have replication at any single temperature, under the assumption that errors are independent with respect to temperature, the linear model approach is appropriate ... noise should swamp signal ... **NICK - COMMENT HERE???**

**Comparative gene expression**

We identified orthologous transcripts between the two assemblies using a `tblastx` reciprocal best hit. Specifically, we translated the fasta sequences into the six possible amino acid frames, performed a BLAST search against the alternate assembly, and retained the transcripts that were each other's best hit. This resulted in YYYYY orthologous transcripts. For these orthologs, we examined differences in expression patterns between the two colonies... 

**Data and analysis availability**

The raw Illumina reads are available at [...], the script for analysis and version history is available at []http://github.com/johnstantongeddes/ApTranscriptome], and version 0.1 of the transcriptome assemblies can be downloaded from [...].  


## Results

The transcriptome assembly for each species resulted in XX,XXX and YY,YYY transcripts (Table 1). 

Of these, we identified 888888 and 77777 transcripts that showed thermally-responsive expression, though only 11111 22222 of these were significant after FDR. 

The classes of genes over-represented included....


## Discussion




---------------------


Understanding the genetic basis of thermal tolerance will improve our ability to predict the effects of climate change on species abundances and distributions. Specifically, genetic diversity of key genes involved in thermal tolerance can be assayed across a species range [...] , giving a proxy for genetic variance in this trait. As the the potential for evolution to prevent populations from going extinct in response to environmental change (i.e evolutionary rescue) depends in part on the genetic variation for thermal tolerance within populations [@Gomuclkiewcz2013].   
 . This is especially true for ectothermic species which are predicted to be highly sensitive to temperature changes ...  


As ectothermic species, the fitness of ant colonies is highly susceptible to variation in temperature. Recent models based on thermal performance curves indicate that many temperate and tropical ectothermic species will be near their physiological limits by 2100 [@Kingsolver2013]. However, these models assume that natural selection will not act to change the thermal limits of this species. 

The study of whole organism gene expression (e.g. transcriptomics) has rapidly increased in the past decade due to improved sequencing technology (RNAseq). However, to date RNAseq studies have focused on differential gene expression between treatment groups (e.g. ANOVA) [...]. While this approach is appropriate for situations where genes are either active or not (e.g. response to infection), for continuous environmental conditions such as temperature gene expression is more likely to be continuous in nature. When only comparing conditions at the ends of a continuum (e.g. 'low' vs 'high' temperatures), genes involved in adaptation to conditions at optimal temperatures will be missed. In this study, we use an alternative approach. We measure gene expression across a thermal gradient and identify those genes that show a significant change in gene expression. This approach is powerful as identify not only genes that have large differential changes in expression between 'low' and 'high' conditions, but also genes with maximum intermediate expression.
