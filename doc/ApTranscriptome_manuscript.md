Thermal reactionome of a temperate ant species
============================================

**Authors**: John Stanton-Geddes, Andrew Nguyen, Rob Dunn, Aaron Ellison, Nate Sanders, Nick Gotelli, Sara Helms-Cahan  

## Introduction

Temperature plays a prominent role in structuring biological diversity from molecular to macroecological scales [@kingsolver2009]. At the population level, thermal tolerance has been assessed in dozens of studies [...] often demonstrating local adaptation [...]. In combination with environmental data, this thermal tolerance data has provided the backbone to dozens of studies examining limits to species geographic ranges [..], and the impacts of climate change on population viability [@deutsch2008; @kingsolver2013]. At the organismal level, studies of thermal tolerance have largely focused on candidate genes such as heat shock proteins [@dahlgaard1998]. However, estimates of the standing genetic variation [@krebs1997; @krebs1997; @williams2012], studies of differential gene expression [@leemans2000, @sorensen2005, @teets2012], quantitative trait loci (QTL) [@morgan2006] and knock-outs [@takahashi2011] in the model system *Drosophila melanogaster* all implicate dozens of additional genes underlying thermal tolerance. Little is known about the diversity of these genes, especially in non-model species. 

Ants are ecologically-important species involved in multiple symbioses and ecosystem-functioning processes. In eastern North America, ants in the *Aphaenogaster rudis-picea* species complex are some of the most abundant found in hardwood forests [@lubertazzi2012] where they disperse seeds for 30% of herbaceous species [@beattie1985]. In addition, ants are key contributors to litter processing and decomposition [@lyford1963; @folgarait1998, @terbough2001, @holec2005]. As ectothermic species, ants experience a wide range of temperature conditions. Throughout the year, ant colonies persist through freezings and hot summers, and at shorter time scales, individual workers can experience substantial thermal shock as they leave the buffered ant nest to forage for food [...]. Thus, at the organismal level we expect that ants utilize multiple molecular pathways to deal with thermal stress in addition to species-level traits such as changes in foraging behavior [@gordon2013] and nest location [...]. 
 
To explore the molecular toolkit used by *Aphaenogaster* to cope with thermal stress, we performed a whole organism RNA sequencing study using the sister species *A. picea* and *A. rudis*. Our specific goals were to 1) generate the set of all transcripts (i.e. transcriptome) expressed by *A. picea* and *A. rudis* across physiologically-relevant temperatures, 2) identify transcripts that have responsive expression across this thermal gradient (i.e. thermal reactionome), 3) explore the functional annotations of these transcripts, and 4) examine differences in the thermal reactionome between the two species. Specifically, we hypothesized that the more southern of the two species, *A. rudis*, would have higher constitutive expression of many of the transcripts that are up-regulated in response to heat shock of the more northern species. We performed our analysis at the transcript level due to limitations of working with a system lacking a reference genome, but refer to these transcripts as genes for simplicity. Finally, this study makes a methodological contribution as we identify responsive genes using a regression approach, rather than simply identifying genes that are differentially-expressed.


## Methods

**Samples**

Colonies of the *Aphaenogaster picea* and *rudis* were collected from Molly Bog, Vermont (University of Vermont Natural Areas; lat, lon) and Durham, North Carolina ( 36.037° N -78.874° W), designated as 'ApVT' and 'ArNC' respectively. The phylogeny of *Aphaenogaster* is currently unresolved, but putative identification based on morphology determined these colonies to be *A. picea* and *A. rudis*, respectively (Bernice Bacon DeMarco, pers. communication). These colonies were maintained in the lab for 6 months prior to experimentation. 

From each species, we haphazardly collected 12 ants at the same time on 12 days. Each day, the 12 ants were placed in glass tubes into a water bath (...) at one of 12 randomly-assigned temperatures, every 3.5C between 0 and 38.5C, for one hour. The minimum and maximum temperatures were selected based on preliminary work showing that these temperatures are close to the critical minimum and maximum temperatures for an *Aphaenogaster* from VT (???) . At the end of the hour, the ants were flash frozen in liquid nitrogen and stored at -80C. mRNA was extracted from three pooled ants using a DNeasy kit (Qiagen, ...). Ants were homogenized in Buffer ATL with zirconium silicate beads in a Bullet Blender (...) and the standard Qiagen protocol was followed from this point. 

**Sequencing, assembly and annotation**

For each species, the 12 samples were barcoded and sequenced in a single lane of 2 x 100bp paired-end reads on an Illumina HiSeq 1500 yielding 200 and 160 million reads for the ApVT and ArNC samples respectively. Reads were filtered to remove Illumina adapter sequences and low quality bases using the program Trimmomatic [@lohse2012]. 

We assembled the sequenced reads into the full set of mRNA transcripts (the `transcriptome`) for the combined data set from both colonies using the Trinity *de novo* transcriptome assembly program [@grabherr2011]. *De novo* transcriptome assembly is prone to falsely identifying alternative transcripts and identifying inaccurate transcripts that are chimeric (e.g. regions of two separate transcripts with enough overlap to assemble into a false, or chimeric, third transcript) [@yang2013] To account for this, we first ran program `CAP3` [@huang1999] to cluster sequences with greater than 90% similarity and merge transcripts with overlaps longer than 100 bp and 98% similar in length, and second ran the program `uclust` which cluster sequences completely contained within longer sequences at greater than 90% similarity (see Supplemental for exact parameters). We used liberal values (90% similarity) to merge transcripts to account for the fact that we were including sequencing data from two colonies of the same species, and wanted to avoid assembling orthologs individually. 

We performed functional annotation of the transcriptome assembly using the web-based tool `FastAnnotator` [@chen2012] which annotates and classifies transcripts by Gene Ontology (GO) term assignment, enzyme identification and domain identification. 


**Thermally-responsive transcripts**

We quantified expression of each transcript using the program `Sailfish` [@patro2013]. `Sailfish` reports transcript expression as transcripts per million (TPM) which is the preferred measure of transcript expression as values are comparable among sequencing runs [@wagner2012]. In addition, `Sailfish` corrects for sequence composition bias and transcript length. As preliminary examination of the data (see supplemental...) indicated that 7deg C samples may have been mis-labeled, we omitted these data from the analysis. The expression values were highly correlated between colonies at each temperature treatment (r^2 > 0.98) indicating that assembling the transcriptome with data from both colonies was justified (Supplemental).

To identify transcripts that have significant changes in expression across the thermal gradient, we fit the linear regression

$$ TPM = \beta_0 + \beta_1(species) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(species * temp) + \beta_5(species * temp^2) + \epsilon $$
    
independently to each transcript. For a continous predictor such as temperature, this regression approach is preferred to an ANOVA approach as it can reveal non-linear responses such as hump-shaped or threshold effects even  [@cottingham2005]. Moreover, as we expect errors in the read count distribution to be independent with respect to temperature, our method is robust to issues of overdispersion in read count data [@anders2010]. To correct for multiple testing, we applied the False Discovery Rate (FDR) approach of Benjamini and Hochberg [-@benjamini1995] at 5% FDR.

Of these significant transcripts, we focused on two subsets. First, to describe the overall molecular toolkit for thermal tolerance in *Aphaenogaster*, we identified the subset of transcripts under our FDR threshold that had significant $\beta_2(temp)$, $\beta_3(temp^2)$, $\beta_4(species * temp)$ or $\beta_5(species * temp^2)$ terms in the linear regression. For this set of transcripts, we specifically looked at the proportion of transcripts that were expressed at both low (< 10C) and high (> 31C) temperatures. To select this group, we required that (a) the transcript had a convex shape and (b) expression at both low and high values was two standard deviations greater than the minimum expression.

Second, to explore differences in thermal responsiveness between the two colonies, we selected the subset of transcripts under our FDR threshold that had significant $\beta_4(species * temp)$ or $\beta_5(species * temp^2)$ terms in the linear regression. We compared levels of constitutive expression for transcripts that showed responsiveness in one species, but not the other. 


**Gene set enrichment analysis**

We explored enrichment of gene set classes in the (1) thermally-responsive and (2) species specific thermally-responsive transcripts using the program `topGO` [@alexa2006; @alexa2010] with the `parentChild` algorithm [@grossman2007]. Briefly, this approach identifies GO terms that are overrepresented in the significant transcripts relative to all GO terms in the transcriptome, accounting for dependencies among the GO terms. 


**Data and analysis availability**

The raw Illumina reads are available at [...], the script for analysis and version history is available at (http://github.com/johnstantongeddes/ApTranscriptome).  


## Results

**Assembly**

The Trinity *de novo* transcriptome assembly included 126,172 transcripts with a total length of 100 million bp. Filtering to remove redundant or chimeric reads resulted in an assembly with 96,253 contigs and a total of 105,536 transcripts. For all transcripts, the total length was 63 million bp with an N_[50] length of 895 bp and a mean contig size of 593 bp. 


**Functional annotation**

Of the 105,536 filtered transcripts, 55,432 had hits to the NCBI-nr database. Of these, 38,711 transcripts mapped to GO terms, 1,659 transcripts were identified to an enzyme and 18,935 transcripts mapped to a domain with >50% coverage. Of these, 5,787 transcripts are annotated to putative genes in the genome sequence of the ant species *Solenopsis invicta* and XXXX are annotated to genes previously identified from organisms within the class Insecta. Over-represented functional words in the transcriptome annotation are shown in figure 1.

![Word cloud of functional annotation terms in the *Aphaenogaster* transcriptome](../wordcloud2.png)


**Thermal responsive transcripts**

We identified 22,582 transcripts with overall model fits of P < 0.05, retaining 8,753 after correcting for FDR, with 5,540 significantly changing expression with temperature. Of these, 822 increase with temperature, 2,058 decrease with temperature, 710 have maximum expression at intermediate temperatures and 1,990 are expressed at both low and high temperatures (Fig. 1).

![Smoothed line plots of scaled expression values for the thermally-responsive transcripts that show bimodal (n = 1,990), intermediate (n = 710), high (n = 822) and low temperature (n = 2,058) expression for the northern (ApVT) and southern (ArNC) *Aphaenogaster* colonies.](../figure/plot_responsive4.png)

Two functional groups of *a priori* interest showed patterns consistent with our expectations. Responsive transcripts annotated to "heat shock proteins" were up-regulated at high temperatures (>31 C) in the northern *A. picea* colony, but not up-regulated in the putatively more heat tolerant *A. rudis* colony (Fig. 3a). Conversely, genes involved in the biological process of dormancy were up-regulated at cold temperatures (< 10 C) in the southern *A. rudis* colony, but not *A. picea* (Fig. 3b). As these genes represent only a handful of the responsive transcripts, the expression patterns of all responsive transcripts can be searched and viewed interactively at [https://johnsg.shinyapps.io/ApRxN-shinyapp/].

![Figure 3. goes here]

**Gene set enrichment analysis**

We identified a total of 106 enriched GO terms. Of these, 12 were enriched at high temperatures, 26 at low temperatures, 52 at both high and low temperatures and 16 at intermediate temperatures. At high temperatures, the enriched GO terms included terms implicated in stress tolerance such as 'regulation of response to stimulus' and 'positive regulation of ubiqutination, while at low temperatures enriched GO terms included 'detection of stimulus' and 'dormancy response' (Table 1, Fig. 3). Terms enriched at both high and low temperatures included 'regulation of response to stimulus', 'positive regulation of protein ubiquitination', and 'negative regulation of programmed cell death'. At intermediate temperatures, enriched GO terms were primarily related to glycosylation and development (e.g. mitosis). 

![Figure 4. Comparison cloud of the GO terms enriched at high, low, intermediate and both high and low (bimodal) temperatures.](../gsea_comparison_cloud.png)

**Differences in expression between colonies**

Of the responsive transcripts, 3,760 had a significant interaction between temperature and species. Overall, expression was higher in the ... species than .... For transcripts that showed increased expression at high temperatures in the *ApVT* species, the main effect coefficient was ... indicating that ... Conversely, for transcripts that showed increased expression at low temperatures in the *ArNC* species, the main effect of the *ApVT* species was ... indicating ... 

Individual genes differed in expression among colonies. For example, both heat shock proteins increased expression at high temperatures in the *ApVT* species, but were constant in the *ArNC* species (Figure ...). 


## Discussion

Thermal tolerance is a key trait for determining the survival of individual organisms, the persistence of populations and the geographic limits of species. While considerable research has focused on a few key genes (e.g. heat shock proteins) involved in thermal stress response [...], little is known about the functional diversity of thermal responsiveness. In this study, we characterized the thermal reactionome of sister ant species, *Aphaenogaster picea* and *Aphaenogaster rudis*. We found that nearly 10% of the transcripts in this species show significant thermal responsiveness. Moreover, we found ... differences in patterns of expression among these two species. 

Responsiveness ... a temperate ant species are thermally responsive. The majority of these transcripts (about 90%) are up-regulated at high or low temperatures. Discuss the 26 "stress response" genes

Previous studies have using *Drosophila* [@sorensen2005] and a temperate flesh fly [@teets2012b] have implicated GO terms .... as responding to heat stress. In this study, we found .... of these GO terms to be enriched at high temperatures for *Aphaenogaster*. However, when we looked at the exact identify of the genes, we only found X overlapping genes.

This result poses a paradox. Why is the heat shock functional response highly conserved at the biological process level, while the exact molecular pathways appear to be highly labile? We suggest two possible explanations. First, regulatory networks controlling gene expression are evolutionarily dynamic. Trade-offs between function in different organisms may constrain which genes are involved in thermal response, thus resulting in different pathways to the same phenotype. Second, variation in molecular pathways may be due to inconsistent selection on thermal stress response. That is, as the environment changes and species shift their distributions, selection may relax or strengthen repeatedly within and across lineages. Different genes may be recruited within each selective episode. Differentiating between these, and alternate, explanations will require greater geographic and population genetic sampling.

A key outcome of this study was developing the first genomic tools for the study of adaptive differentiation in *Aphaenogaster*. Our comparative results are limited to two populations, and an important future step will be expanding this work to eaexamine differences in gene expression among a larger number of populations. Moreover, the potential for species to adapt to climate change will depend in part on their genetic variation for thermal tolerance. The quantitative genetic studies required to estimate this adaptive potential have been out of reach of ant biologists due to the long lifespan of colonies, and inability of ants to mate in the lab. However, future work could leverage genomic information to infer the relatedness among colonies in the absence of a known pedigree, enabling estimation of genetic variation within populations and their adaptive potential. 

This study also provides a methodological contribution to the field. The study of whole organism gene expression (e.g. transcriptomics) has rapidly increased in the past decade due to improved sequencing technology (RNAseq). However, to date RNAseq studies have focused on differential gene expression between treatment groups (e.g. ANOVA) [...]. While this approach is appropriate for situations where genes are either active or not (e.g. response to infection), for continuous environmental conditions such as temperature gene expression is more likely to be continuous in nature. When only comparing conditions at the ends of a continuum (e.g. 'low' vs 'high' temperatures), genes involved in adaptation to conditions at optimal temperatures will be missed. In this study, we use an alternative approach. We measure gene expression across a thermal gradient and identify those genes that show a significant change in gene expression. This approach is powerful as identify not only genes that have large differential changes in expression between 'low' and 'high' conditions, but also genes with maximum intermediate expression.

In summary, we have generated the first transcriptome of a temperate ant species, and have identified approximately 10% of transcripts that show thermal responsiveness. In this group, we find many known stress response genes, but also identify novel genetic pathways that respond to temperature. These results provide key insights on the genetic basis of thermal tolerance, and provide a tremendous resource for the future study of ecological adaptation in ant species.


## References

