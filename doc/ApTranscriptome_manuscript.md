Thermal reactionome of a temperate ant species
============================================

**Authors**: John Stanton-Geddes, Andrew Nguyen, Rob Dunn, Aaron Ellison, Nate Sanders, Nick Gotelli, Sara Helms-Cahan  

## Introduction

Temperature plays a prominent role in structuring biological diversity from molecular to macroecological scales [@kingsolver2009]. At the population level, thermal tolerance has been assessed in dozens of studies often demonstrating adaptive differentiation [@hoffmann2013]. In combination with environmental data, this thermal tolerance data has provided the backbone to dozens of studies examining limits to species geographic ranges [@angert2011; @diamond2013; @gouveia2014], and the impacts of climate change on population viability [@deutsch2008; @kingsolver2013]. At the organismal level, studies of thermal tolerance have largely focused on candidate genes such as heat shock proteins [@dahlgaard1998]. However, estimates of the standing genetic variation [@krebs1997; @krebs1998; @williams2012], studies of differential gene expression [@leemans2000; @sorensen2005; @teets2012], quantitative trait loci (QTL) [@morgan2006] and knock-outs [@takahashi2011] in the model system *Drosophila melanogaster* all implicate dozens of additional genes underlying thermal tolerance. Little is known about the diversity of these genes or their expression patterns, especially in non-model species.

The study of whole organism gene expression (e.g. transcriptomics) has rapidly increased due to improved sequencing technology (RNAseq). To date, RNAseq studies have largely focused on differential expression across few environments [@teets2012; @pespeni2013; @barshis2013; @fu2014] or tissues [@li2013] using an ANOVA-type statistical approach. However, many ecologically-relevant aspects of the environment such as temperature vary continuously in nature. Variation in the phenotype of an organism in response to such a continuous predictor (i.e. phenotypic plasticity) is referred to as a reaction norm. A reaction norm can differ not only in mean values, but also in the shape of the curve as is often found for thermal reaction norms [@hoffmann2013; @murren2014]. An analysis of phenotypic reaction norms found that differences in the shape were greater than differences in mean values in comparisons among species and populations [@murren2014]. The authors concluded that many patterns of biologically-relevant plasticity are likely to be missed in studies with few environments. Similarly, gene expression studies with few environments are likely to miss many biologically-relevant changes in gene expression. To overcome this limitation, we use a novel method to characterize the reaction norm of all genes in the transcriptome (e.g. the reactionome) in response to thermal stress. Specifically, we measure gene expression across a thermal gradient and identify thermally-responsive genes using a regression approach.

We apply this approach in a study of two closely related eastern North American ant species, *Aphaenogaster rudis* and *A. picea* [@umphrey1996], which have more southern and northern distributions, respectively. Ants are ecologically-important species involved in multiple symbioses¸ and are key contributors to seed dispersal [@lengyel2009; @rodriguezcabal2012], litter processing and decomposition [@lyford1963; @folgarait1998; @terborgh2001; @holec2006]. As ectothermic species, ants experience a wide range of temperature conditions. Throughout the year, ant colonies persist through freezings and hot summers, and at shorter time scales, individual workers can experience substantial thermal shock as they leave the buffered ant nest to forage for food. Thus, at the organismal level we expect that ants utilize multiple molecular pathways to tolerate thermal stress.

To explore the molecular toolkit used by *Aphaenogaster* to cope with thermal stress, we performed a whole organism RNA sequencing study. Our specific goals were to 1) generate the set of all genes (i.e. transcriptome) expressed by *Aphaenogaster* across physiologically-relevant temperatures, 2) identify genes that have responsive expression across this thermal gradient (i.e. thermal reactionome), 3) explore the functional annotations of genes that are activated or down-regulated in response to cold and high temperatures. We expected to observe candidate genes in the *response to stress* Gene Ontology (GO) category including *heat shock proteins* and ... in the thermally-responsive group.   

While characterizing the reactionome is useful for description of the molecular pathways for thermal tolerance, a greater understanding of the process of evolution can be gained using a comparative approach [@harvey1991]. Specifically, comparisons of populations or species in contrasting environments can reveal the extent to which geographic differentiation and selection has changed the number of genes involved in thermal tolerance. To accomplish this, we compared the thermal reactionomes between the two species of *Aphaenogaster*. *A. picea*, the northern species, experiences less chronic and acute heat stress than *A. rudis*. Thus, we that *A. picea* would have more genes than *A. rudis* that are up-regulated to high temperatures, while *A. rudis* would have higher constitutive expression of these thermally-tolerant genes. Reciprocally, we expected *A. rudis* to have greater up-regulation of genes in response to cold temperatures, while *A. picea* would have greater constitutive expression of cold-tolerant genes. Finally, phenotypic studies in *Aphaenogaster* [@warren2013] and other ecotherms [@hoffmann2013] have found less variation in critical maximum than minimum temperatures. From this, we hypothesized that if the species had more similar maximum temperatures than minimum temperatures, the northern species *A. picea* would have a broader thermal breadth of gene expression than *A. rudis*. 


## Methods

### Samples

Colonies of the *Aphaenogaster picea-rudis* complex were collected from Molly Bog, Vermont (University of Vermont Natural Areas; lat, lon) and Durham, North Carolina ( 36.037° N -78.874° W), designated as *ApVT* and *ApNC* respectively. The phylogeny of *Aphaeogaster* is currently unresolved [@umphrey1996]. Though preliminary work indicates that the northern part of the distribution (including *ApVT*) to be a distinct clade and possibly species (Bernice Bacon DeMarco, pers. communication), we will refer to these as distinct colonies rather than species. These colonies were maintained in the lab for 5 months prior to experimentation. 

From each species, we haphazardly collected 12 ants at the same time on 12 days. Each day, the 12 ants were placed in glass tubes into a water bath (...) at one of 12 randomly-assigned temperatures, every 3.5°C between 0° and 38.5°C, for one hour. The minimum and maximum temperatures were selected based on preliminary work showing that these temperatures are close to the critical minimum and maximum temperatures for an *Aphaenogaster* from VT (A. Nguyen, unpublished data) . At the end of the hour, the ants were flash frozen in liquid nitrogen and stored at -80°C. mRNA was extracted from three pooled ants using a DNeasy kit (Qiagen Inc; Valencia, CA). Ants were homogenized in Buffer ATL with zirconium silicate beads in a Bullet Blender (Next Advance; Averill Park, NY) and the standard Qiagen protocol was followed from this point. 

### Sequencing, assembly and annotation

For each species, the 12 samples were barcoded and sequenced in a single lane of 2 x 100bp paired-end reads on an Illumina HiSeq 1500 yielding 200 and 160 million reads for the ApVT and ArNC samples respectively. Reads were filtered to remove Illumina adapter sequences and low quality bases using the program Trimmomatic [@lohse2012]. 

We assembled the sequenced reads into the full set of mRNA transcripts (the `transcriptome`) for the combined data set from both colonies using the Trinity *de novo* transcriptome assembly program [@grabherr2011]. *De novo* transcriptome assembly is prone to falsely identifying alternative transcripts and identifying inaccurate transcripts that are chimeric (e.g. regions of two separate transcripts with enough overlap to assemble into a false, or chimeric, third transcript) [@yang2013] To account for this, we first ran program `CAP3` [@huang1999] to cluster sequences with greater than 90% similarity and merge transcripts with overlaps longer than 100 bp and 98% similar in length, and second ran the program `uclust` which cluster sequences completely contained within longer sequences at greater than 90% similarity (see Supplemental for exact parameters). We used liberal values (90% similarity) to merge transcripts to account for the fact that we were including sequencing data from two colonies of the same species, and wanted to avoid assembling orthologs individually. We performed subsequent analyses at the transcript level due to limitations of working with a system lacking a reference genome, but refer to these transcripts as genes for simplicity.

We performed functional annotation of the transcriptome assembly using the web-based tool `FastAnnotator` [@chen2012] which annotates and classifies genes by Gene Ontology (GO) term assignment, enzyme identification and domain identification. 


### Thermally-responsive genes

We quantified expression of each gene using the program `Sailfish` [@patro2013]. `Sailfish` reports gene expression as transcripts per million (TPM) which is the preferred measure of gene expression as values are comparable among sequencing runs [@wagner2012]. In addition, `Sailfish` corrects for sequence composition bias and gene length. As preliminary examination of the data (see supplemental...) indicated that 7°C samples may have been mis-labeled, we omitted these data from the analysis. The expression values were highly correlated between colonies at each temperature treatment (r^2^ > 0.98) indicating that assembling the transcriptome with data from both colonies was justified (Supplemental).

To identify genes that have significant changes in expression across the thermal gradient, we fit the regression model

$$ log(TPM + 1) = \beta_0 + \beta_1(species) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(species * temp) + \beta_5(species * temp^2) + \epsilon $$
    
to each gene. We used $log(TPM + 1)$ as the response to control for skew in the expression data. For a continuous predictor such as temperature, this regression approach is preferred to an ANOVA approach as it can reveal non-linear responses such as hump-shaped or threshold effects even [@cottingham2005]. This method is robust to over-dispersion as we expect errors in the read count distribution [@anders2010] to be independent with respect to temperature. To correct for multiple testing, we calculated adjusted *P* values using the False Discovery Rate (FDR) approach of Benjamini and Hochberg [-@benjamini1995], retaining genes under the 5% FDR threshold.

We identified thermally-responsive genes as the subset under our FDR threshold that had significant $\beta_2(temp)$, $\beta_3(temp^2)$, $\beta_4(species * temp)$ or $\beta_5(species * temp^2)$ terms after step-wise model selection by AIC. All analyses were performed in R [@r_core_team2013] and are fully reproducible (Supplemental file). 

For each thermally-responsive gene, we predicted expression levels across the tested thermal range using the final linear model separately for each colony. We used the predicted gene expression values to partition the thermally-responsive genes into five expression categories for each colony; genes that had greatest expression at high (> 31°) temperatures (*High*), low (< 10°C) temperatures (*Low*), intermediate (10 - 30°C) temperatures (*Intermediate*), both high and low temperatures (*Bimodal*) or were not thermally-responsive in that colony (*NotResp*). For the *Bimodal* group, we required that expression at both low and high temperatures was at least one standard deviation greater than the expression at the mean temperature of 19.5°C.

To describe the overall molecular toolkit for thermal tolerance in *Aphaenogaster*, we performed gene set enrichment analysis (GSEA) on the intersection of each expression category in the two colonies using the program `topGO` [@alexa2006; @alexa2010] with the `parentChild` algorithm [@grossmann2007] implemented in R [@r_core_team2013]. Briefly, this approach identifies GO terms that are overrepresented in the significant genes relative to all GO terms in the transcriptome, accounting for dependencies among the GO terms. 


### Colony-level comparisons

To gain greater insight on the genetic basis of thermal tolerance in *Aphaenogaster*, we performed comparative analyses between the two colonies focusing on the thermally-responsive genes that had a significant colony x temperature interaction. For each gene, we predicted expression levels across temperatures for each colony using the full linear model, and then grouped genes into the four responsive categories, *High*, *Low*, *Intermediate*, and *Bimodal*, as well as a fifth category *Not Expressed*. To focus on differences in gene expression between the two colonies, we performed GSEA for the set difference of each category.

With the same set of genes, we explored the extent to which differences in the thermal reactionome between the colonies were due to changes in the mean or shape of the reaction norm for each individual gene. Following the example of Murren et al [-@murren2014], we tested for overall differences in the (a) mean expression value across temperatures (*M*), (b) slope of expression (*S*), (c) curvature of expression (*C*) and (d) wiggle of expression (*W*), which consists of all higher-order differences in shape not capture by the first three measures (formulas given in Appendix). As genes differed by five orders of magnitude in expression levels, we calculated the mean-standardized difference between the colonies $\Delta M = (Mean_{ApVT} - Mean_{Ar}) / \bar M$, where $\bar M = (Mean_{ApVT} + Mean_{Ar}) / 2$ yielding $\Delta M$, $\Delta S$, $\Delta C$, and $\Delta W$ for the mean, slope, curvature and wiggle of each gene, respectively. Using the mean-standardized values, we performed one-sample two-sided *t*-tests to determine if the overall mean, slope or curvature were greater in *ApVT* than *ApNC*. 

Further, to quantify the relative contribution of changes in the mean, slope, curvature or wiggle on the overall differences in reaction norms, we defined the *Total* difference in the reaction norm for each gene as $\Delta T = \Delta M + \Delta S + \Delta C + \Delta W$. We then partitioned the contribution of each measure by dividing it by $\Delta T$. 

For the subset of genes that responded to temperature in one colony (i.e. *High* expression in the *ApVT* colony)  but not the other, we hypothesized that geographically-divergent selection may have favored constitutive expression. To test this idea, we compared levels of constitutive expression between *A22* and *ApNC* near the optimum temperature of 19.5° C . Specifically, we predicted higher constitutive expression in *ApNC* for genes with *High* expression in *A22*, and conversely, higher constitutive expression  in *A22* for *Low* genes in *ApNC* . 

The *Intermediate* expressed genes are core molecular processes that are expressed at non-stressful temperatures, and shut-off when the organism experiences thermal stress. We hypothesized that if the more southern *ApNC* colony was more thermally-tolerant than *ApVT* colony, genes with 'Intermediate' expression (10-30° C) would be expressed across a wider range of temperatures. To test this with our data, for each *Intermediate* expressed gene in each colony, we randomly sampled 1,000 temperature values weighted by the expression function (e.g. point of maximum expression had highest probability of being sampled). We then calculated the standard deviation of this random draw, and performed a *t*-test comparing the standard deviations of expression between colonies.

In contrast, *Bimodal* expressed genes are those that are activated in response to thermal stress. Following the reasoning for *Intermediate* genes, we hypothesized that a more thermally-sensitive colony would activate these thermally-responsive genes earlier than a more thermally-tolerant colony. We calculated the standard deviation of expression as for *Intermediate* genes. The interpretation is different, as low standard deviation is a gene with nearly flat expression, while the maximum standard deviation of expression is for a highly-convex expression pattern. Thus, finding one colony had overall greater standard deviation of expression for *Bimodal* genes would be indicative of greater thermal sensitivity. 

**Data and analysis availability**

The raw Illumina reads are available at [...], the script for analysis and version history is available at (http://github.com/johnstantongeddes/ApTranscriptome).  



## Results




The Trinity *de novo* transcriptome assembly included 126,172 genes with a total length of 100 million bp. Filtering to remove redundant or chimeric reads resulted in an assembly with 96,253 contigs and a total of 105,536 genes. For all genes, the total length was 63 million bp with an N~50~ length of 895 bp and a mean contig size of 593 bp. Of the 105,536 filtered genes, 55,432 had hits to the NCBI-nr database. Of these, 38,711 genes mapped to GO terms, 1,659 genes were identified to an enzyme and 18,935 genes mapped to a domain with >50% coverage. Of these, 5,787 genes have best hits to putative genes in the genome sequence of the ant species *Solenopsis invicta* and XXXX are annotated to genes previously identified from organisms within the class Insecta.


|    &nbsp;     |  High  |  Low  |  Bimodal  |  Intermediate  |  NotResp  |
|:-------------:|:------:|:-----:|:---------:|:--------------:|:---------:|
|   **ApVT**    |  1179  | 4720  |   1534    |      660       |    399    |
|   **ApNC**    |  1163  | 3457  |    888    |      2376      |    608    |
|  **Overlap**  |  309   | 2417  |    294    |      352       |     0     |

Table: Number of transcripts with maximum expression at high, low, intermediate, both high and low (bimodal) temperatures or are not thermally-responsivefor each colony and their overlap.

% latex table generated in R 3.0.2 by xtable 1.7-3 package
% Fri Apr 11 10:23:28 2014
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrr}
  \hline
 & High & Low & Bimodal & Intermediate & NotResp \\ 
  \hline
ApVT & 1179 & 4720 & 1534 & 660 & 399 \\ 
  ApNC & 1163 & 3457 & 888 & 2376 & 608 \\ 
  Overlap & 309 & 2417 & 294 & 352 &   0 \\ 
   \hline
\end{tabular}
\end{table}


### Thermally-responsive genes

We identified 22,495 genes with overall model fits of P < 0.05, retaining 9,809 after correcting for FDR. Of these transcripts, 8,044 significantly change expression with temperature including 6,727 that have a colony x temperature interaction. The number of transcripts that increased expression with temperature (*High*), decreased expression with temperature (*Low*), increased expression at both high and low (*Bimodal*) temperatures, decreased expression at high and low temperatures (*Intermediate*) or were not expressed differed significantly between the colonies ($\chi^2 = 170, P < 0.0001$; Table 1). For both colonies, the majority of responsive transcripts were in the *Low* category (Fig. 1). More transcripts had *Intermediate* expression in *ApNC* than in *ApVT*, which had nearly twice as many bimodally-expressed transcripts  (Fig. 1, Table 1).  

![Probability density function of peak expression for all thermally-responsive in *ApNC* (blue) and *ApVT* (red) colonies](../results/PDF_expression_all.png)

As the molecular pathways involved in thermal response to cold and hot temperatures may differ, we performed gene set enrichment analysis on each expression category separately. Out of a total of 2,859 GO terms scored, we found 6 enriched (*P* < 0.01) for *High* transcripts, included categories implicated in stress response such as "gene expression" and "tissue regeneration" (Table S1). For *Low* genes, 60 GO terms were enriched including some implicated in stress response (e.g. detection of external stimulus; Table S1). For *Bimodal* genes, 11 GO terms were enriched including 10 "regulation" terms (e.g. regulation of primary metabolic process) indicating that these are general genes involved in activating responses to change (Table S1). Finally, for *Intermediate* genes 6 GO terms were enriched, all some type of metabolic or biosynthetic process (Table S1) consistent with the expectation that these are non-essential genes expressed in optimal conditions that are shut-down under thermal stress.

In addition, we examined the reaction norms for genes that we *a priori* expected to have significant thermal responses. In particular, we looked at the expression patterns of heat shock proteins including ...


### Colony-level comparisons

We quantified differences in the reaction norms of thermally-responsive genes between colonies as changes to due shifts in mean ($\Delta M$), slope ($\Delta S$), curvature ($\Delta W$) and wiggle ($\Delta W$). Overall, we found that there was no significant difference in position of the mean (*t* = -0.95, *P* = 0.34), but changes in the shape due to slope (*t* = 24.67, *P* < 0.001), curvature (*t* = 15.71, *P* < 0.001) and wiggle (*t* = -12.45, *P* < 0.001) were all highly significant. On average, slope (95% CI 0.41 - 0.48) and curvature (95% CI 0.25 - 0.33) were greater in *ApVT* than *ApNC*, indicating a steeper response to temperature. 

Of the total differences ($\Delta T$) in reaction norms between colonies for each gene, about 8% (95% CI 1 - 32%) was due to shifts in mean, 47% (95% CI 14 - 92%) due to shifts in slope, 39% (95% CI 0 - 67%) due to shifts in curvature and 0% (95% CI 0 - 87%) due to shifts in wiggle. 

Plasticity is a well-characterized mechanism by which organism tolerate environmental change. Given our sampling of a warmer (*ApNC*) and colder (*ApVT*) climate species, we hypothesized that genes which are plastically-expressed under heat shock in the colder species would have greater constitutive expression in the warmer species, and vice versa. Consistent with one of our predictions, we found that constitutive expression of *High* genes from *ApVT* was greater in the more southern *ApNC* colony (Wilcoxon V = 911970, *P* < 0.001, 95% CI location shift -0.11 - 0.03). However, for genes with *Low* expression in *ApNC* we also found that constitutive expression was greater in *ApNC* than *ApVT* (Wilcoxon V = 2041149, *P* < 0.001, 95% CI location shift -0.23 - 0.15). 

We tested for a transcriptome-wide pattern of thermal tolerance and sensitivity by examining the standard deviation of expression for the *Intermediate* expressed genes. Consistent with our expectations, the more southern *ApNC* colony had greater thermal breadth for *Intermediate* genes than the *ApVT* colony (*t*~df=1406~ = -9.25, *P* < 0.001, 95% CI difference in means 0.41 - 0.64). In contrast, there was no difference in thermal sensitivity as measured by the standard deviation of expression for *Bimodal* expressed genes (*t*~df=1406~ = 0.31, *P* = 0.76). 

We performed gene set enrichment analysis on the set difference of each expression category to reveal differences among the colonies response to thermal shock. For *High* genes, there were 4 enriched GO terms in *ApVT* and 12 in *ApNC*. For *Low* genes, there were 40 enriched GO terms in *ApVT* and 23 in *ApNC*. For *Intermediate* genes, there were 3 enriched GO terms in *ApVT* and 54 in *ApNC*. For *Bimodal* genes, there were 12 enriched GO terms in *ApVT* and 10 in *ApNC*. 


## Discussion

Thermal tolerance is a key trait for determining the survival of individual organisms, the persistence of populations and the geographic limits of species. While considerable research has focused on a few key genes (e.g. heat shock proteins) involved in thermal stress response [...], little is known about the functional diversity or transcriptome-wide patterns of thermal responsiveness. In this study, we characterized the thermal reactionome of sister ant species, *Aphaenogaster picea* and *Aphaenogaster rudis*. We found that nearly 10% of the genes in this species show significant thermal responsiveness. Moreover, we found ... differences in patterns of expression among these two species. 

Two functional groups of *a priori* interest showed patterns consistent with our expectations. Responsive genes annotated to "heat shock proteins" were up-regulated at high temperatures (>31 C) in the northern *A. picea* colony, but not up-regulated in the putatively more heat tolerant *A. rudis* colony (Fig. 3a). Conversely, genes involved in the biological process of dormancy were up-regulated at cold temperatures (< 10 C) in the southern *A. rudis* colony, but not *A. picea* (Fig. 3b). As these genes represent only a handful of the responsive genes, the expression patterns of all responsive genes can be searched and viewed interactively at [https://johnsg.shinyapps.io/ApRxN-shinyapp/].

![Figure 3. goes here]


Responsiveness ... a temperate ant species are thermally responsive. The majority of these genes (about 90%) are up-regulated at high or low temperatures. Discuss the 26 "stress response" genes

Previous studies have using *Drosophila* [@sorensen2005] and a temperate flesh fly [@teets2012-1] have implicated GO terms .... as responding to heat stress. In this study, we found .... of these GO terms to be enriched at high temperatures for *Aphaenogaster*. However, when we looked at the exact identify of the genes, we only found X overlapping genes.

This result poses a paradox. Why is the heat shock functional response highly conserved at the biological process level, while the exact molecular pathways appear to be highly labile? We suggest two possible explanations. First, regulatory networks controlling gene expression are evolutionarily dynamic. Trade-offs between function in different organisms may constrain which genes are involved in thermal response, thus resulting in different pathways to the same phenotype. Second, variation in molecular pathways may be due to inconsistent selection on thermal stress response. That is, as the environment changes and species shift their distributions, selection may relax or strengthen repeatedly within and across lineages. Different genes may be recruited within each selective episode. Differentiating between these, and alternate, explanations will require greater geographic and population genetic sampling.

A key outcome of this study was developing the first genomic tools for the study of adaptive differentiation in *Aphaenogaster*. Our comparative results are limited to two populations, and an important future step will be expanding this work to examine differences in gene expression among a larger number of populations. Moreover, the potential for species to adapt to climate change will depend in part on their genetic variation for thermal tolerance. The quantitative genetic studies required to estimate this adaptive potential have been out of reach of ant biologists due to the long lifespan of colonies, and inability of ants to mate in the lab. However, future work could leverage genomic information to infer the relatedness among colonies in the absence of a known pedigree, enabling estimation of genetic variation within populations and their adaptive potential. 

In summary, we have generated the first transcriptome of a temperate ant species, and have identified approximately 10% of genes that show thermal responsiveness. In this group, we find many known stress response genes, but also identify novel genetic pathways that respond to temperature. These results provide key insights on the genetic basis of thermal tolerance, and provide a tremendous resource for the future study of ecological adaptation in ant species.


## References

