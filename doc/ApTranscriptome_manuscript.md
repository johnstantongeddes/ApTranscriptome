---
output:
  pdf_document:
    biblatex: yes
biblio-files: ApTranscriptome.bib
---
Thermal reactionome of a temperate ant species
============================================

**Authors**: John Stanton-Geddes, Andrew Nguyen, Mahesh Vangala, Jim Vincent, Rob Dunn, Aaron Ellison, Nate Sanders, Nick Gotelli, Sara Helms-Cahan  

## Introduction

Temperature plays a prominent role in structuring biological diversity from molecular to macroecological scales [@kingsolver2009]. At the population level, thermal tolerance has been assessed in dozens of studies often demonstrating adaptive differentiation [@hoffmann2013]. In combination with environmental data, this thermal tolerance data has provided the backbone to dozens of studies examining limits to species geographic ranges [@angert2011; @diamond2013; @gouveia2014], and the impacts of climate change on population viability [@deutsch2008; @kingsolver2013]. At the organismal level, studies of thermal tolerance have largely focused on candidate genes such as heat shock proteins [@dahlgaard1998]. However, estimates of the standing genetic variation [@krebs1997; @krebs1998; @williams2012], studies of differential gene expression [@leemans2000; @sorensen2005; @teets2012], quantitative trait loci (QTL) [@morgan2006] and knock-outs [@takahashi2011] in the model system *Drosophila melanogaster* all implicate dozens of additional genes underlying thermal tolerance. Little is known about the diversity of these genes or their expression patterns, especially in non-model species.

The study of whole organism gene expression (e.g. transcriptomics) has rapidly increased due to improved sequencing technology (RNAseq). To date, RNAseq studies have largely focused on differential expression across few environments [@teets2012; @pespeni2013; @barshis2013; @fu2014] or tissues [@li2013] using an ANOVA-type statistical approach. However, many ecologically-relevant aspects of the environment such as temperature vary continuously in nature. Variation in the phenotype of an organism in response to such a continuous predictor (i.e. phenotypic plasticity) is referred to as a reaction norm. A reaction norm can differ not only in mean values, but also in the shape of the curve as is often found for thermal reaction norms [@hoffmann2013; @murren2014]. An analysis of phenotypic reaction norms found that differences in the shape were greater than differences in mean values in comparisons among species and populations [@murren2014]. The authors concluded that many patterns of biologically-relevant plasticity are likely to be missed in studies with few environments. Similarly, gene expression studies with few environments are likely to miss many biologically-relevant changes in gene expression. To overcome this limitation, we use a novel method to characterize the reaction norm of all genes in the transcriptome (e.g. the reactionome) in response to thermal stress. Specifically, we measure gene expression across a thermal gradient and identify thermally-responsive genes using a regression approach.

We apply this approach in a study of two closely related eastern North American ant species, *Aphaenogaster carolinensis* and *A. picea* [@umphrey1996], which have more southern and northern distributions, respectively. Ants are ecologically-important species involved in multiple symbioses¸ and are key contributors to seed dispersal [@lengyel2009; @rodriguezcabal2012], litter processing and decomposition [@lyford1963; @folgarait1998; @terborgh2001; @holec2006]. As ectothermic species, ants experience a wide range of temperature conditions. Throughout the year, ant species persist through freezings and hot summers, and at shorter time scales, individual workers can experience substantial thermal shock as they leave the buffered ant nest to forage for food. Thus, at the organismal level we expect that ants utilize multiple molecular pathways to tolerate thermal stress.

To explore the molecular toolkit used by *Aphaenogaster* to cope with thermal stress, we performed a whole organism RNA sequencing study. Our specific goals were to 1) generate the set of all genes (i.e. transcriptome) expressed by *Aphaenogaster* across physiologically-relevant temperatures, 2) identify genes that have responsive expression across this thermal gradient (i.e. thermal reactionome), 3) explore the functional annotations of genes that are activated or down-regulated in response to cold and high temperatures. We expected to observe candidate genes in the *response to stress* Gene Ontology (GO) category including *heat shock proteins* and ... in the thermally-responsive group.   

While characterizing the reactionome is useful for description of the molecular pathways for thermal tolerance, a greater understanding of the process of evolution can be gained using a comparative approach [@harvey1991]. Specifically, comparisons of populations or species in contrasting environments can reveal the extent to which geographic differentiation and selection has changed the number of genes involved in thermal tolerance. To accomplish this, we compared the thermal reactionomes between the two species of *Aphaenogaster*. *A. picea*, the northern species, experiences less chronic and acute heat stress than *A. carolensis*. Thus, we predicted that *A. picea* would have more genes than *A. carolinensis* that are up-regulated to high temperatures, while *A. carolinensis* would have higher constitutive expression of these thermally-tolerant genes. Reciprocally, we expected *A. carolinensis* to have greater up-regulation of genes in response to cold temperatures, while *A. picea* would have greater constitutive expression of cold-tolerant genes. Finally, phenotypic studies in *Aphaenogaster* [@warren2013] and other ecotherms [@hoffmann2013] have found less variation in critical maximum than minimum temperatures. From this, we hypothesized that if the species had more similar maximum temperatures than minimum temperatures, the northern species *A. picea* would have a broader thermal breadth of gene expression than *A. carolinensis*. 


## Methods

### Samples

Colonies of the ant species *Aphaenogaster picea* were collected from Molly Bog, Vermont (University of Vermont Natural Areas; 44.508°N, -72.702°W) and *Aphaenogaster carolinensis* from Durham, North Carolina ( 36.037° N -78.874° W), designated as *ApVT* and *AcNC* respectively. The phylogeny of *Aphaeogaster* is currently unresolved [@umphrey1996], but putative species identification was possible based on morphology (Bernice Bacon DeMarco, pers. communication). Colonies of these species were maintained in the lab for 6 months prior to experimentation. 

From each species, we haphazardly collected 12 ants at the same time on 12 days. Each day, the 12 ants were placed in glass tubes into a water bath at one of 12 randomly-assigned temperatures, every 3.5°C between 0° and 38.5°C, for one hour. The minimum and maximum temperatures were selected based on previous work showing that these temperatures are close to the critical minimum and maximum temperatures for *Aphaenogaster* [@warren2013]. At the end of the hour, the ants were flash frozen in liquid nitrogen and stored at -80°C. We extracted mRNA by homogenizing three pooled ants in 500 uL of RNAzol buffer with zirconium silicate beads in a Bullet Blender (Next Advance; Averill Park, NY), followed by RNAzol extraction (Molecular Research Center Inc; Cincinnati, OH) and then an RNeasy micro extraction (Qiagen Inc; Valencia, CA) following the manufacturer's instructions. 

### Sequencing, assembly and annotation

For each species, the 12 samples were barcoded and sequenced in a single lane of 2 x 100bp paired-end reads on an Illumina HiSeq 1500 yielding 200 and 160 million reads for the *A. picea* and *A. carolinensis* samples respectively. Reads were filtered to remove Illumina adapter sequences and low quality bases using the program Trimmomatic [@lohse2012]. 

We assembled the sequenced reads into the full set of mRNA transcripts, the transcriptome, for the combined data set from both species using the Trinity *de novo* transcriptome assembly program [@grabherr2011]. *De novo* transcriptome assembly is prone to falsely identifying alternative transcripts and identifying inaccurate transcripts that are chimeric (e.g. regions of two separate transcripts that assemble into a false, or chimeric, third transcript) [@yang2013]. We removed potentially false transcripts by first running the program `CAP3` [@huang1999] to cluster sequences with greater than 90% similarity and merge transcripts with overlaps longer than 100 bp and 98% similar in length. Second, we ran the program `uclust` which clusters sequences completely contained within longer sequences at greater than 90% similarity (see Supplemental Information for exact parameters). We used liberal values (90% similarity) to merge orthologous transcripts in the two species that may not have assembled together in the initial *de novo* transcriptome assembly. To identify contaminant sequences, we screened our full transcriptome using the program [DeconSeq](http://deconseq.sourceforge.net/) [@schmeider2011] with the provided bacteria, virus, archae and human [databases of contaminants](ftp://edwards.sdsu.edu:7009/deconseq/db). 

To determine the putative function of the transcripts, we performed functional annotation of the transcriptome assembly using the web-based tool `FastAnnotator` [@chen2012] which annotates and classifies genes by Gene Ontology (GO) term assignment, enzyme identification and domain identification. 

We performed subsequent analyses at the transcript level due to limitations of working with a system lacking a reference genome, but refer to these transcripts as genes for simplicity.

### Thermally-responsive genes

We quantified expression of each gene using the program `Sailfish` [@patro2013] and use the bias-corrected transcripts per million (TPM) [@wagner2012] as our measure of gene expression. We included the contaminant genes identified by DeconSeq at the quantification stage to avoid incorrectly assigning reads to other genes, but removed these from further analyses. As preliminary examination of the data (see supplemental...) indicated that 7°C samples may have been mis-labeled, we omitted these data from the analysis. The expression values were highly correlated between species at each temperature treatment (r^2^ > 0.98) indicating that assembling the transcriptome with data from both species was justified (Supplemental).

To identify genes that have significant changes in expression across the thermal gradient, we fit the regression model

$$ log(TPM + 1) = \beta_0 + \beta_1(species) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(species * temp) + \beta_5(species * temp^2) + \epsilon $$
    
to each gene. We used $log(TPM + 1)$ as the response to control for skew in the expression data. For a continuous predictor such as temperature, this regression approach is preferred to an ANOVA approach as it can reveal non-linear responses such as hump-shaped or threshold effects even [@cottingham2005]. This method is robust to over-dispersion as we expect errors in the read count distribution [@anders2010] to be independent with respect to temperature. To correct for multiple testing, we calculated adjusted *P* values using the False Discovery Rate (FDR) approach of Benjamini and Hochberg [-@benjamini1995], retaining genes under the 5% FDR threshold.

We identified thermally-responsive genes as the subset under our FDR threshold that had significant $\beta_2(temp)$, $\beta_3(temp^2)$, $\beta_4(species * temp)$ or $\beta_5(species * temp^2)$ terms after step-wise model selection by AIC. All analyses were performed in R [@r_core_team2013] and are fully reproducible (Supplemental file). 

For each thermally-responsive gene, we predicted expression levels across the tested thermal range using the final linear model separately for each species. We used the predicted gene expression values to partition the thermally-responsive genes into five expression categories for each species; genes that had greatest expression at high (> 31°) temperatures (*High*), low (< 10°C) temperatures (*Low*), intermediate (10 - 30°C) temperatures (*Intermediate*), both high and low temperatures (*Bimodal*) or were not thermally-responsive in that species (*NotResp*). For the *Bimodal* group, we required that expression at both low and high temperatures was at least one standard deviation greater than the expression at the mean temperature of 19.5°C.

To describe the overall molecular toolkit for thermal tolerance in *Aphaenogaster*, we performed gene set enrichment analysis (GSEA) on the intersection of each expression category in the two species using the program `topGO` [@alexa2006; @alexa2010] with the `parentChild` algorithm [@grossmann2007] implemented in R [@r_core_team2013]. Briefly, this approach identifies GO terms that are overrepresented in the significant genes relative to all GO terms in the transcriptome, accounting for dependencies among the GO terms. To reduce redundancy in the enriched GO terms, we used the web program [ReviGO](http://revigo.irb.hr/).


### Species-level comparisons

To gain greater insight on the genetic basis of thermal tolerance in *Aphaenogaster*, we performed comparative analyses between the two species focusing on the thermally-responsive genes that had a significant species x temperature interaction. For each gene, we predicted expression levels across temperatures for each species using the full linear model, and then grouped genes into the four responsive categories, *High*, *Low*, *Intermediate*, and *Bimodal*, as well as a fifth category *Not Responsive*. To focus on differences in gene expression between the two species, we performed GSEA for the set difference of each category.

With the same set of genes, we explored the extent to which differences in the thermal reactionome between the species were due to changes in the mean or shape of the reaction norm for each individual gene. Following the example of Murren et al [-@murren2014], we tested for overall differences in the (a) mean expression value across temperatures (*M*), (b) slope of expression (*S*), (c) curvature of expression (*C*) and (d) wiggle of expression (*W*), which consists of all higher-order differences in shape not capture by the first three measures (formulas given in Appendix). As genes differed by five orders of magnitude in expression levels, we calculated the mean-standardized difference between the species $\Delta M = (Mean_{A. picea} - Mean_{Ar}) / \bar M$, where $\bar M = (Mean_{A. picea} + Mean_{Ar}) / 2$ yielding $\Delta M$, $\Delta S$, $\Delta C$, and $\Delta W$ for the mean, slope, curvature and wiggle of each gene, respectively. Using the mean-standardized values, we performed one-sample two-sided *t*-tests to determine if the overall mean, slope or curvature were greater in *A. picea* than *A. carolinensis*. 

Further, to quantify the relative contribution of changes in the mean, slope, curvature or wiggle on the overall differences in reaction norms, we defined the *Total* difference in the reaction norm for each gene as $\Delta T = \Delta M + \Delta S + \Delta C + \Delta W$. We then partitioned the contribution of each measure by dividing it by $\Delta T$. 

To test if genes that were responsive to cold temperatures in the southern *A. carolinensis* species or hot temperatures in the northern *A. picea* species had higher constitutive expression in the contrasting species, we performed non-parametric two-sided Wilcoxon tests comparing expression values near the optimum of 19.5° C.

As critical maximum temperatures are more constrained than critical minimum temperatures [@warren2013; @hoffman2013], we predicted that the more northern *A. picea* species would have a broader thermal range of stable gene expression than *A. carolinensis*. To test this with our data, for each *Intermediate* expressed gene in each species, we randomly sampled 1,000 temperature values weighted by the expression function (e.g. point of maximum expression had highest probability of being sampled). We then calculated the standard deviation of this random draw, and performed a *t*-test comparing the standard deviations of expression between species. Conversely, we predicted that the thermal sensitivity of *Bimodal* expressed genes would be greater in *A. carolinensis* than *A. picea*. We performed the same analysis as with the *Intermediate* transcripts, but with low standard deviation being indicative of a relatively un-sensitive gene (nearly flat), while the maximum standard deviation of expression is for a highly-convex expression pattern. Thus, finding one species had overall greater standard deviation of expression for *Bimodal* genes would be indicative of greater thermal sensitivity. 

### Data availability

The raw Illumina reads are available at [...], the script for analysis and version history is available at (http://github.com/johnstantongeddes/ApTranscriptome).  



## Results




The Trinity *de novo* transcriptome assembly included 126,172 genes with a total length of 100 million bp. Filtering to remove redundant or chimeric reads resulted in an assembly with 96,253 contigs and a total of 105,536 genes. For all genes, the total length was 63 million bp with an N~50~ length of 895 bp and a mean contig size of 593 bp. Of the 105,536 filtered genes, 55,432 had hits to the NCBI-nr database. Of these, 38,711 genes mapped to GO terms, 1,659 genes were identified to an enzyme and 18,935 genes mapped to a domain with >50% coverage. Of these, 5,787 genes have best hits to putative genes in the genome sequence of the ant species *Solenopsis invicta* and 47% have best hits to genes from the class Insecta.

### Thermally-responsive genes

We identified 22,495 genes with overall model fits of P < 0.05, retaining 9,809 after correcting for FDR. Of these transcripts, 8,044 significantly change expression with temperature including 6,727 that have a species x temperature interaction. The number of transcripts that increased expression with temperature (*High*), decreased expression with temperature (*Low*), increased expression at both high and low (*Bimodal*) temperatures, decreased expression at high and low temperatures (*Intermediate*) or were not expressed differed significantly between the species ($\chi^2 = 170, P < 0.0001$; Table 1). For both species, the majority of responsive transcripts were in the *Low* category (Fig. 1). More transcripts had *Intermediate* expression in *A. carolinensis* than in *A. picea*, which had nearly twice as many bimodally-expressed transcripts  (Fig. 1, Table 1).  


```
## Error: object 'tt2' not found
```

```
## Error: object 'tt2' not found
```


As the molecular pathways involved in thermal response to cold and hot temperatures may differ, we performed gene set enrichment analysis on each expression category separately. Out of a total of 2,859 GO terms scored, we found 6 enriched (*P* < 0.01) for *High* transcripts, included categories implicated in stress response such as "gene expression" and "tissue regeneration" (Table S1). For *Low* genes, 60 GO terms were enriched including some implicated in stress response (e.g. detection of external stimulus; Table S1). For *Bimodal* genes, 11 GO terms were enriched including 10 "regulation" terms (e.g. regulation of primary metabolic process) indicating that these are general genes involved in activating responses to change (Table S1). Finally, for *Intermediate* genes 6 GO terms were enriched, all some type of metabolic or biosynthetic process (Table S1) consistent with the expectation that these are non-essential genes expressed in optimal conditions that are shut-down under thermal stress.

![Probability density function of peak expression for all thermally-responsive in *A. carolinensis* (blue) and *A. picea* (red) species](../results/PDF_expression_all.png)

As these genes represent only a handful of the responsive genes, the expression patterns of all responsive genes can be searched and viewed interactively at [https://johnsg.shinyapps.io/ApRxN-shinyapp/].

In addition, we examined the reaction norms for genes in the GO category "response to stress" which we *a priori* expected to have significant thermal responses....


### Species-level comparisons

We quantified differences in the reaction norms of thermally-responsive genes between species as changes to due shifts in mean ($\Delta M$), slope ($\Delta S$), curvature ($\Delta C$) and wiggle ($\Delta W$). Overall, we found that there was no significant difference in position of the mean (*t* = -0.95, *P* = 0.34), but changes in the shape due to slope (*t* = 24.67, *P* < 0.001), curvature (*t* = 15.71, *P* < 0.001) and wiggle (*t* = -12.45, *P* < 0.001) were all highly significant. On average, slope (95% CI 0.41 - 0.48) and curvature (95% CI 0.25 - 0.33) were greater in *A. picea* than *A. carolinensis*, indicating a steeper response to temperature. 

Of the total differences ($\Delta T$) in reaction norms between species for each gene, about 8% (95% CI 1 - 32%) was due to shifts in mean, 47% (95% CI 14 - 92%) due to shifts in slope, 39% (95% CI 0 - 67%) due to shifts in curvature and 0% (95% CI 0 - 87%) due to shifts in wiggle. 

Plasticity is a well-characterized mechanism by which organism tolerate environmental change. Given our sampling of a warmer (*A. carolinensis*) and colder (*A. picea*) climate species, we hypothesized that genes which are plastically-expressed under heat shock in the colder species would have greater constitutive expression in the warmer species, and vice versa. Consistent with one of our predictions, we found that constitutive expression was greater in the more southern *A. carolinensis* species for *High* genes in *A. picea* (Wilcoxon V = 1.721 &times; 10<sup>5</sup>, *P* < 0.001, 95% CI location shift -0.105, -0.038). However, for genes with *Low* expression in *A. carolinensis* we also found that constitutive expression was greater in *A. carolinensis* than *A. picea* (Wilcoxon V = 2.041 &times; 10<sup>6</sup>, *P* < 0.001, 95% CI location shift -0.229, -0.148).

We tested for a transcriptome-wide pattern of thermal tolerance and sensitivity by examining the thermal range) for the *Intermediate* expressed genes. Contrary to our expectations based on critical maximum and minimal temperatures, *A. carolinensis* had a transriptome-wide pattern of broader thermal range (*t*~df=1406~ = 9.25, *P* < 0.001, 95% CI difference in means 0.41 - 0.64). In contrast, there was no difference in thermal sensitivity as measured by the standard deviation of expression for *Bimodal* expressed genes (*t*~df=1406~ = 0.31, *P* = 0.76). 

We performed gene set enrichment analysis on the set difference of each expression category to reveal differences among the species response to thermal shock. For *High* genes, there were 4 enriched GO terms in *A. picea* and 12 in *A. carolinensis*. For *Low* genes, there were 40 enriched GO terms in *A. picea* and 23 in *A. carolinensis*. For *Intermediate* genes, there were 3 enriched GO terms in *A. picea* and 54 in *A. carolinensis*. For *Bimodal* genes, there were 12 enriched GO terms in *A. picea* and 10 in *A. carolinensis*. 

#![Expression (scaled to allow comparisons among genes) for genes in GO category "response to stress"](../results/GO_stress_High.png)


## Discussion

Thermal tolerance is a key trait for determining the survival of individual organisms, the persistence of populations and the geographic limits of species. While considerable research has focused on a few key genes (e.g. heat shock proteins) involved in thermal stress response [...], little is known about the functional diversity or transcriptome-wide patterns of thermal responsiveness. In this study, we characterized the thermal reactionome of sister ant species, *Aphaenogaster picea* and *Aphaenogaster carolinensis*. We found that nearly 10% of the genes in this species show significant thermal responsiveness.

### Molecular toolkit of thermal responsiveness

While the sheer number of thermally-responsive transcripts makes generalizations difficult, groupings such as Gene Ontology (GO) can be used to identify categories associated with thermal stress. For heat response, the semantic space of the enriched GO terms included *gene expression*, *regeneration* and *cellular component organization or biogenesis*. ...previous work... For cold response, there were 65 enriched GO terms that grouped by semantic space into *protein modification process*, *sensory perception of mechanical stimulus*, *response to mechanical stimulus*, *folic acid metabolism*, *membrane lipid metabolism* and *protein complex biogenesis*. ...previous work...

*Bimodal* genes are general stress response genes that are activated under both heat and cold stress. The semantic space of the enriched GO terms included *regulation of metabolism*, *biological regulation*, *macromolecule biosynthesis*, ...

*Intermediate* genes represent non-essential molecular functions that are shut-down as conditions become stressful. The semantic space for the enriched GO terms included *lipid metabolic process*, *glycoprotein biosynthetic process*, *membrane lipid metabolic process*, *cellular protein metabolic process*, *membrane lipid biosynthetic process*. ...indicates shutting down of non essential molecule biosynthesis...

**Sara - comparisons of our gene list to other studies**

This result poses a paradox. Why is the heat shock functional response highly conserved at the biological process level, while the exact molecular pathways appear to be highly labile? We suggest two possible explanations. First, regulatory networks controlling gene expression are evolutionarily dynamic. Trade-offs between function in different organisms may constrain which genes are involved in thermal response, thus resulting in different pathways to the same phenotype. Second, variation in molecular pathways may be due to inconsistent selection on thermal stress response. That is, as the environment changes and species shift their distributions, selection may relax or strengthen repeatedly within and across lineages. Different genes may be recruited within each selective episode. Differentiating between these, and alternate, explanations will require greater geographic and population genetic sampling.

### Transcriptome-wide comparisions

Of the 8,492 genes that were thermally-responsive, 80% differed in their responsiveness by colony. The majority of this difference in the reactionomes was due to changes in the slope ($\Delta S$) and curvature ($\Delta C$) which together accounted for nearly 90% of the differences between reaction norms for all transcripts. In a meta-analysis of phenotypic reaction norms, Murren et al. [-@murren2014] similarly found that differences in the shape of the reaction norms accounted for the majority of the differences. Thus, our finding suggests that this phenotypic pattern may be explained by greater genetic lability in the shape than the mean of gene expression. The relative role of plasticity versus adaptation in response to rapid environmental change is widely debated [@vellend...?]. Our data suggest that plasticity through changes in the responsiveness of gene expression are likely to occur more rapidly than changes in mean values. 

An important note is that we would have completely missed this pattern if we had performed a differential expression study using only a few temperatures. Only by measuring gene expression across a thermal gradient and fitting a function to each transcript were we able to characterize the shape of the reaction norm of expression... To our knowledge, this is the first time that a reaction norm approach has been applied to the study of global gene expression. 

Given the large differences in number (Table 1) and responsiveness of genes between the two species, a striking pattern is that 3-4x more genes are upregulated at cold than high temperatures. Intriguingly, this result mirrors the observation in *Aphaenogaster* [@warren2013] and other ectotherms [@hoffmann2013] that critical maximum temperatures are more constrained than critical minimum temperatures. A plausible explanation is that organisms have fewer molecular mechanisms to tolerate heat than cold stress. Future studies ... implications for climate change as the evolutionary potential of ecototherms to climate change may be limited as suggested by ... [kingsolver?]

### Conclusions 

A key outcome of this study was developing the first genomic tools for the study of adaptive differentiation in *Aphaenogaster*. Our comparative results are limited to two populations, and an important future step will be expanding this work to examine differences in gene expression among a larger number of populations. Moreover, the potential for species to adapt to climate change will depend in part on their genetic variation for thermal tolerance. The quantitative genetic studies required to estimate this adaptive potential have been out of reach of ant biologists due to the long lifespan of species, and inability of ants to mate in the lab. However, future work could leverage genomic information to infer the relatedness among species in the absence of a known pedigree, enabling estimation of genetic variation within populations and their adaptive potential. 

In summary, we have generated the first transcriptome of a temperate ant species, and have identified approximately 10% of genes that show thermal responsiveness. In this group, we find many known stress response genes, but also identify novel genetic pathways that respond to temperature. These results provide key insights on the genetic basis of thermal tolerance, and provide a tremendous resource for the future study of ecological adaptation in ant species.


## References

