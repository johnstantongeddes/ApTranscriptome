Dear Dr. Hoffmann,

We appreciate the comments from you and the reviewers. Both reviewers raise important points and concerns, to which we have responded below with changes to the manuscript and/or rebuttal if we disagreed with the reviewer’s suggestion. We feel that this process has helped us to build a stronger case both for the methodology and main results of our study. We especially thank the reviewers for highlighting ambiguities in some of our original language.  While we respond to the reviewers comments in details below, we address the three major concerns here. 


1) Acclimation vs adaptation

Reviewer 1 is concerned that we “cannot tease apart genetic or plastic differences in expression profiles” due to using field-collected ant colonies in this experiment. While a complete-generation common garden experiment is preferable for distinguishing genetic from plastic differences, this is neither possible nor practical for the vast majority of species, including Aphaenogaster. Moreover, rearing a single generation in common conditions does not resolve all these issues (e.g. maternal genetic effects). A key outcome of our work is demonstrating a method to gain insight on the relative role of genetic or plastic differences even when rearing organisms in common gardens is not feasible. We limited differences in developmental acclimation by maintaining colonies in the lab for 6 months prior to sampling. Given that a new cohort is reared in 2-3 months, 6 months allowed for a substantial number of workers to be fully reared under common-garden conditions, such that our samples would have contained a mix of truly common-garden and adult-acclimated workers. This is actually better than has been ascribed to adaptation and published for other long-lived taxa including corals (Barshis et al. 2013 PNAS 110:1387) and fish (Smith et al. BMC Genomics 2013, 14:375). With hindsight, the reviewer’s suggestion of sampling the species from sites where they overlap is great, but unfortunately was not done. However, it is worth noting that in previous work where the species were collected from overlapping sites along an altitudinal gradient in the Smokies (Warren and Chick 2013), both CTmax and CTmin differed significantly, suggesting different underlying physiological adaptations such as we have described here.

We have addressed this concern, noting the points above, in the section “Common Garden Design” in the Methods (lns 165 - 177).

2) Replication of RNA-seq

The standard practice of replication in RNA-seq experiments was developed for ANOVA to identify differences in expression among tissues and discrete conditions (e.g. disease states). For these studies, replication is absolutely essential. In contrast, ecologists are often interested in the response across a gradient, which is more directly modeled by a linear regression analysis (Cottingham et al. 2005). Variance in the expression profiles is estimated from the linear regression model. Replication at any particular temperature is not necessary, as power to detect a relationship is most improved by adding additional temperatures, not by replicating those that have already been assessed.  In fact, this approach results in a more precise variance estimate than if we had replicated only a few temperatures (Cottingham et al 2005). Only transcripts with a significant linear or quadratic relationship with temperature, after correcting for false discovery rate, were retained in the analysis.  Such analyses are just beginning to be used in a transcriptomic context (the only other example we are aware of is Hodgins-Davis et al 2012 GBE). To better clarify our approach for readers, we have added an expanded justification for the regression design in the Introduction  (lns 113-116):

To characterize the thermal reactionome, we characterized the reaction norm for each gene using a regression-based statistical approach to identify temperature-dependent patterns of change in gene expression. We used these response patterns to quantitatively test three mechanistic hypotheses of thermal adaptation.

and explanation in the methods  (lns 178-183):

	Unlike ANOVA-based experimental designs, which derive statistical power from replication within each experimental treatment level, regression designs have greater power when sampling additional values across the range of the continuous predictor variable (Cottingham et al 2005). Thus, we focused our sequencing efforts on maximizing the number of temperatures at which the transcriptome was profiled, rather than on replication at each temperature.


2) Validation with qPCR

This is an important issue to which we gave careful consideration. For two reasons, we did not feel that the results of qPCR would improve this study. 

First, RNA-seq is no longer a new methodology. Many studies have been published in ecological journals using RNAseq data without qPCR validation, such as Barshis et al. (2012 PNAS) and Smith et al. (2013 BMC Genomics).

Second, it is only feasible to perform qPCR on a small subset of transcripts. In our approach, we avoid focusing on the response of individual transcripts (e.g we *do not* provide a list of differentially-expressed genes). Instead, we focus on transcriptome-wide responses to temperature, which groups large numbers of genes into distinct response patterns. While our results are undoubtedly wrong for some individual transcripts, they are likely to be distributed at random in regards to these classifications. Thus, while qPCR may prove that the expression patterns of some transcripts are incorrect, it is unlikely to change the results for whole functional groups.





> Editor Comments to Author:
> 
> This interesting paper has received three careful reviews. The first reviewer recommended rejection, the second was enthusiastic about the paper, and the third sat somewhere in the middle. The first reviewer points out that the work does not tightly control for genetic versus plastic effects and therefore the paper does not precisely test the issues of interest to the authors. I suspect that this can be fixed with reworking of the paper in terms of interpretations, taking care to separate out the different components where appropriate. The reviewer also makes some important points about replication and the measurement of CTmax as well as the nature of control thermal treatments. The second reviewer was complimentary about the overall approach and had only minor comments, which cover some confusion stemming from the text and one of the Figures. Some apparent errors in the wording are also noted. The third reviewer also raises the issue of replication at the different temperatures and has an issue with the terminology in the paper and division into a predefined set of categories. The comments are worth taking on board. The reviewer also points out that some aspects of expression changes are likely to be missed, and calls for new qPCR data as is often carried out. The reviewer points out where results are inconsistent with the literature and novel patterns that appear not to be consistent with expectations, which all need discussion. The authors should take these comments on board when making interpretations of their data, and when interpreting the patterns established from their work. The issue of replication needs consideration along with measurement error and validation of expression patterns established from RNAseq.
> 
> Dr. Ary Hoffmann
> Subject Editor, Molecular Ecology


> Reviewer Comments to Author:
> Reviewer: 1
> 
> Recommendation: Reject
> 
> Comments:
> 
> The current paper by Stanton-Geddes et al generates “reactionomes” by generating RNA seq data  for two different species of ant at 12 different temperatures (0-38.5C).  The authors estimate CTmax in both species and link changes in their expression profiles to these estimates of CTmax and the environments where the species were collected. I think the authors were ambitious, the questions they present relevant and highly topical however, I believe there are significant design issues which limit the interpretation of the results.
> 
> First and fore mostly, if I understand the methods correctly, ants were collected from the field acclimated at 25C for 6 months and then experiments performed. While I understand that many species can be difficult to rear in the laboratory. Without performing a true common garden experiment (i.e both species of ants developed and reared under standard conditions) you simply cannot tease apart genetic or plastic differences in expression profiles. Are the observed differences between the species due to evolved differences related to thermal evolution or simply a consequence of developmental acclimation, by collecting cold adapted species from a cold environment and vice versa only confounds the results even more. One possible way that the authors may have alleviated this issue, to some degree, is to choose 2 species that vary in their distribution but where their distribution overlaps. This may be the case for these two species but the authors have not made that clear in their ms. If their species do overlap in distribution the authors should at the very least have collected these species from the same site. With such a design you could have then assumed that developmental acclimation may have been similar across the two species. Albeit you would still have GxE effects and these are likely to differ between species.

We address this concern as (1) above. 


> My other major concern, which also depends critically on whether I have understood the methods correctly (see below comments) is that for the RNA seq the authors have pooled 3 ants per treatment (not a particular high number) and actually have no replicates per treatment. How do you estimate the variance in your expression profiles?


We address this concern as (2) above.




> Minor comments
> 
> -Although authors have characterized CTmax in the two species, how do these estimates of CTmax relate to the heat treatments that the authors use? Have the authors really examined a phenotype here? The authors CTmax estimates will also be confounded by developmental acclimation.
> 

The maximum temperature we heated ants to was 38.5℃, below the CTmax for both A. carolinensis (~43℃) and A. picea (~42℃) (Fig 3. Warren & Chick 2013), so we have confidence that gene expression changes reflected active responses to thermal stress rather than dysregulation caused by loss of function. We have added this information to the “Common Garden Design” sub-section of the Methods  (lns 191 - 194).

Both colonies were maintained in the lab for 6 months prior to experimentation, so differences in developmental acclimation are as minimized as the biology of the species will allow, but it is of course true that we cannot rule out some residual developmental effect. To reflect this remaining potential source of phenotypic plasticity, we have attempted throughout the ms to make our word choices more precise in regard to whether we have demonstrated evidence for evolutionary differences.

> -The authors need to clearer about where in the distribution of the ant species they are sampling. Are they sampling them from their extremes? Middle?

We have clarified this in the Methods  (lns 154 - 155). 

    These locations are centrally located within each species’ geographic range.

> 
> -Line 40-55 The end of this paragraph reads like methods
> 
> -Line 56-60 also reads like methods
> 

We appreciate the reviewer’s concern, and removed some of the detailed information (lines 56-57 from draft 1), but we do feel feel that some methodological details are necessary here to explain how we tested our hypotheses.  We have attempted to reframe these such that the reason for their inclusion is more clear.

> -Line 57 it would be helpful to have a figure of the temperatures you used here, or at least have it written out in the methods
> 
> -I found the methods section confusing and not detailed enough
> 
> -Line 97 makes it sound like you collected each ant species form Molly Bog Vermont, but I don’t think that’s  what you mean I think you mean you collected a single colony of A. picea from Molly bog and a single population of A. carolinensis from Durham North Carolina

True - we have clarified this in the text  (lns 150 - 153):

	In fall 2012, we collected a single colony of Aphaenogaster picea from Molly Bog, Vermont (University of Vermont Natural Areas; 44.508° N, -72.702° W) and a single colony of Aphaenogaster carolinensis, part of the A. rudis species complex (Umphrey 1996), from Durham, North Carolina (36.037° N -78.874° W). 

> 
> -Why were you not able to determine the CTmax from the original, 

It would have been ideal to have paired measures such as these, but both transcriptome and CTmax measurements require a large number of workers, and the colonies were simply not large enough to do both. We’ve added this information to ln 157 of the manuscript.

could you please describe your CTmax methods, there are many different ways this is calculated. It shouldn’t be necessary for me to have to look up another paper to understand your methods.

We have added this information as requested  (lns 160 - 163):

	We tested the upper and lower critical thermal limits for 5 ants from each of these colonies using a ramp of 1° C per minute, starting at 30° C, and recorded the temperature at which the ants were no longer able to right themselves (Warren and Chick 2013).


> 
> -Line 109 To be clearer could the authors state that this part of the experiment was performed on ants collected in 2012

Done  (ln 150).

> 
> -The sampling and replicates are not clear Line 109-120. The authors collected 10 ants from each colony and exposed them to a single random temperature on a given day for a period of 12 days. Ants were then pooled into 3 which would make me think you had 10 ants per species per treatment you pooled them into 3 to make almost 4 replicates but on line 122 you say for each species the 12 samples were barcoded?

We apologize for the confusion. We collected 10 ants so that we would have extra samples if necessary, though only 3 were used in pooling (per temperature) for sequencing. For clarity, we now simply refer to 3 ants. 

> 
> -Did the authors ensure these treatments did not cause mortality?

Yes. We have added this to the text (lns 191 - 193).

> 
> -Should the authors have a control here i.e 24.5 and look at species specific expression changes?
> 

It’s unclear what the reviewer means here. A major advantage of the regression approach is that we have quantitative information on transcript expression in response to temperature and can predict expression at 24.5, or any other temperature, from the model of each transcripts’ reaction norm.  Species-specific expression changes are identified by a significant species x temperature interaction.  We have adopted the reviewer’s suggestion below to use 24.5 as our baseline value (see response to comment for line 222).

> -Line 162 I am unsure what the authors are suggesting here, that the data was odd so the got rid of it?

As you can see in this figure, which we now include in the Supporting Information, the “Ar_7” (Ar = A. carolinensis) sample clusters with all the other A22 (= A. picea) samples, while the “A22_7” sample clusters with all the other Ar samples. The misplaced samples are marked with red arrows. The most parsimonious explanation for this is a labeling mistake. Removing these samples, rather than relabeling, is the conservative approach. 




> 
> -Line 169 Rather than writing out formula could you at least include a description of your model. Where species considered random or fixed effects? You can run mixed model glm’s in R and maybe this is what you did but your description is unclear. It is also unclear to me what temp and temp2 are in the model, this is not explained anywhere .
> 

The formatting of this equation was lost in the file conversion, temp2 should have been temp2 indicating that we were fitting a second-degree polynomial regression. All effects were fixed as temperature is continuous and we only had two species. We have clarified in the text (lns 250 - 255):

	Temperature and species were both fixed effects, with a quadratic term included for temperature.


> -Line 180-184 I am unsure what exactly you are performing the step wise regression on?
 
The statistical significance of the species and linear and quadratic temperature terms in the model were tested via step-wise model selection.

> -Line 188-192 What about down regulated transcripts?

In a reaction norm context, downregulation is the inverse of up-regulation. That is, transcripts downregulated at low temperatures fall into the *High* category, transcripts downregulated at high temperatures fall into the *Low* category.  Splitting the temperature treatments into two at the midpoint to consider only low or only high temperatures would potentially have allowed us to make distinctions between upregulation and downregulation in these cases, but with only 12 temperatures in total, we would unfortunately be left with too few temperatures in each case to estimate the model.  We have clarified this in the sub-section “Identification of thermally-responsive transcripts”  (lns 279 - 284)

Because expression category was defined by the temperature of maximal expression, both Low and High categories were biased toward transcripts up-regulated at that temperature extreme, but also likely included some transcripts down-regulated at the opposing extreme. The two categories which could unambiguously distinguish up- from down-regulation are Bimodal (up at both extremes) and Intermediate (down at both extremes).

> 
> -The authors grouped their transcripts by expression level and temperature, the temperature categories need better justification

The temperature categories were decided a priori to allow us to apply the results to the biological hypotheses we were testing. This is elaborated in the sub-section “Identification of thermally-responsive transcripts”  (lns 271 - 284):

	We used the predicted transcript expression levels to partition transcripts for each species into five expression categories which were defined a priori to allow us to test predictions derived from three thermal adaptation hypotheses of relative response severity in the two species...


> 
> -Line 222 It seems odd that the authors have used 19.5C as their baseline for constitutive expression when ants were reared for 6 months at 25C. It seem 24.5 C would be a more suitable control.

This is a good point that we had not considered.  As suggested by the reviewer, we changed the baseline for constitutive expression to the rearing temperature of 25℃. See commit afdd69d on github (https://github.com/johnstantongeddes/ApTranscriptome). We use this version of the results in the revised manuscript, though it has no qualitative effect on the major results. 

> 
> -Line 322 the authors start talking about down regulation but as I understood their categories high was high expression at high temperatures, low was high expression at low temperatures etc. Maybe I have misunderstood something here?

See our comment above.

> 
> -Line 424 It is unclear what the authors are trying to say here are they suggesting you can link the number of responsive genes to the level of additive genetic variation?
> 

We were presenting a post-hoc explanation for the observation that many more genes were responsive to cold than warm temperature extremes. As we emphasized in the text (e.g. To the extent that …) this was a speculative explanation for this observation, though consistent with other work (Kellerman et al 2012, Warren and Chick 2013). However, given the uncertainty about the explanation we have chose to remove this section from the manuscript.

> 
> 
> Reviewer: 2
> 
> The authors of this study characterize the thermal reaction norms of all genes in a de novo transcriptome assembly across 12 temperatures spanning the entire thermal breadth of two species of ants: the cool-climate Aphaenogaster picea and the warm-climate Aphaenogaster carolinensis. In the past, studies exploring transcriptomic responses to temperature have mostly used ANOVA-type experimental designs, with only a few temperature levels. But continuous variation is much better characterized with the reaction norm approach used in this study, describing variation among genes and among organisms in the shape of expression responses across temperature. I really enjoyed this paper -- the approach taken here is exciting, and I found the manuscript to be quite well written, with clearly defined hypotheses an expectations. The comparative reaction norm approach taken here is, to my knowledge, unprecedented, and will make an important contribution to the field of thermal biology. I strongly recommend this manuscript for publication (and look forward to its acceptance so I can cite it!). I have only a few minor comments and revisions, outlined below.
> 

We really appreciate these comments. While this study isn't perfect, we also believe that this study is pushing back the boundaries of our knowledge of thermal biology. 

> Lines 61-86. As the authors note, studies exploring transcriptomic responses to stress have sometimes been prone to storytelling, so I appreciate the specific predictions generated by each of the three potential mechanisms of stress adaptation.
> 

This problem of storytelling is rife in the ecological genetic literature, and one of our motivations for writing this paper was to demonstrate how similar methods can be analyzed in an alternate approach to gain novel insight.

> Figure 1
> 
> -Figure 1 is key to the whole paper, but connecting the results to the predictions generated by each of the three hypotheses made my head hurt a little. I wonder if there is a way to make it easier for the reader to draw inferences from this figure? One possibility would be an inset cartoon showing the results that could be expected in support of each of the competing hypotheses? Also maybe include little line cartoons along the axes depicting each of the norms of reaction.
> 

This is a useful suggestion and we have made a figure (Fig. 1) with illustrations of the expected results under the *enhanced response* and the *tolerance* hypotheses that can be compared to the true results. 

> The categories “High” “Low” “Intermediate” and “Bimodal” are defined well in the manuscript, but it’s always nice when key figures can stand alone without reference to the text.
> 

We define these in the caption to Table 1, and refer to this in Figure 1 rather than repeating. We could also switch this if preferred.

> - What are the units in the color scale?

The units are number of transcripts. We have updated the legend to clarify this:

	...the color of each cell (purple = excess, orange = deficit) represents the deviation of the observed from expected number of transcripts. 

> 
> 196-197. I’m missing something about the methods - lines 196-197 state that the analyses focused on transcripts with a significant temperature x species interaction but in Figure S1, 4/5 of the boxes on the diagonal are blue (= excess of transcripts in the same reaction norm category for both species). Shouldn’t these have been excluded by focusing only on transcripts where temperature effects differed between species?
> 

We apologize for being unclear in the text. For the analyses of the how the number of responsive transcripts differs between species (first paragraph under “Statistical Analyses”), we used all responsive transcripts, both those with a significant temperature x species interaction and those with a main effect of temperature only. The justification for including all thermally-responsive transcripts is that the enhanced response hypothesis explicitly predicts that the genes upregulated by the species that sees the stressor infrequently should retain that function in the other, while under the tolerance hypothesis they should not. 

In contrast, for the comparisons of the *temperature of gene activation* (second paragraph under “Statistical Analyses”), only transcripts with a significant temperature x species interaction are used as, by definition, there is no difference in expression levels between species 
We have changed the text (lns 292-295) to clarify this.

	To test whether the temperature at which thermally-responsive transcripts were activated differs between species, we identified the temperature at which there was the greatest change in expression for each transcript in each species, using only the transcripts with a significant species x temperature interaction.


> -207-216 –is t test appropriate here, given that data are categorical rather than continuous? (Expression is observed at 12 distinct temperatures).
>

In this case, we aren't testing expression across the 12 temperatures. For each transcript, we used the predicted expression values from the statistical model to identify the temperature where expression had the greatest change (the temperature of transcript activation). We then compared the vectors of temperature of activation (a continuous variable) between species using the t-test. We now explain this in the text (lns 297 - 299):

	We then performed a t-test to determine if the mean temperature of transcript activation differed between the two species for each group.


> -235-245 Something’s not adding up about the number of transcripts in each category 237-239 says, “Of the Formicidae hits, we identified 9,052 transcripts, or 8% of the total, with expression patterns that were thermally responsive” (implying that the 9052 thermally responsive transcripts all had hits to Formicidae) But then later in the paragraph, “The proportion of responsive transcripts is the same if we focus only on those transcript (4,618 transcripts out of 51,246).” (Implying that only half of the thermally responsive transcripts had blast hits). 
> 

This was our mistake as the phrase “Of the Formicidae hits..." was left in unintentionally from a previous draft. The percentage of responsive transcripts is out of the total.

> 414-418 “In both species, there are 3 -- 4 times more transcripts upregulated at low than high temperatures (Table 1). A potential explanation for this pattern is that selection has been stronger on upper than lower thermal limits, thereby reducing genetic variation at high temperatures (Kingsolver et al. 2013; Kelly et al. 2013).” I would argue that the data show reduced plasticity of gene expression at higher temperatures rather than reduced genetic variation (= variation among individuals).
>


As this section raised concern from multiple reviewers and is not essential to our main conclusions, we have removed it from the manuscript.


> Reviewer: 3
> 
> Recommendation: Accept, minor revisions
> 
> Comments:
> 
> MEC-14-1243- MS TITLE: Thermal reactionomes reveal adaptive responses to thermal extremes in warm and cool-climate ant species
> 
> This paper brings valuable contribution to our understanding of thermal adaptation of two Ant species living in contrasted thermal regions. The authors found that a significant part of the transcriptome was thermally responsive. For both species, more transcripts responded to low than high temperatures, a result that I find counterintuitive and unexpected. The authors also found differential thermal-responsive transcriptional responses between the two species suggesting different thermal adaptation (and putative underlying mechanisms).
> 
> The authors present a set of experiments with insects that were assessed with a range of thermal conditions (12), ranging from lower to upper thermal limits, as well as transcriptional changes in reaction to temperature exposure (RNA seq). The topic is timely with the increasing interest in the molecular and physiological basis of response to temperature/climate change. The nature of cold response/tolerance is highly complex which makes the effort to unravel it very difficult.
> 
> Experiments were carefully conducted under well controlled conditions. The paper is logically arranged and well written, in good English. The methods used in this study were appropriate and the experiments were correctly conducted which makes the data set reliable. I also particularly appreciated that this is an hypothesis-driven study.
> 
> In spite of this, I have a number comments and queries that needs to be answered or at least discussed (see below).
> 
> -Abstract and elsewhere in MS : “enhanced defense hypothesis” . I find the term “defense” a bit confusing as it often used for immune process (defense genes …). “Stronger induced response hypothesis” could be an alternative.
> 

We agree with this comment and have changed the ‘enhanced defense’ to the ‘enhanced response’ hypothesis throughout the paper.

> -L77: please provide factual examples for “shift in alternate suite of tolerance genes”

We have added citations to Neelakanta et al (2012), who showed that inducing expression of an antifreeze glycoprotein can increase the ability of fruit flies to tolerate cold, and Franssen et al. (2011), who showed that a southern population of Zosteria marina is more resilient to heat shock than a northern population (ln 131).

> 
> -L61-86: I am not totally convinced that the whole set of thermal-responsive genes will coordinately change in such way that it will fit to one of the three proposed mechanistic hypotheses of thermal adaptation (1 defense, 2 tolerance or 3 assimilation). I assume that transcriptional response will be more complex and gene-specific. Some genes will fit the first hypothesis, while other will fit the others. Please comment.
> 

The reviewer is almost certainly correct, and it is our goal to determine the relative proportion of the transcriptome that falls into each category. We state this in the final sentence of the paragraph (lns 138-139:

Of course, these hypotheses are not mutually exclusive. By taking a reactionome approach, we can quantify if, and under what conditions, these mechanisms contribute to thermal adaptation.

> -L107: why such a small number of individuals (5) were used to monitor thermal limits?
> 

A previous study (Warren and Chick 2013) used larger sample sizes (10 ants per colony) to demonstrate differences in the thermal limits of these species. Our measurement of thermal limits were simply a confirmation of their results. 

> -L109-120: Ten ants were exposed to 12 different temperatures for each species. I may have misunderstood but it is not clear to me if there were replicates for each temperature?

Please see responses to reviewer 1 on this topic.

> 
> -L109-115: Twelve temperatures were used to establish a reaction norm of transcript expression patterns. Expression was checked at a single time point : after one hour of temperature exposure. It is well known that many thermal-responsive genes show temporal patterns of expression peaking at various time after exposure (see Sorensen 2005 Cell Stress & Chaperones  10, 312–328). Here, interpretations are based on early responsive genes (1h of exposure), which might only explain a small fraction of molecular response to temperature. Instead of taking 12 temperatures, it could have been of interest to check temporal patterns, perhaps with less thermal conditions. In fact, were 12 temperatures truly necessary to establish a reaction norm? My point is that this study should acknowledge that a large suite of thermally-responsive transcripts are overlooked and cannot be captured by the experimental design.


The temporal dynamics of expression response to thermal response are clearly crucial. We acknowledge this complexity by specifying that our reactionome is focused on capturing a wide temperature gradient, at the expense of having to focus only on a single time-point in the temporal response. Specifically, after defining the reactionome in the Intro we state (lns 94-97:

	Although temporal patterns of transcriptional activity (e.g. fast- vs. slow- responding genes) are also important components of an organism's transcriptional response to environmental conditions (Sørensen et al. 2005), we focus here on the response of transcripts across conditions at the same time point.

And in the Methods (ln 197):

	Thus, our reactionome characterized early, but not late, responding genes. 

Finally, in the section on caveats in the Discussion (ln 554-556):

	… species may differ in gene expression along axes which we have not measured here, primarily temporal patterns of gene expression (Sorensen et al. 2005). 


> 
> -Microarray and NGS data are regularly validated by qPCR. I am surprised that authors did not include validation step (on 5-10 genes) in their experiments. Previous studies have shown extremely close correlations between qPCR and RNAseq data [1-4]. Ideally, the authors should re-validate their findings (potentially by qPCR) in a separate cohort of samples.

Please see our general comments above.

> 
> 1.Griffith M, Griffith OL, Mwenifumbo J, Goya R, Morrissy a S, et al. (2010) Alternative expression analysis by RNA sequencing. Nat Methods 7: 843–847. doi:10.1038/nmeth.1503.
> 
> 2.Asmann YW, Klee EW, Thompson EA, Perez E a, Middha S, et al. (2009) 3’ tag digital gene expression profiling of human brain and universal reference RNA using Illumina Genome Analyzer. BMC Genomics 10: 531. doi:10.1186/1471-2164-10-531.
> 
> 3.Wu AR, Neff NF, Kalisky T, Dalerba P, Treutlein B, et al. (2014) Quantitative assessment of single-cell RNA-sequencing methods. Nat Methods 11: 41–46. doi:10.1038/nmeth.2694.
> 
> 4.Shi Y, He M (2014) Differential gene expression identified by RNA-Seq and qPCR in two sizes of pearl oyster (Pinctada fucata). Gene 538: 313–322. doi:10.1016/j.gene.2014.01.031.
> 
> -L197: comparative analyses between species were performed by focusing on transcripts that had species x temperature interaction. What about genes that had unique profile, that is, genes that were species-specific? These might contribute to thermal adaptation. How these genes were treated?

The species-specific genes were included in the analysis. In Table 1, these are the genes that were expressed in one species, but were “Not Responsive” or missing in the other species. We discuss some of these genes in the paragraph in the Discussion that begins “In contrast to cold tolerance, the enhanced thermal limit in A. picea is best explained by the tolerance hypothesis…” (ln 469 - 489).

> 
> -L196-201 & 207-216. It wasn’t clear to me why authors arbitrary defined 5 expression patterns (categories) ? Authors could have use k-means algorithm that can be used to partition the input data set into k partitions (clusters). In such way, main groups of temperature-responsive genes could be identified by K-mean clustering without a priori, and this would also provide statistical measures of over-representation (EASE scores).
> 

The expression categories were not arbitrary, but represent a priori predictions that are biologically meaningful and can be related to the predictions of the three hypotheses (lns 272 - 275). While statistical clustering may be valuable, it would make hypothesis testing difficult. What would one do with a cluster that increased at 26 degrees, decreased at 29 and then leveled off, for example? This could be explored in another paper, but is not a good analytical fit with the goal here. Further, we note that Reviewer 2 found our focus on specific predictions explored by the specific expression patterns to be a strength of the paper.

> -L210. it seems that authors only focused on genes that were upregulated to fall into 5 categories. What about genes that showed downregulation in expression patterns ?
> 

Please see our response to reviewer 1 on this topic.

> -L258: it is rather surprising that majority of transcripts were upregulated at the low temperature, especially since there was no recovery to allow transcriptional machinery to operate. It has repeatedly observed that at low temperature, close to lower thermal limit, transcriptional activity is very weak. For instance, Teets et al. (Physiol. Genomics 44, 764–777, 2012) did not observe any differentially expressed transcript during the cold treatment and Vesala et al. (Insect Mol. Biol. 21, 107–118.2012) only found a few transcripts (3 or 0) modulated by cold in D. montana and D. virilis, respectively. Please comment.
> 

This is an interesting observation that actually strengthens our reactionome approach, as we can quantify not only what is upregulated but the temperature at which that response occurs. We now comment on this in the Discussion (lns 505 - 515):

	In both species, there were 3 -- 4 times more transcripts upregulated at low than high temperatures (Table 1). The degree of upregulation at low temperatures is surprising given previous studies (Vesala et al. 2012; Teets et al. 2012b) that found little transcriptional activity at low temperatures. However, these studies exposed organisms to a few extreme (-10 -- 0°C) temperatures. At these extreme low temperatures we also found few upregulated transcripts (Fig. 3A), whereas the peak of low-temperature transcriptional activation occurred near 10°C (Fig. 3). A potential explanation for this pattern is that increased gene expression functions to support elevated metabolism at moderately cold temperatures, as suggested by the metabolic cold adaptation hypothesis (Addo-Bediako et al. 2002). 

> -L362-379: A.picea is more cold-adapted than A. carolinensis. A.picea also shows more transcripts being expressed at cold than A. carolinensis. Hence, this correlated pattern is assumed to support enhanced defense hypothesis, that is, stress-adapted species show stronger response. This discussion is based on global expression patterns of hundreds of genes that show expression variation according to temperature. The metrics used here is the number of responsive transcripts. These global patterns are assumed to determine differential thermal adaptation between species. It is certainly interesting to look a general patterns and numbers of transcripts, but this bulk of transcripts obviously include a number of genes that have nothing to do with thermal response or are indirect consequences of exposure. My point is that general trends cannot explain all. Some genes are functionally essential for thermal response, while other are not, and I wonder whether looking a GOs and pathways would not be more functionally meaningful. I am not fully convinced that higher number of cold-responsive transcripts directly translates higher cold hardiness.
> 

We completely agree that general trends cannot completely explain anything, but studies that simply report lists of differentially-expressed genes without strong functional work, as are rife in the field, also have substantial explanatory weaknesses. This project focused on broad patterns, testing alternative mechanistic hypotheses (enhanced response, tolerance, genetic assimilation) that did not depend on the response of any specific GO terms or pathways. However, as suggested, we do look at GO terms that are enriched for the global responses - see the section “Gene set enrichment analysis” in the methods, and subsequent discussion of GO terms (metabolism, transcript regulation, etc) in the Results and Discussion. The results of these analyses, which again are prone to storytelling, are used to support the tests of our specific hypotheses. As as example, we find the result that in response to thermal extremes, the warmer-climate species shuts down metabolism-related genes while the cool-climate species increases metabolism and damage repair responses, including  7 heat-shock proteins, to be quite compelling (third paragraph of Discussion, lns 476 - 482). 


> -L388-390: More Hsp transcripts are expressed in A. carolinensis compared to A.picea when exposed to high thermal limit. This is assumed to support tolerance hypothesis. Interpretations depend on whether Hsp expression is viewed as stress reactive or protective effect.
> 

Our main line of evidence for the tolerance hypothesis in A. carolinensis comes from requiring higher temperatures to get a transcriptional response than A. picea. That is, A. carolinensis tolerates greater heat before mounting a transcriptional response. This indicates that there must be other mechanisms that support this tolerance. We note that A. carolinensis activates fewer heat shock proteins than A. picea, and suggest other mechanisms (proteome stability, behavioral quiescence) that may underlie A. carolinensis’ increased thermal tolerance.

> -L413-418 : More transcripts are expressed at lower temperature. This pattern is explained by greater genetic and phenotypic variation at the lower temperature extreme.  Does more genetic and phenotypic variation necessary translates into more genes being responsive?
> 

See our response above - we have removed this section.

> -L492-496:  It is mentioned that A. Picea responds by increasing metabolism and activating heat shock proteins and other protective molecules. I am not sure that transcriptomics can reach such a level of interpretation. RNAseq shows molecular correlates but can certainly not show direct evidence that proteins abundance and activity or more generally pathways and active molecules are actually altered.

We have changed this statement to (lns 593 - 596):

	A. picea responds by increasing expression of transcripts related to metabolism, heat shock proteins and other protective molecules, whereas A. carolinensis decreases expression of transcripts related to metabolism and likely relies on other mechanisms for thermal tolerance. 


> Transcriptomics targets genes expression which is situated upstream of the functional molecules. Moreover, it does not account for downstream regulations, and genes and proteins expression might operate with considerable time lag. The molecular machinery provide ample possibilities for post-transcriptional modifications and regulation, so that a direct relationship between mRNA and protein levels cannot be assumed or expected. In addition, mRNA fluctuations (birth/death) are considered as an important source of stochasticity in gene expression, especially when mRNA abundance is low.

This point is absolutely correct, and we briefly acknowledged it in the line

	Finally, the mapping of changes in gene expression to organismal fitness is far from direct (Feder & Walser 2005), and large differences in patterns of gene expression may have only small effects on organismal fitness. 

We have added the reviewer’s caveats (lns 566 - 569):

	In particular, functional protein levels cannot be expected to be fully linked to mRNA abundance due to post-transcriptional modification, regulation, mRNA fluctuations and protein stability (Feder & Walser 2005).  
