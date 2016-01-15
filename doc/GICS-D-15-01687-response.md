Dear Dr. Tagu,


We appreciate your comments and the feedback from the reviewers. We have revised the manuscript carefully to address the issues that were raised. While we give a detailed response to each reviewer below, we address the main concerns here.


We apologize that the annotation file was not available; this error has been resolved. The annotation file (hf-113-38), transcriptome (hf113-41) and gene expression files (hf-113-42) are now available from the LTER Data Portal: http://dx.doi.org/10.6073/pasta/05ea6464df30efa2f1e2c7439366bf47. We reiterate that all the analyses are available in the supplemental Technical Report and on Github (https://github.com/johnstantongeddes/ApTranscriptome) , so it is possible to recreate the complete analysis, including generating the transcriptome and quantifying gene expression. 


We followed the suggestions of Reviewer #2 to perform a BUSCO analysis of our transcriptome. This suggestion has greatly improved the paper as we can report that the transcriptome captures the majority of genes (~80%) in these species. This number is better than many BUSCO results for transcriptomes of non-model species reported by Simao et al (2015) in supplementary table (http://busco.ezlab.org/files/BUSCO-SOM.pdf), especially given that some genes are not likely to be expressed in the temperature conditions we tested. In addition, the number of duplicated BUSCOs is not excessive, indicating that divergence between the species is not great enough to cause problems with assembly of the joint transcriptome. 


Finally, in a number of places, the reviewers request a more complete description of the number of pathways and genes involved in temperature-specific responses. We have addressed this by including a Supplemental Table (Table S1) with all the responsive transcripts and their annotation. We do not include additional analyses as they would not produce data that would allow us to differentiate among our hypotheses of thermal adaptation, and thus would only make the manuscript more confusing. 


We hope that you now find this manuscript suitable for publication in BMC Genomics.


Sincerely,


John Stanton-Geddes on behalf of all the authors




Reviewer reports:


Associate editor:


Dear authors
Your manuscript has been analyzed by 5 reviewers. This helps me a lot to give you feedback for improving the manuscript. I consider your work as original and the data are consistent. The main issue concerns what I call “annotation of the transcriptomes”, as suggested in different ways by several reviewers. You thus need to add knowledge on the putative annotation of RNA-Seq, without getting too far in functional hypotheses (this is not your main goal, and you don’t have functional demonstration). But these annotations will give an added value to the paper. Concerning reviewer 4-2nde issue, I consider that you need to keep on your tracks of an hypotheses-driven paper, but you might follow reviewer’s comment by leveraging may be the way to present theses 3 hypotheses (that are may be not exclusive, that are bases for structuring the work…). Of course, it is crucial that you take into account all the different issues mentioned in the 5 reports. I hope you will be able to
provide us a new version of this manuscript. "Major revision" will leave you more time to proceed to modifications.


Sincerely


Denis Tagu




Reviewer #1: Stanton-Geddes et al. report a fascinating comparison between related warm and cool-climate adapted ant species, exposing them to increasing temperatures in a shared garden (for six months) paradigm. The abstract as well as the performed statistical analysis allow a good distinction between three hypothesis about adaptation to higher temperatures: enhanced response, genetic assimilation and tolerance. The conclusion stresses also the challenges to adapt tolerance within a short time and the discussion is well balanced.


I have only one but a major point of criticism: Annotation, gene descriptions


--1. Certainly, there are numerous caveats about annotation and even more if only large-scale data are considered. However, for this comparison,
at least the annotation should be available. This is not the case, I could read the technical report which is again well written including all required scripts and the statistical analyses involved are strong and well done. However, the server given to access the actual annotation was not reachable for me.


We apologize for this problem. The annotation file is now available through the LTER Data Portal: http://dx.doi.org/10.6073/pasta/05ea6464df30efa2f1e2c7439366bf47


--2. However, this is a more general point: Only a certain fraction of all genes was analyzed, probably only based on statistical significance though also some GO categories are mentioned.
The analysis on adaptation to temperature using transcriptome data is however incomplete for the reviewer if not at least some clear and transparent part of analysis regarding
categories and pathways involved and affected by the temperature change is done and results on this are given including the original annotation data to back-check the claimed results.
As it stands, one would have to do the reassembly and annotation again as only the really raw files are deposited. This is not good.




The assembly, annotation and gene expression files, as well as the scripts that generated them, are now all available, as noted above. We concur that the information on the annotation for the individual thermally-responsive transcripts is useful information for readers of the manuscript, and now provide this in an additional supplementary Table … In conjunction with the complete annotation file, provided as noted above, any subsequent analyses and checks should be possible.




--3. Furthermore, make crystal clear how you selected the limited number of genes involved, either by stating involved GO categories, or, if you do not like this criterion, then how the threshold used translates into number of genes compared.


As described in the manuscript (lns …) we  “selected the limited number of genes involved” using a statistical-approach that identified the genes that showed a thermally-responsive beyond the null expectation generated by random reassignment of the data. The GO categories that were enriched were determined by performing a gene set enrichment analysis, as described in the Methods, lns 37-45 of page 15. The full code to reproduce the analysis are available Technical Report.




--4. In addition, for me the analysis is incomplete if not at least there is an attempt to compare pathways activated by higher temperature looking at the cold adapted species versus the warm adapted one. For instance, the cold adapted species should express more genes involved in stress response and it would increase confidence in the analysis to find this confirmed and confirmed such that I can verify the conclusion as a reader which again really requires the concrete genes and their annotation not just the raw data.


While we concur that this would be a valid analysis, it is not the approach we choose in this paper. As we state in the paper, in the absence of functional proof of the validity of the annotation, such attempts are prone to “storytelling” and it was our goal to avoid this problem. 




--5. Best would be to underlie the evolutionary hypothesis on the three modes of adaptation by looking at the actual pathways responding to the thermal changes and in addition at selected individual genes whether their function supports the hypothesis.


As the authors should have their annotation files there (even if I can not access the server), they should be able to easily produce about 2 figures including this valuable information about involved pathways and how many genes in both species respond to which temperature change and species-specific differences in numbers and kind (Figure 1, detailed Table in supple material) and pinpoint their specific hypothesis by looking at a few selected concrete genes (Figure 2 or Table).


The existing Table 1 contains the number of genes in both species that respond to each temperature change.


Again this is a nice paper but to add these points is a MUST in this reviewer´s opinion, otherwise one is not able to really reproduce your claims.


While we agree that pathway analysis is a worthy goal for future work, it would not provide any additional data to test our explicit hypotheses for this study. For this reason, we do not feel that it is a ‘must’ at this point. 




Reviewer #2: This manuscript by J. Stanton-Geddes et al. presents a transcriptomic study aimed at understanding the processes by which two ants from the same genus, may adapt to different temperature clines.
For this, they created a reference transcriptome, and used it to study gene responses of A. picea (adapted to colder climates) and A. carolinensis (living in more temperate regions). Under their experimental setup, workers are placed at different temperature conditions before RNA extraction. Instead of performing standard differential expression analysis, the authors use their data to study differential vector responses, which is in my opinion a much appropriate methodology. This methods allows the classification of genes into different responses categories, which allow the author to test plausible hyoptheses for temperature adaptation (enhanced response, tolerance, etc.).
Overall, the paper is very nice to read and I found the experiments well conducted and well presented.
I have no major concern precluding its publication.
Two minor points could be taken into consideration in order to improve the paper.
1/ The use of the term « reactionome » may be considered a little overhyped. There has been many concerns in the litterature about the use and abuse of -omic names in the recent decade. It's usually aimed at describing an exhaustive catalogue of a particular molecular or physiological component. I wonder whether this is really the case here.


We appreciate the reviewer’s comments and their concern about the term “reactionome”. However, there is a growing interest in how organisms respond to temperature variation given climate change projections, so having an umbrella term within which these studies can be grouped will be useful for the research community, similar to the way “reaction norm” has been useful to previous research.


2/ The reference transcriptome is well detailed in the methods and in the supp. material. I was expecting to find metrics of quality using the CEGMA and/or BUSCO analysis. Analysis using these gene sets should ideally retrieve 100% of full length genes if the transcriptome is complete and exhaustive. That could be presented in order to give people more insight into the resource developed.


We performed the BUSCO analysis and retrieved nearly 80% of BUSCOs, indicating that our transcriptome is largely complete. This number is better than many BUSCO results for transcriptomes of non-model species reported by Simao et al (2015) in their supplementary table (http://busco.ezlab.org/files/BUSCO-SOM.pdf), especially given that we only exposed ants to varying temperatures and would not expect some genes to be expressed at all in these conditions. 


In addition, this analysis nicely addresses another issue; that of the impact of assembly the transcriptomes of the two species together. If this were to cause a problem, we would expect there to be an excess of duplicated BUSCO blocks. However, only 8% of BUSCOs are duplicated which is also in line with other transcriptome results in Simao et al. Thus it does not appear that we have incomplete assembly of the two species. 




Reviewer #3: Review of Stanton-Geddes et al. Thermal reactionomes reveal divergent responses to thermal extremes in warm and cool-climate ant species.


In the presented manuscript the authors analyze transcriptomics data to investigate the molecular transformations occurring during adaption to changing environmental conditions. In particular, the authors aim to address the limitations of transcriptomic approaches by implementing a reaction norm approach to the study of gene expression.


The present work is of interest for the ever-growing community interested in gene expression profiles, and even more for climatic changes purposes. This manuscript does make some important contributions to these two fields because it uses a different methodology across a wide range of temperatures to compare how gene expression profiles vary across an environmental gradient.


I like the technical approaches use to demonstrate/refute some important hypothesis and think that this manuscript represents an interesting and important advance. Also, I would like to thank the authors for their transparency as they have made available all scripts/analysis/data to the readers. The manuscript is also extremely well written and a great pleasure to read.


Overall, I think the methods and arguments are strong and I only have minor suggestions for the authors. Most of the potential issues that I did have when I first read the manuscript were already answer by the authors who clearly know the limitations of their study. For example, I believe that more replicates would have been favorable to the study and to strengthen the results. However, the authors did address this potential issue and their answer is satisfactory for me.


Two points that the authors may wish to address a bit more specifically, 1) How did the authors control for age/behavior variation in their study design? 


This point was addressed in two ways. First, all colonies were maintained in the lab for multiple months prior to sampling, and only active workers were used in the heat shock treatment. Second, we pooled 3 ants from each species in each temperature treatment for RNAseq analysis which should limit differences due to age/behaviour variation. 


2) Did the de novo transcriptome assembly performed on both species together affect their results?


Certainly. In our first pass, we assembled each transcriptome separately but found an extremely high degree of overlap. Assembling the transcriptomes together and then mapping reads separately for each species to the combined assembly, containing the shared transcripts as well as those distinct to each species, was the more conservative approach when testing for divergent expression patterns between species, so it was the method we selected. As noted in the methods, we used a number of post-assembly steps (CAP3 and UCLUST) to ensure that the combined transcriptome was as de-duplicated as possible.


Abstract -


Page 2 L 22: The term thermally responsive can be confusing - Can the authors rephrase this sentence?


We have rephrased this “changed expression with temperature.”


Background -


The introduction focuses on explaining the hypothesis tested and the methods. It is well written, however I think more explanations on why ants were chosen as a model organism for this particular topic may be beneficial for the readers. How can climate change affect ant in particular? Do we even know?


There is abundant research on how climate change can impact ants in particular (Pelini et al. 2012; Diamond et al. 2012; Bewick et al. 2014) and ectotherms (Kingsolver et al. 2013; Sunday et al. 2010) in general. Some of this information was contained in the first paragraph of the Methods, but we concur that it may be more relevant earlier. We have expanded on this information in the fourth paragraph of the Background section.  


Since how many years have these two species diverged? How closely related are they? Can the authors speculate on whether there is any consequence or bias introduced (e.g. some gene expression profiles may change for non-climatic reasons but rather phylogenetic)?




A dated analysis of species divergence has not been done and is beyond the scope of this paper, but the species likely diverged recently given the history of glaciation in this region; A. picea occupies area that was glaciated ~16k years ago.  The close phylogenetic relationship between these two species is documented in DeMarco & Cognato (2015), which we have now cited in the manuscript. Neutral drift in gene expression profiles is absolutely a concern, and we address this in the 8th paragraph of the Discussion, beginning, “Moreover, by characterizing response… the reactionome approach can help to distinguish selection from neutral drift in gene expression…”. Briefly, we tested if the number of transcripts that switched expression categories (e.g. from Low to High) differed by more than would be expected by chance.


Page 4, L34: That missing? Rephrase


It’s unclear what the reviewer is referring to here; where is “that” missing? 


It is a personal recommendation, but I do appreciate a couple of sentences at the end of the introduction restating clearly what the aims of the manuscripts are. Your introduction abruptly ends with the descriptions of the final hypothesis, and we are left a bit in the open.


We appreciate this suggestion and have added a final summary paragraph to the Background section: 


“To summarise, in this project we generated the transcriptomes of two closely-related temperate ant species, and quantified their gene expression across a wide range of thermal conditions. We then evaluated three non-mutually exclusive hypotheses (enhanced response, tolerance and genetic assimilation) of the evolution of thermal adaptation by comparing the number and expression patterns of transcripts between species in response to extreme low and extreme high temperatures. Finally, we used gene ontology information to determine which gene products and pathways are involved in thermal adaptation in the two species.”






Results -


Page 5, L4: fit? fitted?


Changed, though this is acceptable American English: http://grammarist.com/usage/fit-fitted/


Page 5, L17: Add the percentage for this value so we can compare with above.


Done.


Page 5, L9: relationshIp


Fixed, thanks.


Page 5, L9; retained? For what? For which analysis?


We clarified this “..as true positive transcripts for all further analyses.”


Page 5, L51: confusing, can you rephrase this sentence?


This sentence is expanded on in the following paragraph regarding changes in the overlap in response patterns between species, so we simply deleted it here.


Page 6, L 27: I agree with the authors that GO terms are usually good arguments for story telling, however in some cases they may be useful, thus I appreciate their efforts to include them in this manuscripts. However, the paragraph is a bit repetitive and could be made more attractive by adding a supplementary figure showing the enriched GO term for each category. Also, please cite the Supplementary tables with the GOterm already in this paragraph (I could only find it in the discussion).


We have added the suggested citation.




Discussion -


L 8, Page 8: about twice


Fixed, thanks.


L 36, Page 8: what do you mean by growing season?


The time period when temperatures are suitable for foraging and colony development.


Page 8, L 57: "the" twice


Fixed, thanks.


Page 9, L 13: This sentence is very long and difficult to understand, rephrase it




What is the life span of these species? Can the authors speculate on whether there are any consequences or bias introduced from the sampling design? The age of the sampled workers may have varied across temperatures and species and somehow affected the expression profiles. How did the authors control for age and behavior (nurse/forager) variations.


Please see our response above.


Page 9, L 58: Can the authors expend their reasoning presented here. How would drift solely affect their results?


Phylogenetic drift alone may have caused divergence in whether transcripts are thermally-responsive in each species, and the critical temperature of each transcript. If this were the case, we would expect no systematic patterns in changes between expression categories (random changes from one expression category to another), and no overall change in the overall critical temperature thresholds between species (some up, some down). In contrast, we see evidence for an excess of changes from some expression categories to others (Fig 1), and changes in critical temperature thresholds in line with predictions (Fig 2).




Methods -


Page 11, L 22: this part is interesting and could be moved to the introduction, as it presents some arguments on why/how ants are a good model organisms


We have addressed this comment above.


Are the species geographical ranges overlapping?


The distributions of A. rudis, the species-complex including A. carolinensis, overlaps with A. picea in Pennsylvania and New York, and along altitudinal gradients further south. However, it’s unknown the extent to which both species are found at a site in a given region. We have added this information to this section.




Page 12, L 23: I am not really getting this point, why was the sampling done over 12 days? One temperature each day?


To control for any differences in circadian rhythm that influence gene expression. We have added a note to this in this section.


I am really surprised that the de novo assembly done for both species together yield such great results, which brings me back to my initial question, how long has it been since these two species diverged? Would this study design influence your ability to detect for example novel genes that are expressed in one species compared to the other and related to environmental adapations?  Did the authors try to assemble the de novo transcriptomes separately and compare them?


We have addressed this comment above. 


Page 15, L 53: repetitive (twice in the same paragraph).








Reviewer #4: In this paper, the authors examine how two ant species that live in different thermal habitats respond to exposure to a broad gradient of temperatures. By doing this, the authors are able to model the reaction norm of each transcript in relation to temperature. They then classify transcripts in each species into different categories, depending on how transcript abundance changes across temperatures. Counts for these different classifications are used to test mechanistic hypotheses about thermal adaptation (i.e., enhanced response hypothesis, tolerance hypothesis, and genetic assimilation hypothesis).


In general, I found the problem addressed by this paper interesting and important. Also, the data were novel, the methods were largely appropriate, and the ant species system seemed well suited to the main problem. I was largely enthusiastic about this study, but had some significant concerns regarding the present version of this manuscript:


We appreciate these comments.


1) The authors note that all fully sequenced ant genomes have around 18,000 genes, yet the analyses in this paper are focused on nearly 100,000 transcripts.  This implies that many genes are represented multiple times among the transcripts. This problem is a major potential confounder to the study, as the authors partition the transcripts into different reaction norm classes. Numbers in these classes are then used to make inferences about broader mechanisms of thermal evolution. The potential for misinference here seems high. I am not sure this concern is mitigated by simply acknowledging the complexities of de novo transcriptome assembly. It seems like more effort and/or a more conservative strategy towards which de novo assembled transcripts are used in the study is important.  Perhaps even BLAST against genes in the existing ant genomes might be useful in grouping together current transcripts into a smaller number of entities.


Transcripts can represent alternative splicing of the same gene, which may have truly different responses to temperature and thus not confounding at all. That said, there are certainly problems for mis-inference due to issues we discuss in the paper. While we could spend many more hours trying alternative strategies for analysis, a better approach is to generate better data using a genome assembly, which is planned. However, at this time we have to do the best analysis possible with the data we have. 




2) I thought the way the three hypotheses were presented up front and imposed on the data seemed a little forced. I agree with the authors that it is worthwhile to get at general mechanisms of thermal evolution. However, the evolution of thermal response may not possible to infer through the strategy the authors propose (i.e., by looking at the distribution of reaction norms across individual genes). Thermal response is likely controlled by many genes that may act combinatorially or may function as components of complicated gene regulatory networks. The approach here seems like an oversimplification. I felt this would be a much stronger paper if these hypotheses were proposed in the conclusion/discussion, instead of serving as a general framing for the paper. The data are very interesting and don't need such a framing to be compelling.


We appreciate this opinion, but based on previous feedback, we felt that it was necessary to present the hypotheses upfront to bring the reader along. It’s absolutely true that the evolution of thermal tolerance may be due to complex interactions or regulatory networks, but in this paper we ask the simple question if patterns of expression response of individual genes are indicative of thermal evolution. The three hypotheses form specific predictions that we can test and reject. Where we reject our hypotheses, it is likely due to more complex factors.


3) For a paper targeting BMC Genomics, I thought there could be more genomics to this paper. E.g., the paper only very briefly describes the GO Enrichment for particular classes of thermal response in each species. I recognize the challenge associated with a non-model species and a de novo transcriptome assembly, but it seems like something should be possible. I was also a little concerned that the observed GO enrichments might be impacted by the issue I describe above in (1).


In our experience in following the literature, transcriptomic studies fall well within the purview of this journal, and quantitative studies are not inherently less “genomic” than functional gene-centered analyses.  The purpose of this paper is to present an alternative framework for analysing transcriptome data to understand thermal adaptation in a non-model species, not to repeat what others have done before. 




Reviewer #5: The authors of this paper examine transcriptional responses across a thermal gradient in cold and warm adapted ant species to test mechanistic hypotheses of thermal adaptation.  One novel contribution of the paper is that they propose a new statistical method for examining transcriptional reaction norms, which I think could prove to be a useful tool for other researchers interested in the role of transcriptional plasticity in adaptation (although admittedly I do not feel qualified to critique some of the statistical methods). Overall the paper is well written, but I do have some suggestions that I think could help improve the paper.


By stating that they "extend the reaction norm approach to RNAseq analysis and introduce the reactionome, which we define as a characterization of the reaction norm for all genes in an organism's transcriptome across an environmental gradient" it makes it sound like the authors are the first to consider this type of approach.  They seem to overlook literature that either suggests or even applies this approach (e.g. Genomic reaction norms: using integrative biology to understand molecular mechanisms of phenotypic plasticity, and several recent papers on Rhagoletis, killfish, Drosophila, and spider mites in Molecular Ecology examine genomic reaction norms and their potential contribution to adaptation).  I agree that the analytical approach the authors use appears to be novel, but that should be distinguished from the general idea of looking at genomic reaction norms.


We appreciate the reminder of Aubin-Horth paper, which we inadvertently left out of our Discussion. We’ve added this citation to the section and highlighted how our study differs from these (paragraph beginning, “A major contribution of this study…”)


“Similar approaches have been used in other species (Aubin-Horth 2009, Hodgins-Davis 2012), but to our knowledge, none have applied a regression approach to identify a complete list of responsive transcript across an environmental gradient. “


I have some questions about the "tolerance hypothesis." First, I don't really understand why "tolerance" describes this scenario but does not also apply to an enhanced reduction response.  What exactly is meant by "tolerance?" Second, I don't really understand why structural changes necessarily imply a "shift to an alternate suite of tolerance genes or pathways." Couldn't you have structural changes to induced genes?  


The enhanced response hypothesis posits that thermally-responsive genes are the same between species, but genes respond at a less extreme temperature. In contrast, the tolerance hypothesis posits that genes becomes less responsive by shifting to ‘non-responsive’ or another expression category. By structural changes, we’re referring to changes to proteins or other structures that do not respond on the time-frame (1 hour) tested in this experiment. 


It seems like the main question the authors are addressing has to do with whether plasticity facilitates subsequent adaptive evolution, in which case you would predict that induction of plastic genes specifically would be enhanced (thus you would have considerable overlap in induction patterns between species) vs. the alternative that other genes or pathways would be involved (and you would not expect overlap).  It is just a suggestion, but maybe the others would want to consider framing the study along these lines, or at least explain their rationale a little
bit more.


This is exactly what the genetic assimilation hypothesis predicts - that changes that are induced to a temperature extreme in one species become constitutively expressed in the other species that experiences those conditions more frequently. We tested this hypothesis and found no evidence for it (Fig. 3). 


How genetically divergent are these species?  If known this would be helpful information, particularly with reference to the transcriptome assembly as from my understanding the transcriptome was built using combined reads from both species.  How likely would it be, for example, that Trinity would assemble alternative alleles from each species into separate components due to genetic divergence? More description and/or justification for this approach may be necessary, especially if the genomes are pretty divergent.


Please see our response to Reviewer #3. 


I think a generic figure that shows the predicted shape of the curves for the different categories (high, low, intermediate, bimodal) would be helpful.


We have added this figure as suggested.


If I understand correctly the best fitting model for each gene was used to predict count levels for the various categories that the authors were interested in.  It seems like the validity of this approach would depend heavily on how well the model actually fits the data.  Is there a way to provide information that would justify this approach?  Also, to clarify the thermal responsive transcripts were found by a model that included terms for species and temperature.  Predicted counts were then determined separately for each species (so any species term drops out of the model)?  This is the way I understand it, but I think it could be clarified a little bit in the text.


This is generally correct, but the species-specific predictions included the species term, so any species-specific effect is included in the model. This is challenging to describe concisely, and we’d refer anyone truly interested in the details to the Technical Report.


On a few occasions the authors mention that they corrected for the false discovery rate, but