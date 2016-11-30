---
title: "Comparison of FROGS/MOTHUR/UPARSE on simulated and synthetic microbial communities"
author: "Mahendra Mariadassou, Géraldine Pascal"
date: "17 juin 2016"
output:
  html_document:
    number_sections: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: true
      smooth_scroll: true
  pdf_document:
    toc: yes
---









# Introduction

This is an R markdown document intended to compare the performances of FROGS, MOTHUR and UPARSE in terms of accuracy on both simulated and synthetic microbial communities. 

## Metrics used for comparison 

The results of FROGS, MOTHUR and UPARSE are compared using three different metrics: 
* **Divergence**: Bray-Curtis distance (expressed in percent) between the true taxonomic composition of the community and the one inferred by the otu-picking tool. The divergence is measured at all taxonomic ranks from Phylum to either Genus (utax) or Species (Silva). 
* **FN**: Number of false negative taxa (*i.e.* present in the original bacterial community but not discovered by the otu picking method); 
* **FP**: Number of false positive taxa (*i.e.* discovered by the otu picking method but not present in the original bacterial community)

## Experimental design

The experimental design differed for the simulated communities (for which a full-factorial design was used) and the synthetic communities. 

### Simulated bacterial communities

The simulated communities were built according to the following design:
* **Databank**: Biobank from which taxa were drawn to construct theoretical communities, either Silva (*silva*) or Utax (*utax*). 
* **Number of OTUs**: 20, 100, 200, 500 and 1000;
* **Abundance distribution**: abundances of OTUs were either uniform (*uniform*) or sampled from a power distribution (*power_law*);
* **Dataset**: Theoretical communities. Each dataset (5 for each combinaison of abundance distribution and number of OTUs) correspond to a unique **ideal** bacterial community specified by its own taxa set and corresponding vector of relative abundances. 
* **Set Number**: Biological replicates (10 for each dataset), *i.e.* communities created by sampling organisms with replacement in the theoretical communities.  
* **Amplicon**: variable region of the 16S rRNA used to produce the ampicon sequences, either the V3-V4 (*V3V4*) variable region or the V4-V4 (*V4V4*) variable region

This resulted in a total of 2 databanks $\times$ 5 community sizes $\times$ 2 abundance distribution $\times$ 10 theoretical communities $\times$ 10 replicates for each theoretical community $\times$ 2 amplicons $=$ 2000 samples (1000 per databank). 



### Synthetic communities

The experimental design used for the synthetic community was sligthly different:

|nb_OTU |amplicon |abundance_law | count|
|:------|:--------|:-------------|-----:|
|20sp   |V3V4     |even          |     4|
|20sp   |V3V4     |staggered     |     4|
|20sp   |V4V4     |even          |     1|
|20sp   |V4V5     |even          |     1|
|4sp    |V3V4     |uneven        |    19|

Three samples corresponding to communities of size 20 with abundance distribution *even* were used to compare the amplicon *V3V4*, *V4V4* and *V4V5* (1 per amplicon). 8 samples corresponding to communities of size 20 and amplicon *V3V4* were used to compared abundance distribution *even* and *staggered* (4 per distribution) and finally, 19 samples (community size 4, amplicon *V3V4* and distribution *even*) were used to compare the accuracy of the different otu picking methods. 


# Analysis of the results

## Statistical Analysis: Material and Methods

For each of the three metrics (divergence, FN and FP) we performed two-sided paired test, either parametric (paired t-test) or non-parametric (signed rank test, also known as paired mann-whitney test) to assess the difference in accuracy between FROGS and each of the competitors. 

The tests were peformed at the theoretical community levels (*dataset*) using biological replicates (*set_number*) as replicates. We chose to compare the methods at this level because it the finest one for which we have replication. Pooling different theoretical communities and/or abundance distributions to compare the method at higher levels (*e.g* community size $\times$ amplicon) will blur the signal as a method may be  outclass the others for even abundances but perform worse on different abundance disrtibutions. 

For each theoretical community, we declared FROGS better (resp. worse) than its competitor when the test was significant at the 0.05 level and FROGS had a lower (resp. higher) metric than its competitor. When the test was not significant, the methods were declared tied. Finally, we aggregated the results to count for each condition (community size $\times$ abundance distribution $\times$ amplicon) the number of theoretical communities favoring one or none of the methods. 

Before presenting the statistical analysis per se, we first present the results graphically for each of the databank and for the synthetic communities. 

## Silva databank

[//]: # (Import and format data)









### Vizualisation 

#### Divergence 

The comparisons of divergence at the sample level in the scatterplots shows that on average, FROGS has comparable but better performances than MOTHUR and UPARSE: most samples end up in the upper left corner (corresponding to the region "divergence FROGS < divergence competitor") but no too far away from the first diagonal (grey line). 


![plot of chunk plot-data-raw](Figures_20160722/Fig-english/plot-data-raw-1.png)

A more traditional representation using boxplot of the excess divergence of FROGS, with samples from all theoretical communities pooled together, confirms the results: FROGS has similar (compared to UPARSE) or lower (compared to MOTHUR) divergence for the vast majority of samples. Note that the y-range was reduced from $[-40, 3]$ to $[-15, 3]$ in order to exclude outliers (communities with low FROGS but very high MOTHUR divergence) and zoom in on the boxplots. 

![plot of chunk plot-data-raw-boxplot](Figures_20160722/Fig-english/plot-data-raw-boxplot-1.png)

Finally a focus on the accuracy of FROGS alone shows that divergence levels vary between 0 and 10% and as expected, is higher for fine classification (Species) than for coarse ones (Phylum). Unsurprisingly, the V3V4 amplicon gives less distorded view of communities than the V4V4. 
![plot of chunk plot-frogs-only-divergence](Figures_20160722/Fig-english/plot-frogs-only-divergence-1.png)

#### False Positive and False Negative OTUs

We repeat the graphical exploration of the resutls with False Positive and False Negative OTUs. A first representation shows that use of the V3V4 amplicon leads to more false postive and less false negative than the V4V4. The graphics also highlight the gigantic number of false positive inferred by mothur (up to 20 times more than the real community size).  

![plot of chunk plot-fp-fn-boxplot](Figures_20160722/Fig-english/plot-fp-fn-boxplot-1.png)

A focus on FROGS and UPARSE leads to similar patterns: FROGS always produces less false negatives than UPARSE but produces a bit more false positive under power law abundance distribution and a bit less under uniform abundance distribution. 

The lower number of false positive in under power law abundances could be due to the abundance based filters used in UPARSE. 

![plot of chunk plot-fp-fn-boxplot-no-mothur](Figures_20160722/Fig-english/plot-fp-fn-boxplot-no-mothur-1.png)

### Statistical Analysis

#### Divergence

We present the results of the paired tests, either parametric (t-test, top) or non-parametric (signed rank test, bottom). Both tests show that FROGS perform as well or better as UPARSE and MOTHUR in most conditions. The only condition in which FROGS does worse than UPARSE is small community size (20). 

The real strength of FROGS lies in its ability ot give a more accurate view of large communities (size > 200) at fine scales (Species or Genus level). 



![Comparison of FROGS to competing methods using paired t-test.](Figures_20160722/Fig-english/batch-t-test-plot-1.png)



![Comparison of FROGS to competing methods using signed rank tests.](Figures_20160722/Fig-english/batch-wilcox-test-plot-1.png)

#### False Positive and False Negative OTUs

The same paired test as in the previous section reveal that FROGS strictly outperforms MOTHUR in terms of both FP and FN taxas. It also produces less FN than UPARSE. Additionnally, it produces less FP than UPARSE for uniform distributions and more power law ones. Overall, FROGS produces less FP and less FN than either of UPARSE and MOTHUR for high community sizes (>200 for uniform distributions, >1000 for power law distributions). 



![plot of chunk batch-t-test-otus-plot](Figures_20160722/Fig-english/batch-t-test-otus-plot-1.png)



![plot of chunk batch-wilcox-test-plot-otus](Figures_20160722/Fig-english/batch-wilcox-test-plot-otus-1.png)


## Utax databank









### Vizualisation 

#### Divergence 

The comparisons of divergence at the sample level in the scatterplots shows that on average, FROGS has comparable but better performances than MOTHUR and UPARSE: most samples end up in the upper left corner (corresponding to the region "divergence FROGS < divergence competitor") but no too far away from the first diagonal (grey line). 


![plot of chunk utax-plot-data-raw](Figures_20160722/Fig-english/utax-plot-data-raw-1.png)

A more traditional representation using boxplot of the excess divergence of FROGS, with samples from all theoretical communities pooled together, confirms the results: FROGS has similar (compared to UPARSE) or lower (compared to MOTHUR) divergence for the vast majority of samples. Note that the y-range was reduced from $[-40, 3]$ to $[-15, 3]$ in order to exclude outliers (communities with low FROGS but very high MOTHUR divergence) and zoom in on the boxplots. As expected, all methods perform quite similarly up to the order level and the main differences appear at the *Family* and *Genus* levels, where MOTHUR and MOTHUR_SOP produces much larger divergences than competing methods.  

![plot of chunk utax-plot-data-raw-boxplot](Figures_20160722/Fig-english/utax-plot-data-raw-boxplot-1.png)

Finally a focus on the accuracy of FROGS alone shows that divergence levels vary between 0 and 10% and as expected, is higher for fine classification (Genus) than for coarse ones (Phylum). Unsurprisingly, the V3V4 amplicon gives less distorded view of communities than the V4V4. Overall, FROGS recover community compositions very well expect at the genus level for complex communities (size > 200), with uniform abundances and sequenced using the V4V4 region. 
![plot of chunk utax-plot-frogs-only-divergence](Figures_20160722/Fig-english/utax-plot-frogs-only-divergence-1.png)

#### False Positive and False Negative OTUs

We repeat the graphical exploration of the resutls with False Positive and False Negative OTUs. A first representation shows that use of the V3V4 amplicon leads to more false postive and less false negative than the V4V4. The graphics also highlight the gigantic number of false positive inferred by mothur (up to 20 times more than the real community size).  

![plot of chunk utax-plot-fp-fn-boxplot](Figures_20160722/Fig-english/utax-plot-fp-fn-boxplot-1.png)

A focus on FROGS and UPARSE leads to similar patterns: FROGS always produces less false negatives than UPARSE but produces a bit more false positive under power law abundance distribution and a bit less under uniform abundance distribution. 

The lower number of false positive in under power law abundances could be due to the abundance based filters used in UPARSE. 

![plot of chunk utax-plot-fp-fn-boxplot-no-mothur](Figures_20160722/Fig-english/utax-plot-fp-fn-boxplot-no-mothur-1.png)

### Statistical Analysis

#### Divergence

We present the results of the paired tests, either parametric (t-test, top) or non-parametric (signed rank test, bottom). Both tests show that FROGS perform as well or better as UPARSE and MOTHUR in most conditions. The only condition in which FROGS does worse than UPARSE is small community size (20). 

The real strength of FROGS lies in its ability ot give a more accurate view of large communities (size > 200) at fine scales (Species or Genus level). 



![Comparison of FROGS to competing methods using paired t-test.](Figures_20160722/Fig-english/utax-batch-t-test-plot-1.png)



![Comparison of FROGS to competing methods using signed rank tests.](Figures_20160722/Fig-english/utax-batch-wilcox-test-plot-1.png)

#### False Positive and False Negative OTUs

The same paired test as in the previous section reveal that FROGS strictly outperforms MOTHUR in terms of both FP and FN taxas. It also produces less FN than UPARSE. Additionnally, it produces less FP than UPARSE for uniform distributions and more power law ones. Overall, FROGS produces less FP and less FN than either of UPARSE and MOTHUR for high community sizes (>200 for uniform distributions, >1000 for power law distributions). 



![plot of chunk utax-batch-t-test-otus-plot](Figures_20160722/Fig-english/utax-batch-t-test-otus-plot-1.png)



![plot of chunk utax-batch-wilcox-test-plot-otus](Figures_20160722/Fig-english/utax-batch-wilcox-test-plot-otus-1.png)


## Synthetic Communities 

Due to the different design used for synthetic communities, we going to perform focused comparisons of the samples. 







### Amplicon Effect 

We study the amplicon effect on the even community with 20 species as it is the only one for which different amplicons were used. We represent the trends observed at different taxonomic levels. UPARSE seems to do better than FROGS and MOTHUR but there is only one sample per amplicon so we can't assess the significance of that trend. 

Note that all methods have high divergences compared to simulated datasets. This may reflect experimental limitations (sequencing and amplification bias, copy number variations, etc) rather than intrinsic complexity of the synthetic community and/or differences between the methods. 

![plot of chunk amplicon-effect-paired-plot](Figures_20160722/Fig-english/amplicon-effect-paired-plot-1.png)

### Abundance Distribution Effect

We study the abundance distribution effect on the community with 20 species sequenced with the V3V4 amplicon region. FROGS is better than MOTHUR and UPARSE on the staggered community and worse on the even one... 

![plot of chunk distribution-effect-paired-plot](Figures_20160722/Fig-english/distribution-effect-paired-plot-1.png)

Although the base divergence is not very satisfying on either distribution. 
![plot of chunk distribution-effect-boxplot](Figures_20160722/Fig-english/distribution-effect-boxplot-1.png)

As we have 4 replicates for this community, we can compare all methods using a paired t-test (there are not enough replicates for the non parametric signed rank test to reach significance). The statistical analysis confirm that FROGS outperforms MOTHUR and UPARSE on communities with staggered abundances. 



![plot of chunk real-distribution-t-test-plot](Figures_20160722/Fig-english/real-distribution-t-test-plot-1.png)

### Comparison on a Toy Community

Finally, we compare the three methods on a toy community with 20 species and even abundances. We have  19 replicates for this community, enough to use either the signed rank test and the paired t-test. Once again FROGS base divergence level is not fantastic
![plot of chunk real-method-effect-4-species-boxplot](Figures_20160722/Fig-english/real-method-effect-4-species-boxplot-1.png)

but in line with divergences obtained by competitors
![plot of chunk real-method-effect-4-species-paired-plot](Figures_20160722/Fig-english/real-method-effect-4-species-paired-plot-1.png)

The statistical analyses confirms the graphical diagnostic of the boxplot: FROGS is generally better than MOTHUR and tied with UPARSE, except at the Genus rank where it outperforms both of them. It also performs worse than UPARSE at the Phylum rank (for the paired t-test)


```
## Source: local data frame [6 x 9]
## Groups: databank, nb_OTU, dataset, amplicon, abundance_law, rank [3]
## 
##   databank nb_OTU   dataset amplicon abundance_law   rank    measure
##     (fctr)  (dbl)    (fctr)   (fctr)        (fctr) (fctr)      (dbl)
## 1     utax      4 dataset_1     V3V4        uneven  Genus -0.3847674
## 2     utax      4 dataset_1     V3V4        uneven  Genus -1.0759372
## 3     utax      4 dataset_1     V3V4        uneven Family  0.1877899
## 4     utax      4 dataset_1     V3V4        uneven Family -0.9641701
## 5     utax      4 dataset_1     V3V4        uneven  Order  0.2005283
## 6     utax      4 dataset_1     V3V4        uneven  Order -0.7309186
## Variables not shown: pval (dbl), best (fctr)
```

![plot of chunk real-method-t-test-4-species-plot](Figures_20160722/Fig-english/real-method-t-test-4-species-plot-1.png)





![plot of chunk real-method-wilcox-test-4-species-plot](Figures_20160722/Fig-english/real-method-wilcox-test-4-species-plot-1.png)
