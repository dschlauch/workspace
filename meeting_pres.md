Biweekly Update
========================================================
author: Dan Schlauch
date: March 7, 2016

Agenda
========================================================

- Objectives for Manuscript
    - Better communicate the story
    - ??
    - Biological insight
        - Send to Ed, Craig when ready for their input.



The story (Part 1: Setup)
========================================================
- There is a challenge in Network Inference
    - NI methods can be effective, but are weak at context specific results.
    - Existing NI methods can often perform well when measured against gold standards such as ChIP, but context specific validation is lacking.  Differential regulation.
    - Inferring GRNs from expression data and then computing network differences do not yield repeatable results (see later slides)

The story (Part 2: Solution)
========================================================
-  Our method tackles this problem
    1. Our method produces **repeatable results**
        - Strong correlation across independent studies of different tissues.
        - Strongest correlation across same tissues in independent studies.
    2. Our method's results are **biologically relevent** and **independently verified**.
        - Mitochondrial function, other biopoetry needed here.
        - CONDOR
    3. These results **could not have been obtained** via conventional methods.
        - Differential expression, comparing edges...
        
Additionally,
========================================================

- For supplemental arguments...
    1. Our method has a strong intuitive theoretical basis.
    2. Our method works in simulations with known drivers of state transitions.
    3. Our method is fast for high throughput scale data.
    

Edgeweight differences do not validate across studies!
========================================================
<div align="center">
WGCNA <br>
<img src="plots/ECLIPSE_COPDGENE_WGCNA_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LTCOPD_WGCNA_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LGRC_WGCNA_edgeweight_difference_comparison.png" width=200 height=200> <br>
<img src="plots/COPDGENE_LGRC_WGCNA_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/COPDGENE_LTCOPD_WGCNA_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/LGRC_LTCOPD_WGCNA_edgeweight_difference_comparison.png" width=200 height=200> 
</div>

Edgeweight differences do not validate across studies!
========================================================
<div align="center">
Context Likelihood of Relatedness (CLR) <br>
<img src="plots/ECLIPSE_COPDGENE_CLR_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LTCOPD_CLR_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LGRC_CLR_edgeweight_difference_comparison.png" width=200 height=200> <br>
<img src="plots/COPDGENE_LGRC_CLR_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/COPDGENE_LTCOPD_CLR_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/LGRC_LTCOPD_CLR_edgeweight_difference_comparison.png" width=200 height=200> 
</div>

Edgeweight differences do not validate across studies!
========================================================
<div align="center">
Algorithm for the Reconstruction of Gene Regulatory Networks (ARACNE)<br>
<img src="plots/ECLIPSE_COPDGENE_ARACNE_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LTCOPD_ARACNE_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LGRC_ARACNE_edgeweight_difference_comparison.png" width=200 height=200> <br>
<img src="plots/COPDGENE_LGRC_ARACNE_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/COPDGENE_LTCOPD_ARACNE_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/LGRC_LTCOPD_ARACNE_edgeweight_difference_comparison.png" width=200 height=200> 
</div>

Edgeweight differences do not validate across studies!
========================================================
<div align="center">
BERE<br>
<img src="plots/ECLIPSE_COPDGENE_bere_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LTCOPD_bere_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/ECLIPSE_LGRC_bere_edgeweight_difference_comparison.png" width=200 height=200> <br>
<img src="plots/COPDGENE_LGRC_bere_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/COPDGENE_LTCOPD_bere_edgeweight_difference_comparison.png" width=200 height=200>
<img src="plots/LGRC_LTCOPD_bere_edgeweight_difference_comparison.png" width=200 height=200> 
</div>


========================================================


```r
summary(cars)
```

```
     speed           dist       
 Min.   : 4.0   Min.   :  2.00  
 1st Qu.:12.0   1st Qu.: 26.00  
 Median :15.0   Median : 36.00  
 Mean   :15.4   Mean   : 42.98  
 3rd Qu.:19.0   3rd Qu.: 56.00  
 Max.   :25.0   Max.   :120.00  
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](meeting_pres-figure/unnamed-chunk-2-1.png)