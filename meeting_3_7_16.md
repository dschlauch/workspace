Biweekly Update
========================================================
author: Dan Schlauch
date: March 7, 2016

Agenda
========================================================

- Major Objectives for Manuscript
    - Better communicate the story
    - Biological insight
        - Send to Ed, Craig when ready for their input.
- Address list of concerns with last draft



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
    

New supporting results:
========================================================
The following plots demonstrate that commonly used NI methods to not produce reproducible results.

Edgeweight differences for COPD and Smoker Controls were taken across studies.

Edgeweight differences are essentially uncorrelated across studies.

This holds for all tested NI methods.

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

Network differences do not validate across studies!
========================================================
Network differences between cases and controls are **completely uncorrelated** between COPD studies using commonly the used WGCNA, CLR, and ARACNE. (Far less than 1% of the variation in one study is explained by the variation in any other)

### So, are NI methods for context-specific networks useless?

Maybe not.  Is there any way to observe drivers of state transitions at the network level using this information?

Yes! Our [NAME] method focuses on regulatory patterns rather than individual edges and DOES validate across studies. This is our "story".

MONSTER validation
========================================================
<embed src="TM_manuscript/figures/Magnitude_comparison_same_tissue.pdf" width="450" height="400" type='application/pdf'>
<embed src="TM_manuscript/figures/Magnitude_comparison_diff_tissue.pdf" width="450" height="400" type='application/pdf'>

Differential expression methods cannot find these TFs
========================================================
<embed src="TM_manuscript/figures/dTFIvsLIMMA.pdf" width="800" height="500" type='application/pdf'>


Misc. manuscript tasks...
========================================================
JQ Comments: Figure for transition heatmap <br>
<embed src="TM_manuscript/figures/TM_heatmap.pdf" width="450" height="400" type='application/pdf'>


Misc. manuscript tasks...
========================================================
<img src="TM_manuscript/figures/Venn.png" width=500 height=500> 

Misc. manuscript tasks...
========================================================
Name suggestion:

**MO**deling **N**etwork **S**tate **T**ransitions from **E**xpression and **R**egulatory data

(_MONSTER_)

Misc. manuscript tasks...
========================================================
Figure 1.

JQ comments:
"the matrices should represent the networks" <br>
<img src="TM_manuscript/figures/figure1a.png" width=300 height=300> vs
<img src="TM_manuscript/figures/figure1.png" width=300 height=300> 


Misc. manuscript tasks...
========================================================
Equations here, not rendering

Misc. manuscript tasks...
========================================================
<img src="TM_manuscript/figures/Venn.png" width=500 height=500> <br>
Do we want/need this?  Can table replace it? Or keep table in supplement (it's big).


Misc. manuscript tasks...
========================================================
More supporting literature validation

*Autophagy in chronic obstructive pulmonary disease: Homeostatic or pathogenic mechanism?*
Ryter (2008 and 2010), others
- E2F4 (one of the top hits in all studies) binds with EGR-1 (top 5 hit in ECLIPSE and LTCOPD, also towards the top of LGRC and COPDGENE) in response to stress (such as smoke).

*Egr-1 Regulates Autophagy in Cigarette Smoke-Induced Chronic Obstructive Pulmonary Disease*
Chen et al. (2008)

Misc. manuscript tasks...
========================================================
ELK1, ELK4 implicated in asthma.

*Emerging role of MAP kinase pathways as therapeutic targets in COPD*
Mercer (2006)
ELK1 phosphorylated by ERK1/2
