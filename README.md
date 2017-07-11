# branchpointer

##Prediction of intronic splicing branchpoints

###Introduction

The majority of human genes are spliced, forming a mature mRNA following intron removal and subsequent exon ligation in the spliceosome complex [1]. During this reaction, U1 snRNP and SF1 bind the 5’ splice site (5’ SS) and the branchpoint respectively, and a trans-esterification reaction forms an intron lariat intermediate. A subsequent trans-esterification reaction between the 5’ SS and the 3’ SS removes the intron lariat and forms a spliced RNA product from the flanking exons. Recognition of the sequence-based splicing elements – the 3’ SS, 5’ SS and the branchpoint – by small ribonucleoprotein particles (snRNPs) is a critical step in defining exon boundaries, and therefore the mature mRNA product formed [2].
Branchpoint elements have traditionally been difficult to identify, and their role in splicing regulation and disease poorly understood. Branchpointer is an R package for predicting splicing branchpoints using primary genome sequence and exon annotations alone. Branchpointer uses a machine-learning model, trained with empirical branchpoint annotations, to identify branchpoint elements in introns with best-in-class sensitivity (60.7%) and specificity (97.8%). In addition, branchpointer can evaluate the impact of mutations on branchpoint architecture to inform functional interpretation of genetic variants.

1. Will CL, Lührmann R. Spliceosome structure and function. Cold Spring Harb. Perspect. Biol. 2011;3. 
2. De Conti L, Baralle M, Buratti E. Exon and intron definition in pre-mRNA splicing. Wiley Interdiscip. Rev. RNA. 2013. p. 49–60. 

###Installation

R-package branchpointer can be installed:

    library(devtools)
    install_github("betsig/branchpointer")

After installation, the package can be loaded into R.

    library(branchpointer)

For details of how to use this package, please see the vignette.

###Dependencies

branchpointer requires [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) to be installed if using a local .fa file for sequence retrieval.
Alternatively, an R BSgenome object can now be used (see vignette).

**Package**: branchpointer

**Type**: Package

**Title**: Prediction of intronic splicing branchpoints

**Version**: 1.3.0

**Date**: 2017-07-10

**Author**: Beth Signal

**Maintainer**: Beth Signal <b.signal@garvan.org.au>

**Description**: Predicts branchpoint probability for sites in intronic branchpoint windows. 
queries can be supplied as intronic regions; or to evaluate the effects of mutations, SNPs.
