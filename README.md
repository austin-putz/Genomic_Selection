# Genomic_Selection

Anything to do with genomic selection. 

## Single-Step GBLUP

Single-step GBLUP was an addition to the traditional GBLUP (VanRaden, 2008) by including non-genotyped animals. This method blends the <b>A</b> and **G** matrices into one **H** matrix. The form of H is intense, but the form of H^{-1} is very simple. Legarra et al. (2009) derived the **H** matrix at the same time Christensen and Lund (2010) did. There is a nice review by Legarra (2014) in Livestock Science. 

### APY

Algorithm for Proven and Young (APY) was created by I. Misztal (2014, 2016). This algorithm uses the fact that there is limited information in the G matrix. You only need to do the direct inverse of a subset of animals in G and then use recursion to get the rest of the inverse. This has proven to be as accurate (or more so) than the raw inverse. 
