<h1>Genomic Selection</h1>

<hr />

Anything to do with genomic selection. 

## Linux Scripts

Anything related to genomic selection that I've written in bash. 

`fmtGenoR.sh` is a script to convert BLUPF90 format to space delimited format for R/Julia/Python/etc. 

<h2>Single-Step GBLUP</h2>

Single-step GBLUP was an addition to the traditional GBLUP (<cite>VanRaden, 2008</cite>) by including non-genotyped animals. This method blends the <b>A</b> and <b>G</b> matrices into one <b>H</b> matrix. The form of H is intense, but the form of H<sup>-1</sup> is very simple. <cite>Legarra et al. (2009)</cite> derived the <b>H</b> matrix at the same time <cite>Christensen and Lund (2010)</cite> did. There is a nice review by <cite>Legarra (2014)</cite> in Livestock Science. 

<h3>APY</h3>

Algorithm for Proven and Young (APY) was created by <cite>Misztal (2014, 2016)</cite>. This algorithm uses the fact that there is limited information in the <b>G</b> matrix. You only need to do the direct inverse of a subset of animals in <b>G</b> and then use recursion to get the rest of the inverse. This has proven to be as accurate (or more so) than the raw inverse. 

## Creating the G matrix with Julia

`createG.jl` is a scipt for Julia (probably out of date at the pace that Julia is developing) 
for creating the <b>G</b> matrix "by hand". 



