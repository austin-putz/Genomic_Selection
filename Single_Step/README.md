# Single Step GBLUP

This page is for things related to Single-step GBLUP.

The `createH.R` and `createH_test.R` are only examples. 
They have not been tested really. But they should help you 
understand the H matrix. It uses VanRaden (2008) basic method
ZZ'/sum2pq. 

You can find the software here: 
[UGA Home Page](http://nce.ads.uga.edu/wiki/doku.php)
or 
[Downloads Page](http://nce.ads.uga.edu/wiki/doku.php?id=distribution)

## How to use BLUPF90 software

1. You have to first run the renumf90 program [here](http://nce.ads.uga.edu/wiki/doku.php?id=readme.renumf90)
  * run with `renumf90 <<< file.par`
2. Then run the application program [here](http://nce.ads.uga.edu/wiki/doku.php?id=application_programs)
  * run with `blupf90 <<< renf90.par` (or choose your needed program from the list)
  * `renf90.par` is the general parameter file written by `renumf90` 

## APY

Algorithm for Proven and Young (APY) was created by I. Misztal (2014, 2016). This algorithm uses the fact that there is limited information in the **G** matrix. You only need to do the direct inverse of a subset of animals in **G** and then use recursion to get the rest of the inverse. This has proven to be as accurate (or more so) than the raw inverse. It started by using proven (high accuracy animals), but now has been shown to work well with enough animals using a random sample of the population. His 2014 paper introduces it and the 2016 paper explains more about why it works. 
