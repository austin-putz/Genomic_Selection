# Single Step GBLUP

This page is for things related to Single-step GBLUP.

The `createH.R` and `createH_test.R` are only examples. 
They have not been tested really. But they should help you 
understand the H matrix. It uses VanRaden (2008) basic method
ZZ'/sum2pq. 

You can find the software here: 
[UGA Home Page](http://nce.ads.uga.edu/wiki/doku.php){:target="_tab"}
or 
<a href="http://nce.ads.uga.edu/wiki/doku.php?id=distribution" target="_blank">New Tab</a>

## How to use BLUPF90 software

1. You have to first run the renumf90 program [here](http://nce.ads.uga.edu/wiki/doku.php?id=readme.renumf90)
  * run with `renumf90 <<< file.par`
2. Then run the application program [here](http://nce.ads.uga.edu/wiki/doku.php?id=application_programs)
  * run with `blupf90 <<< renf90.par` (or choose your needed program from the list)
  * `renf90.par` is the general parameter file written by `renumf90` 

## APY

Algorithm for Proven and Young (APY) was created by I. Misztal (2014, 2016). This algorithm uses the fact that there is limited information in the G matrix. You only need to do the direct inverse of a subset of animals in G and then use recursion to get the rest of the inverse. This has proven to be as accurate (or more so) than the raw inverse. 
