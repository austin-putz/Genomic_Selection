# Linux Scripts for Genomic Selection

I'll put some of my scripts here for anyone who wants them. 

## fmtGenoR.sh

BLUPF90 needs the formatted genotypes to start at the same 
character and are usually in one string. 
However, when you read into R, Python, or Julia, you need them
to be separated to more easily read them in. 

Will change the format from:

  Sire001     1020201020 
  
  Sire002     1010202010

to:

  Sire001     1 0 2 0 2 0 1 0 2 0 
  
  Sire002     1 0 1 0 2 0 2 0 1 0
  
This makes it possible to read it into your favorite language easier (R, Python, Julia, etc). 





