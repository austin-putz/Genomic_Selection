#!/bin/bash

#================================================================================#

# fmtGenoR
# This script will take genotypes for BLUPF90 and convert to be
# read into R or another program

#================================================================================#

# CURRENT FORMAT
# ID1     20120120
# ID2     20120201
# ID3     20102021

# DESIRED FORMAT (easily read into R with fread from data.table)
# ID1 2 0 1 2 0 1 2 0
# ID2 2 0 1 2 0 2 0 1
# ID3 2 0 1 0 2 0 2 1

#================================================================================#

# In R
# library(data.table)
# genotypes <- fread("my_file.snp", header=F, sep=" ")

#================================================================================#

# Set file you want to convert
  file=$1

# separate them into 2 separate files
  awk ' { print $1 } ' $file > animals
  awk ' { print $2 } ' $file > genotypes

# now separate the genotypes with a space
  sed 's/\(.\{1\}\)/\1 /g' genotypes > genotypes_sep

# paste them back together (animal + ID)
  paste animals genotypes_sep > $file.reformatted

# remove left over files
  rm -f animals genotypes genotypes_sep


