#==============================================================================#
# Test the createH Function
#==============================================================================#

# set working directory
  setwd("~/Documents/Programming/R/Animal_Breeding/Relationship_Matrices/")
  
# library
  library(MASS)

#==============================================================================#
# Read in Data
#==============================================================================#
  
# read in data
	animal  <- 13:26
  
# create data frame 
	data.11.1 <- data.frame(animal, 
			sire  = c(0,0,13,15,15,14,14,14,1,14,14,14,14,14),
			dam   = c(0,0,4,2,5,6,9,9,3,8,11,10,7,12),
			mean  = rep(1,length(animal)),
			EDC   = c(558,722,300,73,52,87,64,103,13,125,93,66,75,33),
			fat_DYD = c(9.0,13.4,12.7,15.4,5.9,7.7,10.2,4.8,7.6,8.8,9.8,9.2,11.5,13.3),
			SNP1  = c(2,1,1,0,0,1,0,0,2,0,0,1,0,1),
			SNP2  = c(0,0,1,0,1,1,0,1,0,0,1,0,0,0),
			SNP3  = c(1,0,2,2,1,0,1,1,0,0,1,0,0,1),
			SNP4  = c(1,0,1,1,2,1,1,0,0,1,0,0,1,1),
			SNP5  = c(0,0,1,0,0,0,0,0,0,1,0,1,1,0),
			SNP6  = c(0,2,0,1,0,2,2,1,1,2,1,1,2,2),
			SNP7  = c(0,0,0,0,0,0,0,0,2,0,0,0,0,0),
			SNP8  = c(2,2,2,2,2,2,2,2,2,2,2,2,2,1),
			SNP9  = c(1,1,1,2,1,2,2,2,1,0,2,0,1,0),
			SNP10 = c(2,0,2,1,2,1,0,0,2,0,1,0,0,0))
	rm(list="animal")
  
# A matrix and inverse
	animal <- 1:26
	sire   <- c(rep(0,12), data.11.1$sire)
	dam    <- c(rep(0,12), data.11.1$dam)
	ped    <- data.frame(animal, sire, dam)
	rm(list=c("animal","dam","sire"))
	
# display pedigree
  print(ped)
  
# subset to get M
  M <- data.11.1[6:14, c(1, 7:16)]
  
# write pedigree to file
  # write.table(data.11.1[,1:5], file="data_11_1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
  # write.table(ped, file="pedigree_11_1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
  # write.table(M, file="M_11_1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
#==============================================================================#
# Use create H funcion
#==============================================================================#

# Source the createH() function
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                  "Relationship_Matrices/createH.R"))
  
# use function
  iH <- createH(ped=ped, wts=c(0.95, 0.05, 1, 1))
  round(iH, 2)
  

  
# scramble ped and M matrix to try to get the same answer
  ped2 <- ped[c(15:20, 1:5, 21:26, 6:14), ]
  M2   <- M[c(7:9, 4:6, 1:3), ]
  
# use function
  iH2 <- createH(ped=ped2, M=M2, wts=c(0.95, 0.05, 1, 1))

  
  
  
  
# write out iH
  write.matrix(iH, file="InverseH.txt", sep=" ")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  