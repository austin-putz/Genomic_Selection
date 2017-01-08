#==============================================================================#
# APY Algorthim
#==============================================================================#

# source A function
	source(paste0("/Users/austinputz/Documents/Programming/R/Animal_Breeding/", 
                    "Gota_Morota/Pedigrees/createA.R"))

# library(data.table)
  library(pedigreemm)

#------------------------------------------------------------------------------#
# Get data from Mrode's Chapter 11
#------------------------------------------------------------------------------#

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
	
# create A matrix
	A      <- createA(ped)
	A.inv  <- solve(A)
  
#------------------------------------------------------------------------------#
# Create the G matrix
#------------------------------------------------------------------------------#
	
# get SNP matrix
	M     <- data.11.1[, 7:16]
	M     <- as.matrix(M)
  
# get average allele frequencies
	all.freq       <- apply(M, 2, mean) / 2
	P.allele.freqs <- matrix(rep(all.freq, nrow(M)), byrow=T, ncol=ncol(M))
	P              <- P.allele.freqs * 2
	rm(P.allele.freqs)
  
# get Z matrix
# only use first three SNPs and first 8 as reference population
  Z <- M - P
  
# create G
  G <- (Z %*% t(Z)) / sum(2*all.freq*(1-all.freq))
	print(round(G, 2))

# A sub of only genotyped
	A22 <- A[13:26, 13:26]
	
# set alpha and beta
	alpha <- 0.05
	beta  <- (1-alpha)
	
# get G to be invertable
	Gstar  <- (alpha*A22) + (beta*G)
	Gstar2 <- G + diag(ncol(G))*0.0001
	
# invert G
	Ginv  <- solve(Gstar)
	Ginv2 <- solve(Gstar2)
  print(round(Ginv, 2))
  
#------------------------------------------------------------------------------#
# Create the Ginv matrix with APY
#------------------------------------------------------------------------------#
  
# row and col names
  rownames(G) <- 1:ncol(G)
  colnames(G) <- 1:ncol(G)

# split G
	Gcc <- G[1:4, 1:4]
	Gcn <- G[1:4, 5:14]
	Gnc <- G[5:14, 1:4]
	Gnn <- G[5:14, 5:14]
  
# set n1 (number core) and n2 (number non-core)
  n1 <- ncol(Gcc)
  n2 <- ncol(Gnn)
	
# invert core G
	Gccinv <- solve(Gcc)

# solve top part of matrix
	top    <- (-1*Gccinv) %*% Gcn
	bottom <- diag(n2)
	Right <- rbind(top, bottom)

# initiate M
	M <- matrix(0, ncol=n2, nrow=n2)

# set row and column names for G
	rownames(M) <- rownames(Gnn)
	colnames(M) <- colnames(Gnn)
	
  # get Gii matrix
	Gii <- diag(diag(Gnn))

# initialize vector
	Mvec <- matrix(0, ncol=n2, nrow=n2)
	
# loop through to get elements of M
	for (i in 1:n2) {
	  
	  # subset column of Gcn
	  vec = matrix(Gcn[, i], ncol=1)
	  print(vec)
	  
	  # do matrix multiplication
	  print(t(vec) %*% Gccinv %*% vec)
	  Mvec[i, i] = (t(vec) %*% Gccinv %*% vec)
	  
	}

# add to Gii
	M = Gii + Mvec

# pad Gccinv with 0's
	topright   <- matrix(0, ncol=n2, nrow=n1)
	bottomleft <- matrix(0, ncol=n1, nrow=n2)
	bottomright <- matrix(0, ncol=n2, nrow=n2)

# combine
	top <- cbind(Gccinv, topright)
	bot <- cbind(bottomleft, bottomright)
  Left <- rbind(top, bot)

# calculate APY
  GinvAPY <- Left + Right %*% solve(M) %*% t(Right)
  round(GinvAPY, 1)
  
#------------------------------------------------------------------------------#
# Create the Ginv matrix with APY
#------------------------------------------------------------------------------#
  
# compare to Ginv
  cor(c(Ginv), c(GinvAPY))
  plot(c(Ginv), c(GinvAPY))
  
  cor(c(Ginv2), c(GinvAPY))
  plot(c(Ginv2), c(GinvAPY))
  
#------------------------------------------------------------------------------#
# Use the function I created
#------------------------------------------------------------------------------#
  
# source
  source("~/Documents/Programming/R/Functions/APY_function.R")
  
# APY function
  GinvAPY2 <- APY(G, core=c(rep(1,4), rep(0,10)))
  
  cor(c(GinvAPY), c(GinvAPY2))
  
#------------------------------------------------------------------------------#
# Plot out the correlations
#------------------------------------------------------------------------------#
  
  library(ggplot2)
  
# create dataframe
  Gdata <- data.frame(A22=c(A22), G=c(G), 
                        Gstar=c(Gstar), Gstar2=c(Gstar2), 
                        Ginv=c(Ginv), Ginv2=c(Ginv2),
                        GinvAPY=c(GinvAPY2))
  
# plot between G and Gstar (combine with pedigree)
  ggplot(Gdata, aes(x=c(A22), y=c(G))) +
    geom_point() +
    ggtitle("Correlation between A and G") +
    xlab("A") +
    ylab("G") +
    theme_bw()
  
# plot between G and Gstar (combine with pedigree)
  ggplot(Gdata, aes(x=c(G), y=c(Gstar))) +
    geom_point() +
    ggtitle("Correlation between G and Gstar") +
    xlab("Regular G") +
    ylab("Gstar by combining with Pedigree") +
    theme_bw()
  
# plot between G and Gstar2 (add small value to diagonal)
  ggplot(Gdata, aes(x=c(G), y=c(Gstar2))) +
    geom_point() +
    ggtitle("Correlation between G and Gstar") +
    xlab("Regular G") +
    ylab("Gstar by adding constant to diagonal") +
    theme_bw()
  
# plot between both Gstars
  ggplot(Gdata, aes(x=c(Gstar), y=c(Gstar2))) +
    geom_point() +
    ggtitle("Correlation between Gstar and Gstar2") +
    xlab("Gstar by combining with Pedigree") +
    ylab("Gstar by adding constant to diagonal") +
    theme_bw()
  
# plot between both Gstars
  ggplot(Gdata, aes(x=c(Gstar), y=c(Gstar2))) +
    geom_point() +
    ggtitle("Correlation between Gstar and Gstar2") +
    xlab("Gstar by combining with Pedigree") +
    ylab("Gstar by adding constant to diagonal") +
    theme_bw()
  
# plot between both Gstars inverses
  ggplot(Gdata, aes(x=c(Ginv), y=c(Ginv2))) +
    geom_point() +
    ggtitle("Correlation between Gstar and Gstar2 inverses") +
    xlab("Gstar inverse (Pedigree)") +
    ylab("Gstar inverse (Diagonal)") +
    theme_bw()
  
# plot between both Gstars inverses
  ggplot(Gdata, aes(x=c(Ginv), y=c(GinvAPY))) +
    geom_point() +
    ggtitle("Correlation between Gstar and G APY inverses") +
    xlab("Gstar inverse (Pedigree)") +
    ylab("G inverse (APY)") +
    theme_bw()
  
  library(reshape2)
  
  rownames(Ginv) <- 1:14
  colnames(Ginv) <- 1:14
  
  rownames(GinvAPY2) <- 1:14
  colnames(GinvAPY2) <- 1:14
  
  Ginv.melt <- melt(Ginv)
  Gapy.melt <- melt(GinvAPY2)
  
# plot regular Ginv (Pedigree)
  ggplot(Ginv.melt, aes(x=as.factor(Var2), y=Var1)) +
  geom_tile(aes(fill=value), colour = "black") +
  scale_y_reverse() +
    geom_text(aes(label = sprintf("%1.2f",value)), vjust = 1, size=3) +
    scale_fill_gradient(low = "white", high = "steelblue")
    ggtitle("G inverse (Pedigree)") +
    theme_bw() 

# plot Ginv APY
  ggplot(Gapy.melt, aes(x=as.factor(Var2), y=Var1)) +
  geom_tile(aes(fill =value), colour = "black") +
  scale_y_reverse() +
    geom_text(aes(label = sprintf("%1.2f",value)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "steelblue")
    ggtitle("G inverse (APY)") +
    theme_bw() 

  
  
  
  
  
  
  
  
  



















