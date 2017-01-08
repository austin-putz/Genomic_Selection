#==============================================================================#
# APY function
#==============================================================================#

# Author: Austin M. Putz
# Created: Jan 7, 2017
# Modified: Jan 7, 2017
# License: GpLv2

# As always, ABSOLUTELY NO WARRANTY
# This code is simply for learning/illustrative purposes

# INPUTS
#   G    = Genomic relationship matrices
#   core = vector indicating which animals you want to be core

# BEGIN
APY <- function(G, core){
  
  # set n
  n  = dim(G)[1]
  n1 = sum(core)
  n2 = n - n1
  
  core = as.logical(core)
  noncore = !core
  
  # find dimentions of G
  cat("The dimentions of G: ", n, "\n")
  
  cat("The number core: ", n1, "\n")
  cat("The number non-core: ", n2, "\n")
  
  # Split G
  Gcc = G[core, core]
  Gcn = G[core, noncore]
  Gnc = G[noncore, core]
  Gnn = G[noncore, noncore]
  
  # invert Gcc
  Gccinv = solve(Gcc)
  
  # create LEFT side of APY inv
  zeros12 = matrix(0, nrow=n1, ncol=n2)
  zeros21 = matrix(0, nrow=n2, ncol=n1)
  zeros22 = matrix(0, nrow=n2, ncol=n2)
  
  # put together LEFT
  LEFTtop = cbind(Gccinv, zeros12)
  LEFTbot = cbind(zeros21, zeros22)
  LEFT    = rbind(LEFTtop, LEFTbot)
  
  # put together RIGHT
  RIGHTtop = (-1*Gccinv) %*% Gcn
  RIGHTbot = diag(n2)
  RIGHT = rbind(RIGHTtop, RIGHTbot)
  
  # Create M
  Gii = diag(n2) * Gnn
  
  # initialize vector
	Mvec <- matrix(0, ncol=n2, nrow=n2)
	
  # loop through to get elements of M to add to Gii
	for (i in 1:n2) {
	  
	  # subset column of Gcn
	  vec = matrix(Gcn[, i], ncol=1)
	  
	  # do matrix multiplication
	  Mvec[i, i] = (t(vec) %*% Gccinv %*% vec)
	  
	}
	
  # create M by adding Gii to Mvec (diagonal)
	M = Gii + Mvec
  
	# get inverse of M
	Minv = solve(M)
  
	# calculate final APY inverse
	GinvAPY = LEFT + (RIGHT %*% Minv %*% t(RIGHT))
	
	# return value of APY inverse
  return(GinvAPY)
  
}














