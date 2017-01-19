#==============================================================================#
# createH Function
#==============================================================================#

# start function
  createH <- function(ped, M, wts=c(0.95, 0.05, 1, 1), verbose=FALSE){
    
    #----------------------------------------#
    # INPUTS
    #----------------------------------------#
    
    # M = marker covariate matrix with first column IDs
    # ped = pedigree with 0 or NA as missing
    # wts = vector of (alpha, beta, tau, omega)
    
    #----------------------------------------#
    # HELP
    #----------------------------------------#
    
    # Make sure the pedigree is 3 columns (Animal, Sire, Dam)
    # wts must be a vector of 4
        
    #----------------------------------------#
    # Load libraries
    #----------------------------------------#
    
    # library(data.table)
      library(pedigreemm)
    
    #----------------------------------------#
    # Check and set inputs
    #----------------------------------------#
    
    # check if wts is a 4 length
    if (length(wts) != 4) stop("You need 4 wts (alpha, beta, tau, omega) in that order")
    
    # set inputs
    alpha = wts[1]
    beta  = wts[2]
    tau   = wts[3]
    omega = wts[4]
    
    # show wts
    cat("\nWeights:\n")
    cat("\talpha = ", alpha ,"\n")
    cat("\tbeta = ",  beta ,"\n")
    cat("\ttau = ",   tau ,"\n")
    cat("\tomega = ", omega ,"\n")
    
    # dimentions of M
    nGeno <- nrow(M)
    cat("\nNumber of genotyped individuals: ", nGeno, "\n")
    
    # dimentions of M
    nMarkers <- ncol(M) - 1
    cat("\nNumber of markers: ", nMarkers, "\n\n")
    
    # set M
    list.geno   <- M[, 1]
    rownames(M) <- M[, 1]
    M <- as.matrix(M[, -1])
    
    cat("\nM matrix\n")
    print(M)
    
    # check for duplicates in M matrix
    # if (any(duplicated(list.geno))) stop("\nThere are duplicated genotype IDs\n")
    
    # check dimentions of pedigree
    nPedInd <- nrow(ped)
    nPedCol <- ncol(ped)
    
    # print output
    cat("\nNumber of individuals in the pedigree", nPedInd, "\n")
    cat("\nNumber of columns in the pedigree", nPedCol, "\n\n")
    
    # check to make sure there is only 3 columns in the pedigree
    if (nPedCol != 3) stop("\nThere are not 3 columns in the pedigree, it should be Animal, Sire, Dam\n")
    
    # set column names
    names(ped) <- c("Animal", "Sire", "Dam")
    
    # replace 0's with NA's in Pedigree
    ped[ped==0] <- NA
    
    # add column for genotyped or not in pedigree
    ped$Genotyped <- as.numeric(ped$Animal %in% rownames(M))
    print(ped)
    
    # number of genotyped and non-genotyped
    nGeno    <- sum(ped$Genotyped)
    nNonGeno <- length(ped$Genotyped) - nGeno
    
    # print output
    cat("\nNumber genotyped", nGeno, "\n")
    cat("\nNumber non-genotyped", nNonGeno, "\n\n")
    
    #----------------------------------------#
    # Use pedigree function
    #----------------------------------------#
    
    # editPed() to add parents to top of pedigree
    ped.edit <- editPed(sire=ped$Sire, 
                      dam=ped$Dam, 
                      label=ped$Animal)
    cat("\neditPed from pedigreemm\n")
    print(ped.edit)
    
    # pedigree() function to create pedigree S4 object
    ped.complete <- pedigree(sire= ped.edit$sire, 
                              dam= ped.edit$dam, 
                            label= ped.edit$label)
    cat("\nComplete ped from pedigreemm\n")
    print(ped.complete)
    
    # calculate A
    A <- getA(ped.complete)
    cat("\nA matrix\n")
    print(A)

    # add row and column names to A
    rownames(A) <- ped.complete@label
    colnames(A) <- ped.complete@label
    
    # list of genootyped and non-genotyped animals
    # list.geno    <- ped$Animal[ped$Genotyped==1]  # did up above from M
    list.nongeno <- ped$Animal[ped$Genotyped!=1]
    
    cat("\nList of genotyped animals\n")
    print(list.geno)
    cat("\nList of non-genotyped animals\n")
    print(list.nongeno)
    
    # subset A matrix into parts for non-genotyped (2) and genotyped (1)
    A_11 <- A[rownames(A) %in% list.nongeno, colnames(A) %in% list.nongeno]
    A_22 <- A[rownames(A) %in% list.geno,    colnames(A) %in% list.geno]
    A_12 <- A[rownames(A) %in% list.nongeno, colnames(A) %in% list.geno]
    A_21 <- t(A_12)
    
    cat("\nA11 matrix\n")
    print(A_11)
    cat("\nA12 matrix\n")
    print(A_12)
    cat("\nA2 matrix\n")
    print(A_22)
    
    #----------------------------------------#
    # Create A-1 and subsets
    #----------------------------------------#
    
    # get A inverse
    iA <- solve(A)
    rownames(iA) <- rownames(A)
    colnames(iA) <- colnames(A)
    cat("\nInverse of A\n")
    print(iA)
    
    # subset A matrix into parts for non-genotyped (2) and genotyped (1)
    iA_11 <- iA[rownames(iA) %in% list.nongeno, colnames(iA) %in% list.nongeno]
    iA_12 <- iA[rownames(iA) %in% list.nongeno,    colnames(iA) %in% list.geno]
    iA_22 <- iA[rownames(iA) %in% list.geno,    colnames(iA) %in% list.geno]
    iA_21 <- t(iA_12)
    
    cat("\nInverse A11 matrix\n")
    print(iA_11)
    cat("\nInverse A12 matrix\n")
    print(iA_12)
    cat("\nInverse A22 matrix\n")
    print(iA_22)
    
    #----------------------------------------#
    # Create G
    #----------------------------------------#
    
    cat("\nM matrix\n")
    print(M)
    
    # allele freqs
    allelefreq <- apply(M, 2, mean) / 2
    cat("\nAllele frequencies\n")
    print(allelefreq)
    
    # get P
    P <- matrix(allelefreq*2, byrow=TRUE, ncol=ncol(M), nrow=nrow(M))
    cat("\nP matrix\n")
    print(P)
    
    # Z
    Z <- (M - P)
    cat("\nZ matrix\n")
    print(round(Z, 2))
    
    # G
    sum2pq <- (2*sum(allelefreq*(1-allelefreq)))
    G <- (Z %*% t(Z)) / sum2pq
    # set rownames
    rownames(G) <- rownames(M)
    colnames(G) <- rownames(M)
    cat("\nG matrix\n")
    print(round(G, 3))
    
    # get Gstar
    Gstar <- (alpha*G) + (beta*A_22[order(rownames(G)), order(colnames(G))])
    cat("\nGstar matrix\n")
    print(round(Gstar, 3))

    #----------------------------------------#
    # Sort A inverse and get lower right for H-1
    #----------------------------------------#
    
    # sort A inverse
    iGiA22 <- (tau*solve(Gstar)) - 
                (omega*iA_22[order(rownames(Gstar)), order(colnames(Gstar))])
    rownames(iGiA22) <- rownames(Gstar)
    colnames(iGiA22) <- colnames(Gstar)
    cat("\niGiA22 matrix\n")
    print(round(iGiA22, 3))
    
    # sort A-1
    iH_11 <- iA_11
    iH_12 <- iA_12
    iH_21 <- t(iH_12)
    iH_22 <- iA_22[order(rownames(iGiA22)), order(colnames(iGiA22))] + 
                iGiA22[order(rownames(iGiA22)), order(colnames(iGiA22))]
    
    # top half of inverse H
    topH <- cbind(iH_11, iH_12)
    cat("\nTop H matrix\n")
    print(round(topH, 3))
    
    # bottom half of inverse H
    botH <- cbind(iH_21, iH_22)
    cat("\nBottom H matrix\n")
    print(round(botH, 3))
    
    # inverse H
    iH   <- rbind(topH, botH)
    cat("\nInverse H matrix\n")
    print(round(iH, 3))

    #----------------------------------------#
    # Return the inverse of H
    #----------------------------------------#
    
    # return
    return(iH)
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  







