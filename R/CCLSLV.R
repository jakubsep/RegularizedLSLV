##1 Cardinality constrained LSLV
CCLSLV_base <- function(DATA, R, CARD, MaxIter=200, eps=10^-6) {
  J <- dim(DATA)[2]
  I <- dim(DATA)[1]
  CARD = rep(CARD, R)
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc        #Used to compute convergence criterium
  loadings <- INITLOADINGS(DATA, R, CARD)
  A <- loadings
  scores <- matrix(nrow = I, ncol = R)
  if (sum(length(CARD)) > 1){
    CARDc <- J-CARD #Number of loadings -per PC- that have to be zero
    while (convAO == 0) {
      #1. Update component scores
      XP <- DATA %*% loadings/ (sqrt(I))
      scores <- sqrt(I) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
      A <- t(DATA) %*% scores/ I
      #Calculate loss
      Loss <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Loss)
      #2. Update loadings
      loadings <- A
      for (r in 1:R){
        ind <- sort(abs(A[,r]), index.return = TRUE)
        loadings[ind$ix[1:CARDc[r]],r] <- 0
      }
      #Calculate loss
      Lossu <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Lossu)
      if (iter > MaxIter) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }} else {
      CARDc <- J*R-CARD #Number of loadings -over all PCs- that have to be zero
      while (convAO == 0) {
        #1. Update component scores
        XP <- DATA %*% loadings / (sqrt(I))
        scores <- sqrt(I) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A <- t(DATA) %*% scores / I
        #Calculate loss
        Loss <- RESIDUAL(DATA,scores,loadings)
        Lossvec <- c(Lossvec,Loss)
        #2. Update loadings
        loadings <- A
        ind <- sort(abs(A), index.return = TRUE)
        loadings[ind$ix[1:CARDc]] <- 0
        #Calculate loss
        Lossu <- RESIDUAL(DATA,scores,loadings)
        Lossvec <- c(Lossvec,Lossu)
        if (iter > MaxIter ) {
          convAO <- 1
        }
        if (abs(Lossc-Lossu) < eps) {
          convAO <- 1
        }
        iter <- iter + 1
        Lossc <- Lossu
      }
    }
  uslpca <- list('scores' = scores, 'loadings' = loadings, 'Lossvec' = Lossvec, 'Residual' = Loss)
  return(uslpca)
}

##2. LOSS
RESIDUAL <- function(DATA, SCORES, LOADINGS){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))

  return(res)
}

##3. INTIAL LOADINGS
INITLOADINGS <- function(DATA, R, CARD){
  J <- dim(DATA)[2]
  I <- dim(DATA)[1]
  svd1 <- svd(DATA, R, R)
  P1 <- matrix(rnorm(J*R), ncol = R, nrow = J)
  P2 <- svd1$v %*% diag(svd1$d[1:R])/sqrt(I)
  P <- P1*.3 + P2*.7
  if (length(CARD)>1){
    CARDc <- J-CARD
    for (r in 1:R){
      ind <- sort(abs(P[,r]), index.return = TRUE)
      P[ind$ix[1:CARDc[r]],r] <- 0
    }
  } else {
    CARDc <- J*R-CARD
    ind <- sort(abs(P), index.return = TRUE)
    P[ind$ix[1:CARDc]] <- 0
  }
  loadings <- P
  return(loadings)
}

###function for recovery rate
num_correct <- function (TargetP, EstimatedP){
  total_vnumber <- dim(TargetP)[1] * dim(TargetP)[2]
  TargetP[which(TargetP != 0)] <- 1
  sum_select <- sum(TargetP)
  sum_zero <- total_vnumber - sum_select
  EstimatedP[which(EstimatedP != 0)] <- 1
  total_correct <- sum(TargetP == EstimatedP) # this is the total number of variables correctedly selected and zeros correctly retained
  prop_correct <- total_correct/total_vnumber
  return(prop_correct)
}

#########################################################################################
#####                             Multistart procedure                              #####
#########################################################################################
CCLSLV <- function(DATA, R, CARD, MaxIter=200, eps=10^-6, nstarts=25){

  if(R == 1){
    stop("Parameter R = 1 is not allowed.")
  }

  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()

  for (n in 1:nstarts){
    result <- CCLSLV_base(DATA, R, CARD, MaxIter, eps)

    Pout3d[[n]] <- result$loadings
    Tout3d[[n]] <- result$scores
    LOSS[n] <- result$Residual
    LOSSvec[[n]] <- result$Lossvec
  }
  # choose solution with lowest loss value
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }

  return_varselect <- list()
  return_varselect$loadings <- Pout3d[[k]]
  return_varselect$scores <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  return_varselect$Loss <- LOSS[k]

  return(return_varselect)
}


# Index of sparsity
IS_CC <- function(data, R, MaxIter=200, eps=10^-6, nstarts=10){
  IS_all <- matrix(ncol = 5, nrow = ncol(DATA), dimnames=list(c(), c('card', 'isvalue', 'vaf', 'propzero', 'propunique'))) # Matrix to keep solutions for all cardinalities
  variables <- R*ncol(data)
  for (j in 1:ncol(data)){
    #1. Obtain CCLSLV solution
    tested_card = ncol(data)-j+1
    CCLSLV <- CCLSLV(data, R, CARD = tested_card, MaxIter, eps, nstarts)
    #2. Calculate BIC
    ssrespca <- 0
    J <- dim(data)[2]
    va <- 0
    nrzero <- 0
    vzero <- 0
    #calculation SSRES PCA
    vzero <- sum(rowSums(data^2))
    d <- svd(data)$d
    ssrespca <- ssrespca+sum((d[(R+1):length(d)])^2)
    loading <- c(round(CCLSLV$loadings,3))
    va <- va+sum((d[1:R])^2)
    nrzero <- nrzero+sum(loading==0)
    dfjspca <- 0
    for (jr in 1:(J*R)){
      dfjspca <- dfjspca+length(unique(loading[jr][loading[jr]!=0]))
    }
    bicvalue <- CCLSLV$Loss/ssrespca+log(dfjspca)
    vs <- vzero-CCLSLV$Loss
    nrcoef <- J*R
    nrzeqcoef <- nrcoef-dfjspca
    IS_all[j,] <- c(tested_card, va*vs/vzero^2*nrzeqcoef/nrcoef, 1-(CCLSLV$Loss/vzero), dfjspca/(nrcoef-nrzero), nrzero/nrcoef)
  }
  return(list(value=IS_all[which.max(IS_all[, 2])], matrix=IS_all))
}
