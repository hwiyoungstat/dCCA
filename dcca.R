dCCA <- function(X,Y,Z,lambda){
  
  Max.Ite <- 1000
  EPS <- 1e-4
  step <- 1
  
  Group1 <- which(Z==0)
  Group2 <- which(Z==1)
  
  P <- ncol(X)
  Q <- ncol(Y)
  N <- dim(X)[1]
  
  # Initialize
  covX <- diag(P)
  covY <- diag(Q)
  Sig <- cov(X,Y)
  Sig1 <- cov(X[Group1,],Y[Group1,])
  Sig2 <- cov(X[Group2,],Y[Group2,])
    
  # Sxeig <- eigen(covX,symmetric=TRUE)
  # invsqrtcovX <- Sxeig$vectors %*% diag(1/sqrt(Sxeig$values)) %*% t(Sxeig$vectors)
  # 
  # Syeig <- eigen(covY,symmetric=TRUE)
  # invsqrtcovY <- Syeig$vectors %*% diag(1/sqrt(Syeig$values)) %*% t(Syeig$vectors)
  # 
  # 
  # 
  # Xmat <- invsqrtcovX %*% Sig %*% solve(covY) %*% t(Sig) %*% invsqrtcovX
  # Ymat <- invsqrtcovY %*% t(Sig) %*% solve(covX) %*% Sig %*%invsqrtcovY
  # 
  # Xeig <- eigen(Xmat,symmetric = TRUE)
  # Yeig <- eigen(Ymat,symmetric = TRUE)
  # 
  # Xc <- invsqrtcovX %*% Xeig$vectors
  # Yc <- invsqrtcovY %*% Yeig$vectors
  # 
  # u.initial <- Xc[,1]
  # v.initial <- Yc[,1]

  
  
  u.initial <- rep(1,P)
  v.initial <- rep(1,Q)
  
  # Iteration dCCA
  
  u.old <- u.initial
  v.old <- v.initial
  a.old <- 1
  b.old <- 1
  
  for(i in 1:Max.Ite){
#    print(i)
   sgn <- as.numeric(sign(t(u.old)%*%(Sig1-Sig2)%*%v.old))
   u.new <- u.old + step*(c((Sig+lambda*sgn*(Sig1-Sig2))%*%v.old)+a.old*u.old)
   v.new <- v.old + step*(c((t(Sig)+lambda*sgn*t(Sig1-Sig2))%*%u.old)+b.old*v.old)
   
   a.new <- a.old + step/2*(sum(u.old^2)-1)
   b.new <- b.old + step/2*(sum(v.old^2)-1)
  
  # projection
   u.new <- u.new/sqrt(sum(u.new^2))
   v.new <- v.new/sqrt(sum(v.new^2))
  
  # Convergence Evaluation
   if(sum((u.new-u.old)^2)<EPS & sum((v.new-v.old)^2)<EPS){break}
   
   u.old <- u.new
   v.old <- v.new
   
   a.old <- a.new
   b.old <- b.new
  }
  
  result <- list("u"=u.new,"v"=v.new)
  
}