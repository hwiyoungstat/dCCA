Screen <- function(w,thres){
  
  Lambda.N <- 100
  Lambda <- seq(0,1,length.out=Lambda.N)
  KL.Div <- rep(0,length(Lambda))
  #Thres <- 0.2
  Thres <- thres
  A <- 1*(abs(w)>Thres)
  p <- dim(w)[1]
  q <- dim(w)[2]
  pi <- sum(A)/(p*q)
  
  
  Sx <- seq(1,p)
  Sy <- seq(1,q)
  
  Sx.new <- seq(1,p)
  Sy.new <- seq(1,q)
  
  Sx.list <- list(c(p+q-2))
  Sy.list <- list(c(p+q-2))
  
  Sx.opt.list <- list(Lambda.N)
  Sy.opt.list <- list(Lambda.N)
  
  density <- matrix(0,ncol=(p+q-2),nrow=length(Lambda))
  
  for(t in 1:(p+q-2)){
    
    cat("\r", t ,"of", (p+q-2))
    KL.X <- rep(0,length(Sx.new))
    KL.Y <- rep(0,length(Sy.new))
    #    A.new <- matrix(A[Sx.new,Sy.new],nrow=length(Sx.new),ncol=length(Sy.new))
    A.new <- abs(matrix(w[Sx.new,Sy.new],nrow=length(Sx.new),ncol=length(Sy.new)))
    
    if(length(Sx.new)==1){Degx<-Inf}else{
      for(i in 1:length(Sx.new)){
        KL.X[i] <- sum(A.new[i,])
      }
      v <- which.min(KL.X)
      Degx <- mean(A.new[v,])
    }
    
    if(length(Sy.new)==1){Degy<-Inf}else{
      for(j in 1:length(Sy.new)){
        KL.Y[j] <- sum(A.new[,j])
      }
      u <- which.min(KL.Y)
      Degy <- mean(A.new[,u])
    }
    
    ratio <- length(Sx.new)/length(Sy.new)
    ratio <- 1
    
    
    ifelse(ratio*Degx <= 1/ratio*Degy, Sx.new <- Sx.new[-v], Sy.new <- Sy.new[-u])
    Sx.list[[t]] <- Sx.new
    Sy.list[[t]] <- Sy.new
    
    
    for(l in 1:length(Lambda)){
      lambda <- Lambda[l]
      density[l,t] <- sum(abs(w[Sx.new,Sy.new]))/(length(Sx.new)*length(Sy.new))^lambda
      density[l,t] <- sum(abs(A[Sx.new,Sy.new]))/(length(Sx.new)*length(Sy.new))^lambda
    }
  }
  
  for(l in 1:length(Lambda)){
    opt <- which.max(density[l,])
    Sx.hat <- Sx.list[[opt]]
    Sy.hat <- Sy.list[[opt]]
    
    w1 <- c(A[Sx.hat,Sy.hat])
    w0 <- c(c(A[,-Sy.hat]),c(A[-Sx.hat,Sy.hat]))
    
    pi1 <- sum(w1)/length(w1)
    pi0 <- sum(w0)/length(w0)
    
    
    
    temp1 <- sum(w1*pi1*log(pi1/pi)+(1-w1)*(1-pi1)*log((1-pi1)/(1-pi)))
    temp2 <- sum(w0*pi0*log(pi0/pi)+(1-w0)*(1-pi0)*log((1-pi0)/(1-pi)))
    KL.Div[l] <- temp1 + temp2
    
    Sx.opt.list[[l]] <- Sx.hat
    Sy.opt.list[[l]] <- Sy.hat
  }
  
  opt.L <- which.max(KL.Div)
  Lambda.opt <- Lambda[opt.L]
  
  
  Sx.hat <- Sx.opt.list[[opt.L]]
  Sy.hat <- Sy.opt.list[[opt.L]]
  
  
  
  mylist <- list("density"=density,"Lambda"=Lambda.opt,"Lambdas"=Lambda,"KL"=KL.Div,"Sx"=Sx.hat,"Sy"=Sy.hat,"Sx.list"=Sx.opt.list,"Sy.list"=Sy.opt.list)
  return(mylist)
}