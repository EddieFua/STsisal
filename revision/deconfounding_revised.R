source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/apply.constraints.C.R')
source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/apply.constraints.S.R')
source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/eps.R')
source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/lsqnonneg.col.R')
source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/lsqnonneg.R')
source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/norm.R')
source('/Users/fuyinghao/Documents/异质性分解/deconf package/deconf/R/normalize.col.in.matrix.R')
`deconfounding_revised` <-
  function(I, S.start, n.cell.types,
           n.iterations=1000, error.threshold=0)
  {
    
    #column normalization
    I.noised<-normalize.col.in.matrix(I,col.value=1,method="mean")
    #I <- I.noised$Matrix
    
    #generate start values which fulfil the constraints
    # S.start <- matrix(
    #   sample(seq(from=1,to=(nrow(I) * n.cell.types),by=1),
    #          replace=TRUE,size=nrow(I)* n.cell.types),
    #   nrow=nrow(I),ncol=n.cell.types)

    C.start <- matrix(
      sample(seq( from=1,to=(n.cell.types * ncol(I)),by=0.1),
             replace=TRUE,size= n.cell.types * ncol(I)),
      nrow=n.cell.types,ncol=ncol(I))
    
    return.list <- list()
    #print(n.iterations)
    # alternating: either take S as fix and calculate C or take C
    # as fix and calculate S.
    # start: S is considered as fix
    
    return.list[["S"]] <- apply.constraints.S(S.start)
    return.list[["C"]] <- apply.constraints.C(C.start)
    return.list[["nsim"]] <- 0
    return.list[["error"]] <- Inf
    
    return.list[["error.S"]] <- Inf
    return.list[["error.C"]] <- Inf
    return.list[["min.error"]] <- Inf
    
    return.list[["min.C"]] <- 0
    return.list[["stuck.min.C"]] <- 0
    return.list[["min.S"]] <- 0
    return.list[["stuck.min.S"]] <- 0
    
    exit<- FALSE
    i <- 0
    
    while(exit == FALSE)
    {	i <- i +1
    print("")
    print(paste("iteration",i,"min.error",return.list[["min.error"]]))
    return.list[["nsim"]] <- i
    
    
    if (i %% 2 == 1){	#S fix
      
      x <- lsqnonneg(return.list[["S"]]$Matrix,I.noised$Matrix,return.list[["C"]]$Matrix)
      C2 <- apply.constraints.C(x$x)
      
      if (	x$resnorm >= return.list[["error"]] |
           x$resnorm == return.list[["error.C"]]  ) #no changes, we stuck in local minima
      {
        print("Found minimum for C.")
        return.list[["min.C"]] <- 1
        #but maybe this minimum is not the global minimum ?
        if ( return.list[["error"]] > error.threshold)
        {
          #add some noise
          #return.list[["C"]]$Matrix <- return.list[["C"]]$Matrix[c(2,1),]
          return.list[["stuck.min.C"]] <- 1
          
          return(return.list)
          
        }
      } else {
        return.list[["error.C"]] <- x$resnorm
      }
      return.list[["C"]]$Matrix <- C2$Matrix
      
    }
    
    if (i %% 2 == 0){   #C fix
      
      #switch the C and G matrizes
      x <- lsqnonneg(t(return.list[["C"]]$Matrix), t(I.noised$Matrix),
                     t(return.list[["S"]]$Matrix))
      
      S2 <- apply.constraints.S(t(x$x))
      
      if (	x$resnorm >= return.list[["error"]] |
           x$resnorm == return.list[["error.S"]] ) #no changes, we stuck in local minima
      {
        print("Found minimum for S.")
        return.list[["min.S"]] <- 1
        #but maybe this minimum is not the global minimum ?
        if ( return.list[["error"]] >  error.threshold) {
          #add some noise
          #return.list[["S"]]$Matrix <- return.list[["S"]]$Matrix[,c(2,1)]
          #apply.constraints.S(add.noise(S2$Matrix))
          return.list[["stuck.min.S"]] <- 1
          
          return(return.list)
        }
      } else {
        return.list[["error.S"]] <- x$resnorm
      }
      return.list[["S"]]$Matrix <- S2$Matrix
    }
    
    #determine the error
    return.list[["error"]] <- x$resnorm
    print(paste("error=",return.list[["error"]]))
    
    if (return.list[["min.error"]] > x$resnorm )
      return.list[["min.error"]] <- x$resnorm
    
    if (return.list[["error"]] < error.threshold | i >= n.iterations)
    {
      exit<- TRUE
    }
    }
    return (return.list)
  }

