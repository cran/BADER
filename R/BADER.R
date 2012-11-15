BADER <-
function (
    kA, kB,
    sizeFactors=TRUE,
    start=NULL,
    burn=1e3, reps=1e4,
    printEvery=1e2, saveEvery=1,
    t0 = 1e1,
    mode = 0 )
{
    #### check input ####
    if ( is.null(dim(kA)) | is.null(dim(kB)))
        stop("kA and kB should be matrices")

    if ( dim(kA)[1] > dim(kA)[2] | dim(kB)[1] > dim(kB)[2])
        warning("kA and kB have more rows than columns. In a typical usage example where
            rows represent samples and columns genes the contrary should be the case")

    if ( dim(kA)[2] != dim(kB)[2])
        stop("kA and kB need to have the same number of columns")

    if ( !is.null(start) & typeof(start) != "list")
        stop("The 'start' parameter needs to be either NULL or a list of starting parameters")

    #### initialize helper variables ####
    nA <- dim(kA)[1]
    nB <- dim(kB)[1]
    m <- dim(kA)[2]
   
    #### estimate size factors ####
    if ( sizeFactors == FALSE )
    {
        sA <- rep(1, nA)
        sB <- rep(1, nB)
    }
   
    if ( sizeFactors == TRUE )
    {
        data <- rbind(kA,kB)
        geoMeans <- exp(colMeans(log(data)))
        sfac <- apply(data,1,function(x)
            median( (x / geoMeans)[geoMeans>0]))

        sA <- sfac[1:nA]
        sB <- sfac[(nA+1):(nA+nB)]
    }

    #### guess starting parameters for MCMC ####    
    ## use simple method-of-moments estimators ##  
    lambdaA <- ifelse(kA/sA > 0, log(kA/sA), NA)
    lambdaA[is.na(lambdaA)] <- min(lambdaA, na.rm=TRUE)

    lambdaB <- ifelse(kB/sB > 0, log(kB/sB), NA)
    lambdaB[is.na(lambdaB)] <- min(lambdaB, na.rm=TRUE)

    muA <- colMeans(lambdaA)
    muB <- colMeans(lambdaB)

    gam <- muB - muA

    pi0 <- 0.1
    ind <- rep(0,m)
    cutoff <- sort(abs(gam))[round((1-pi0)*m)] 
    ind[abs(gam) > cutoff] <- 1
    gam[ind == 0] <- 0
    sigmaGamma <- var(gam[gam!=0])
  
    varA <- apply(kA/sA,2,var)
    varB <- apply(kB/sB,2,var)

    alphaA <- rep(NA,m)
    alphaB <- rep(NA,m)
    alphaA[varA-exp(muA) > 0] <- log(((varA - exp(muA))*exp(-2*muA))[varA-exp(muA) > 0])
    alphaB[varB-exp(muB) > 0] <- log(((varB - exp(muB))*exp(-2*muB))[varB-exp(muB) > 0])

    tau <- var(c(alphaA,alphaB), na.rm=TRUE)
    alphaA[is.na(alphaA)] <- mean(alphaA, na.rm=TRUE)
    alphaB[is.na(alphaB)] <- mean(alphaB, na.rm=TRUE)


    psi0 <- mean(c(alphaA,alphaB))
    
    #### use start parameters if present ####
    if ( typeof(start) == "list")
    {
        params <- c("lambdaA", "lambdaB", "muA", "gam", "ind", 
            "alphaA", "alphaB", "psi0", "tau", "sigmaGamma", "pi0")

        for ( param in intersect(names(start),params) )
        {
            if ( !is.null(dim(get(param))) && dim(get(param)) != dim(start[[param]]))
                stop(paste("Starting parameter '",param,"' does not have the correct dimensions",sep=""))
            if ( length(get(param)) != length(start[[param]]))
                stop(paste("Starting parameter '",param,"' does not have the correct length",sep=""))
            
            assign(param,start[[param]])        
        }
    }

    #### perform MCMC ####
    if ( mode == 0 )
    {
        results <- .C("rnaseq", as.double(kA), as.double(kB), as.double(sA), as.double(sB),
            as.integer(nA), as.integer(nB), as.integer(m), as.integer(burn), as.integer(reps),
            as.integer(saveEvery), as.integer(printEvery), as.integer(t0),
            as.double(lambdaA), as.double(lambdaB), as.double(ind), as.double(muA),
            as.double(gam), as.double(alphaA), as.double(alphaB), as.double(pi0),
            as.double(sigmaGamma), as.double(psi0), as.double(tau)
        )
    }
    if ( mode == 1 )
    {
        distGamma <- rep(0,floor(reps/saveEvery)*m)
        results <- .C("rnaseq_post_dist", as.double(kA), as.double(kB), as.double(sA), as.double(sB),
            as.integer(nA), as.integer(nB), as.integer(m), as.integer(burn), as.integer(reps),
            as.integer(saveEvery), as.integer(printEvery), as.integer(t0),
            as.double(lambdaA), as.double(lambdaB), as.double(ind), as.double(muA),
            as.double(gam), as.double(alphaA), as.double(alphaB), as.double(pi0),
            as.double(sigmaGamma), as.double(psi0), as.double(tau),
            as.double(distGamma)
        )
    }
    if ( mode == 2 )
    {
        distGamma <- distMuA <- distAlphaA <- distAlphaB <- rep(0,floor(reps/saveEvery)*m)
        distPi0 <- distSigmaGamma <- distPsi0 <- distTau <- rep(0,floor(reps/saveEvery))
        
        results <- .C("rnaseq_verbose", as.double(kA), as.double(kB), as.double(sA), as.double(sB),
            as.integer(nA), as.integer(nB), as.integer(m), as.integer(burn), as.integer(reps),
            as.integer(saveEvery), as.integer(printEvery), as.integer(t0),
            as.double(lambdaA), as.double(lambdaB), as.double(ind), as.double(muA),
            as.double(gam), as.double(alphaA), as.double(alphaB), as.double(pi0),
            as.double(sigmaGamma), as.double(psi0), as.double(tau),
            as.double(distGamma), as.double(distMuA), as.double(distAlphaA), as.double(distAlphaB)
        )
    }

    #### return results ####
    if ( mode == 0 )
    {
        obj <- list (
            logMeanA = results[[16]],
            logMeanB = (results[[16]] + results[[17]]),
            logFoldChange = results[[17]],
            logDispA = results[[18]],
            logDispB = results[[19]],
            diffProb = results[[15]]
        )
        class (obj) <- c(class(obj),"BADER") 
    }
    if ( mode == 1 )
    {
        obj <- list (
            logMeanA = results[[16]],
            logMeanB = (results[[16]] + results[[17]]),
            logFoldChange = results[[17]],
            logDispA = results[[18]],
            logDispB = results[[19]],
            diffProb = results[[15]],
            logFoldChangeDist = matrix(results[[24]],ncol=m)
        )
        class (obj) <- c(class(obj),"BADER") 
    }
    if ( mode == 2 )
    {
        obj <- list (
            logMeanA = results[[16]],
            logMeanB = (results[[16]] + results[[17]]),
            logFoldChange = results[[17]],
            logDispA = results[[18]],
            logDispB = results[[19]],
            diffProb = results[[15]],
            logFoldChangeDist = matrix(results[[24]],ncol=m),
            logMeanADist = matrix(results[[25]],ncol=m),
            logMeanBDist = matrix(results[[24]] + results[[25]],ncol=m),
            logDispADist = matrix(results[[26]],ncol=m),
            logDispBDist = matrix(results[[27]],ncol=m)
        )
        class (obj) <- c(class(obj),"BADER") 
    }

    return(obj)
}
