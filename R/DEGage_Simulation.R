
#'@title QC_Simulation
#'@description
#'Performs quality control on user input for DEGage_simulation
#'@param ngenes is the number of genes to be simulated total
#'@param ndegs is the number of genes to be differentially expressed
#'@param cellgroups is factor with integers corresponding to cell conditions
#'@param lfc is either an integer or vector with a length of ndegs which dictates log fold change value to adjust means for differentially expressed genes
#'@param prop.zeros The proportion of 0's to introduce to simulate dropouts. Must be a value between 0 and 1.
#'@param seed sets a seed for random generation. If user does not input one, a random seed is generated
#'@param ncores the number of cores to allocated for parallel computing
QC_simulation <- function(ngenes, ndegs, cellgroups, lfc = 1, prop.zeros = .3, seed = NULL, ncores = 4){
  if(!inherits(ngenes, "numeric")) stop("'ngenes' must be numeric")
  if(!inherits(ndegs, "numeric")) stop("'ndegs' must be numeric")
  if(!inherits(cellgroups, "factor")) stop("'cellgroups' must be a factor. Try factor(cellgroups) instead.")
  if(prop.zeros > 1 || prop.zeros < 0) stop("'prop.zeros' must be a value between 0 and 1")
}

#'@title gen_rownames
#'@description
#'A helper function that generates row names for simulation.
#' @param ngenes the number of total genes
#' @param ndegs the number of desired differentially expressed genes
#' @param ncores the number of cores to allocated for parallel computing
gen_rownames <- function(ngenes, ndegs, ncores){
  cl.cores <- detectCores()
  cl <- makeCluster(getOption("cl.cores", ncores))  #allocate cores to use

  Degnames <- parSapply(cl = cl,
                        X = c(1:ndegs),
                        FUN = function(x) paste("DEG",x, sep = ""))
  nondegnames <- parSapply(cl = cl,
                           X = c((ndegs+1): ngenes),
                           FUN = function(x) paste("GENE",x, sep = "" ))
  rowtitles <- append(Degnames, nondegnames)
  stopCluster(cl)
  return(rowtitles)
}

#'@title gen_colnames
#'@description
#' A helper function that generates a name for each cell for simulation. Names indicate which condition/group each cell belongs to
#' @param ncells the number of cells to be generated
#' @param ngroup1 the number of cells in group one
#' @param ncores the number of cores to allocated for parallel computing
gen_colnames <- function(ncells, ngroup1, ncores){
  cl.cores <- detectCores()
  cl <- makeCluster(getOption("cl.cores", ncores))  #allocate cores to use

  group1.titles <- parSapply(cl = cl,
                             X = 1:ngroup1,
                             FUN = function(x) paste("Cell",x,".Group1", sep = "" ))
  group2.titles <- parSapply(cl = cl,
                             X = (ngroup1+1):ncells,
                             FUN = function(x) paste("Cell",x,".Group2", sep = "" ))
  celltitles <- append(group1.titles, group2.titles)
  stopCluster(cl)
  return(celltitles)
}

#'@title gen_colnames
#'@description
#' A helper function that introduces dropouts at a predetermined volume
#' @param finaldf simulated count matrix
#' @param ncells the number of cells simulated
#' @param ndropouts the number of dropouts to be introduced
introduce_dropouts <- function(finaldf, ncells, ndropouts){
  dropout.index <- sample(c(1:ncells),ndropouts)
  finaldf[dropout.index] = 0
  return(finaldf)
}

#'@title DEGage_Simulation
#'@description
#'Simulates counts. DEG's are simulated under two conditions
#'@param ngenes is the number of genes to be simulated total
#'@param ndegs is the number of genes to be differentially expressed
#'@param cellgroups is factor with integers corresponding to cell conditions
#'@param lfc is either an integer or vector with a length of ndegs which dictates log fold change value to adjust means for differentially expressed genes
#'@param prop.zeros The proportion of 0's to introduce to simulate dropouts. Must be a value between 0 and 1.
#'@param seed sets a seed for random generation. If user does not input one, a random seed is generated
#'@param ncores the number of cores to allocated for parallel computing
#'@param dispersions a list of dispersions to use to override default generation
#'@param means a list of mean to use to override default generation
#'@export
DEGage_Simulation <- function(ngenes, ndegs, cellgroups, lfc = 1, prop.zeros = .3, seed = NULL, ncores = 4, dispersions = NULL, means = NULL){

  QC_simulation(ngenes, ndegs, cellgroups, lfc = 1, prop.zeros = .3, seed = NULL, ncores = 4)

  #Setting the seed and other useful variables
  if(is.null(seed)){
    seed <- sample(1:100000, 1)
  }
  set.seed(seed)

  ncells = length(cellgroups)
  ngroup1 <- sum(cellgroups == levels(cellgroups)[1])
  ngroup2 <- sum(cellgroups == levels(cellgroups)[2])
  finaldf <- data.frame()

  #Detecting cores for parallel computing
  cl.cores <- detectCores()
  if(cl.cores < ncores){
    stop("'ncores' is greater than the number of available cores. Change ncores to a smaller value")
  }


  #Name generation
  rowtitles <- gen_rownames(ngenes+ndegs, ndegs, ncores)
  celltitles <- gen_colnames(ncells, ngroup1, ncores)

  cl <- makeCluster(getOption("cl.cores", ncores))  #allocate cores to use

  #Parameter Sampling. means are sampled from a gamma distrbution,
  #then multiplied by a random scale factor.
  if(is.null(dispersions)){
    size.list <- rgamma(ngenes+ndegs,1)
  }else{
    size.list <- dispersions
  }

  if(is.null(means)){
    scale.factor.list <- sample(0:100,(ngenes+ndegs), replace = TRUE)
    mean.list <- floor(rgamma(ngenes+ndegs, 1) * scale.factor.list)
  }else{
    if(length(means) == 1){
      mean.list <- rep(means, ngenes+ndegs)
    }else{
      mean.list <- means
    }
  }

  #DEG simulation.
  DEGmeans <- mean.list[1:ndegs]
  if(length(size.list == 1)){
    size.list <- rep(size.list, ndegs+ngenes)
  }
  DEGmeans <- data.frame(means = DEGmeans, disp = size.list[1:ndegs])
  #Group one simulation: simulates counts based off a negative binomial distrubution with predetermined means.
  #Dispersions are sampled from a gamma distrbution

  DEGgroup1df <- data.frame()
  for(i in 1:nrow(DEGmeans)){
    x <- rnbinom(n = ngroup1, size = DEGmeans[i, ]$disp, mu = DEGmeans[i, ]$means)
    DEGgroup1df <- rbind(DEGgroup1df, x)
  }

  #Group two simulation: simulated in the same manner as group one simulation, except means are multiplied by specified
  #LFC value
  Group2means.adj <- DEGmeans$means * (2^lfc)
  Group2means.adj <- data.frame(means = Group2means.adj, disp = size.list[1:ndegs])

  DEGgroup2df <- data.frame()
  for(i in 1:nrow(Group2means.adj)){
    x <- rnbinom(n = ngroup2, size = Group2means.adj[i, ]$disp, mu = Group2means.adj[i, ]$means)
    DEGgroup2df <- rbind(DEGgroup2df, x)
  }

  finaldf <- cbind(DEGgroup1df, DEGgroup2df)

  #Non DEG simulation: genecounts are simulated as one large group with no change in means.

  nondegdf <- data.frame(means = mean.list[(1:ngenes)],
                         disp = size.list[(ndegs + 1):(ngenes+ndegs)])

  nondegcounts <- data.frame()
  for(i in 1:nrow(nondegdf)){
    x <- rnbinom(n = ngroup2 + ngroup1, size = nondegdf[i, ]$disp, mu = nondegdf[i, ]$means)
    nondegcounts <- rbind(nondegcounts, x)
  }

  #Output format
  colnames(nondegcounts) <- celltitles
  colnames(finaldf) <- celltitles
  finaldf <- rbind(finaldf, nondegcounts)
  rownames(finaldf) <- rowtitles

  #dropout management
  ndropouts = ncells*prop.zeros
  finaldf <- parApply(cl=cl, X = finaldf,
                      MARGIN = 1,
                      FUN = introduce_dropouts,
                      ncells = ncells,
                      ndropouts = ndropouts)
  finaldf <- as.data.frame(t(as.data.frame(finaldf)))

  #Closing parallel computing cluster
  stopCluster(cl)

  return(finaldf)
}

