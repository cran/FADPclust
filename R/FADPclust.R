##' @title Functional Data Clustering Using Adaptive Density Peak Detection
##' @description Clustering of univariate or multivariate functional data by finding cluster centers from estimated density peaks. FADPclust is a non-iterative procedure that incorporates KNN density estimation algorithm. The number of clusters can also be selected by the user or selected automatically through an internal clustering criterion.
##' @param fdata for univariate functional data clustering: a functional data object produced by fd() function of fda package; for multivariate functional data clustering: a list of functional data objects produced by fd() function of fda package.
##' @param cluster integer, or a vector of integers specifying the pool of the number of clusters in automatic variation. The default is 2:10.
##' @param method character string specifying the method used to calculate the pseudo functional k-nearest neighbor density. Valid options of are 'FADP1' and 'FADP2' (see details in references). The default is 'FADP1'.
##' @param proportion numeric, a number or numeric vector of numbers within the range [0,1], specifying to automatically select the smoothing parameter k in density estimation (see details). The default is 0.1, 0.2, ... ,1.
##' @param f.cut numeric, a number within the range [0,1], specified to automatically select cluster centroids from the decision plot. The default is 0.15.
##' @param pve numeric, a number within the range [0,1], the proportion of variance explained: used to choose the number of functional principal components. The default is 0.99. When the method is chosen to be 'FADP1', there is no need to specify parameter 'pve' for univariate functional data clustering.
##' @details Given n functional objects or curves, FADPclust() calculates f(x) and delta(x) for each object based on the semi-metric distance (see details in references), where f(x) is the local density calculated by the functional k-nearest neighbor density estimator of curve x, and delta(x) is the shortest semi-metric distance between sample curve x and y for all samples y such that f(x) <= f(y). Functional objects or curves with large f and large delta values are labeled class centroids. In other words, they appear as isolated points in the upper right corner of the f vs delta plot (the decision plot, see details in FADPplot). After cluster centroids are determined, other obejects are clustered according to their semi-metric distances to the closes centroids.
##'
##' The smoothing parameter k in functional k-nearest neighbor density estimation must be explicitly provided. Following Lauter (1988)'s idea, suggest that the optimal size of k satisfies a certain proportion, k = a*n^(4/5), where a is a parameter about the optimal proportion to be determined. Here, users enters variable 'proportion' to specify the parameter a.
##' @return An 'FADPclust' object that contains the list of the following items.
##' \itemize{
##' \item{nclust:}{ number of clusters. }
##' \item{para:}{ smoothing parameter k selected automatically by KNN estimation.}
##' \item{method:}{ character string introducing the method used to calculate the smoothing parameter. }
##' \item{clust:}{ cluster assignments. A vector of the same length as the number of observations. }
##' \item{density:}{ final density vector f(x). }
##' \item{delta:}{ final delta vector delta(x). }
##' \item{center:}{ indices of the clustering centers. }
##' \item{silhouette:}{ silhouette score from the final clustering result. }
##' }
##' @references
##' \itemize{
##' \item Lauter, H. (1988), "Silverman, B. W.: "Density Estimation for Statistics and Data Analysis.," Biometrical Journal, 30(7), 876-877.
##' \item Wang, X. F., and Xu, Y. (2016), "Fast Clustering Using Adaptive Density Peak Detection," Statistical Methods in Medical Research.
##' \item Rodriguez, A., and Laio, A. (2014), "Machine learning. Clustering by fast search and find of density peaks," Science, 344(6191), 1492.
##' \item Yaohui, L., Zhengming, M., and Fang, Y. (2017), "Adaptive density peak clustering based on K-nearest neighbors with aggregating strategy," Knowledge-Based Systems, 133(oct.1), 208-220.
##' }
##' @seealso \code{\link{FADPsummary}}, \code{\link{FADPplot}}.
##' @examples
##' ###univariate functional data
##' data("simData1")
##' plot(simData1, xlab = "x", ylab = "y")
##' FADP1.ans <- FADPclust(fdata = simData1, cluster = 2:10, method = "FADP1",
##'                        proportion = seq(0.02, 0.2, 0.02))
##' FADP2.ans <- FADPclust(fdata = simData1, cluster = 2:10, method = "FADP2",
##'                      proportion = seq(0.02, 0.2, 0.02), pve = 0.9)
##' FADPsummary(FADP1.ans); FADPplot(FADP1.ans)
##' FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
##'
##' \donttest{
##' ###multivariate functional data
##' data("simData2")
##' FADP1.ans <- FADPclust(fdata = simData2, cluster = 2:10, method = "FADP1",
##'                        proportion = seq(0.02, 0.2, 0.02), pve = 0.9)
##' FADP2.ans <- FADPclust(fdata = simData2, cluster = 2:10, method = "FADP2",
##'                      proportion = seq(0.02, 0.2, 0.02), pve = 0.9)
##' FADPsummary(FADP1.ans); FADPplot(FADP1.ans)
##' FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
##' }
##' @import fda
##' @import cluster
##' @import MFPCA
##' @import funData
##' @importFrom stats dist
##' @importFrom fda.usc metric.lp
##' @importFrom fda.usc semimetric.basis
##' @export

FADPclust <- function(fdata, cluster = 2:10, method = "FADP1", proportion = seq(0.1,1,0.1), f.cut = 0.15, pve = 0.99 )
{
  if( class(fdata) != "fd" & class(fdata) != "list" ){
    stop( "Error in 'fdata'! The input of 'fdata' should be a functional data object produced by fd() function of fda package for univariate FADP, a list of functional data objects for multivariate FADP.", sep="\n" )
  }

  if( class(proportion) != "numeric" ){
    stop( "Error in 'proportion'! The input of 'proportion' should be a number or numeric vector of numbers within the range [0,1].", sep="\n" )
  }

  if( min(proportion) < 0 | max(proportion) > 1){
    stop( "Error in 'proportion'! The input value of 'proportion' should be between 0 and 1.", sep="\n" )
  }

  if( method != "FADP1" & method != "FADP2" ){
    stop( "Error in 'method'! The input of 'method' is wrong, valid options of are 'FADP1' and 'FADP2'.", sep="\n" )
  }

  if( class(f.cut) != "numeric" ){
    stop( "Error in 'f.cut'! The input of 'f.cut' should be a number between 0 and 1.", sep="\n" )
  }

  if( f.cut <= 0 | f.cut > 1){
    stop( "Error in 'f.cut'! The input value of 'f.cut' should be between 0 and 1.", sep="\n" )
  }

  if( pve <= 0 | pve > 1){
    stop( "Error in 'pve'! The input value of 'pve' should be between 0 and 1.", sep="\n" )
  }

  ### Related functions ###
  knn_density1 <- function(distance, k){
    distmat    <- as.matrix(distance)
    n          <- nrow(distmat)
    r.k        <- apply(distmat, 1, sort)[k+1,]
    den        <- k / (n * r.k)
    return(den)
  }

  knn_density2 <- function(score, k) {
    n <- nrow(score); p <- ncol(score)
    if (k >= n) {stop("k is not a reasonable and valid number!")}
    distance <- dist(score, method = "euclidean", upper = TRUE)
    distmat <- as.matrix(distance)
    v.d <- pi^(p/2) /gamma(p/2+1)
    r.k <- apply(distmat, 1, sort)[k+1,]
    f_hat <- k / (n * v.d * r.k)
    eta.k <- apply(distmat, 1, sort)[k + 1, ]
    mu.k  <- mean(eta.k)
    h     <- mu.k + sqrt(sum((eta.k - mu.k)^2)/(n - 1))
    den   <- c()
    k   <- c()
    for (i in 1:nrow(distmat)) {
      k[i] <- sum(as.vector(distmat[i, ]) <= h) - 1
    }
    den <- k/(n * f_hat * (h)^p)
    return(den)
  }

  adp1 <- function(distance, clusters) {
    distance  <- as.matrix(distance)
    n.curve   <- nrow(distance)
    k.list    <- unique(ceiling(proportion * n.curve^(4/5)))
    k.len     <- length(k.list)
    clu.len   <- length(clusters)
    den_del   <- array(data = NA, dim = c(2, n.curve, k.len))
    sil_temp  <- matrix(data = NA, nrow = k.len, ncol = clu.len)
    for (i in 1:k.len) {
      den     <- knn_density1(distance = distance, k = k.list[i])
      del     <- c()
      den.max <- which(den == max(den), arr.ind = TRUE)
      for (j in 1:n.curve) {
        if (j %in% den.max) {
          del[j] <- max(distance[j, ])
        } else {
          del[j] <- min(distance[j, which(den[j] < den)])
        }
      }
      den_del[, , i] <- rbind(den, del)
      ### Selection Criterion
      clu.order <- 1
      alpha     <- c()
      for (m in clusters) {
        if(round(n.curve*f.cut) < m){
          alpha[clu.order] <- m/n.curve
        }else{ alpha[clu.order] = f.cut }
        den.index              <- which(den %in% sort(den, decreasing = T)[1:round(n.curve*alpha[clu.order])])
        center.temp            <- which(del %in% sort(del[den.index], decreasing = T)[1:m])
        result.temp            <- apply(distance[center.temp, ], 2, which.min)
        sil                    <- summary(silhouette(result.temp, dmatrix = distance))[["avg.width"]]
        sil_temp[i, clu.order] <- sil
        clu.order              <- clu.order + 1
      }
    }
    index          <- which(sil_temp == sil_temp[which.max(sil_temp)], arr.ind = T)
    index          <- index[1, ]
    den            <- den_del[1, , index[1]]
    del            <- den_del[2, , index[1]]
    nclust         <- clusters[index[2]]
    if(round(n.curve*f.cut) < nclust){
      alpha_opt <- nclust/n.curve
    }else{ alpha_opt = f.cut }
    den.index      <- which( den %in% sort(den, decreasing = T)[1:round(n.curve*alpha_opt)] )
    center.index   <- which(del %in% sort(del[den.index], decreasing = T)[1:nclust])
    res_clustering <- as.numeric(apply(distance[center.index, ], 2, which.min))
    result         <- list(length(unique(center.index)), round(k.list[index[1]]), method, res_clustering, den, del, center.index, max(sil_temp))
    names(result)  <- c("nclust", "para","method", "clust", "density", "delta", "center", "silhouette")
    return(result)
  }

  adp2 <- function(distance, clusters, score) {
    distance  <- as.matrix(distance)
    n.curve   <- nrow(distance)
    k.list    <- unique(ceiling(proportion * n.curve^(4/5)))
    k.len     <- length(k.list)
    clu.len   <- length(clusters)
    den_del   <- array(data = NA, dim = c(2, n.curve, k.len))
    sil_temp  <- matrix(data = NA, nrow = k.len, ncol = clu.len)
    for (i in 1:k.len) {
      den     <- knn_density2(score = score, k = k.list[i])
      del     <- c()
      den.max <- which(den == max(den), arr.ind = TRUE)
      for (j in 1:n.curve) {
        if (j %in% den.max) {
          del[j] <- max(distance[j, ])
        } else {
          del[j] <- min(distance[j, which(den[j] < den)])
        }
      }
      den_del[, , i] <- rbind(den, del)
      ### Selection Criterion
      clu.order <- 1
      alpha     <- c()
      for (m in clusters) {
        if(round(n.curve*f.cut) < m){
          alpha[clu.order] <- m/n.curve
        }else{ alpha[clu.order] = f.cut }
        den.index              <- which(den %in% sort(den, decreasing = T)[1:round(n.curve*alpha[clu.order])])
        center.temp            <- which(del %in% sort(del[den.index], decreasing = T)[1:m])
        result.temp            <- apply(distance[center.temp, ], 2, which.min)
        sil                    <- summary(silhouette(result.temp, dmatrix = distance))[["avg.width"]]
        sil_temp[i, clu.order] <- sil
        clu.order              <- clu.order + 1
      }
    }
    index          <- which(sil_temp == sil_temp[which.max(sil_temp)], arr.ind = T)
    index          <- index[1, ]
    den            <- den_del[1, , index[1]]
    del            <- den_del[2, , index[1]]
    nclust         <- clusters[index[2]]
    if(round(n.curve*f.cut) < nclust){
      alpha_opt <- nclust/n.curve
    }else{ alpha_opt = f.cut }
    den.index      <- which( den %in% sort(den, decreasing = T)[1:round(n.curve*alpha_opt)] )
    center.index   <- which(del %in% sort(del[den.index], decreasing = T)[1:nclust])
    res_clustering <- as.numeric(apply(distance[center.index, ], 2, which.min))
    eta.k.opt      <- apply(distance, 1, sort)[k.list[index[1]] + 1, ]
    mu.k.opt       <- mean(eta.k.opt)
    h.opt          <- mu.k.opt + sqrt(sum((eta.k.opt - mu.k.opt)^2)/(n.curve - 1))
    result         <- list(length(unique(center.index)), round(k.list[index[1]]), method, res_clustering, den, del, center.index, max(sil_temp))
    names(result)  <- c("nclust", "para", "method", "clust", "density", "delta", "center", "silhouette")
    return(result)
  }

  if( class(fdata) == "fd" & method == "FADP1" ){
    nbasis      <- fdata$basis$nbasis
    type        <- fdata$basis$type
    distance_L2 <- semimetric.basis(fdata1 = fdata, fdata2 = fdata, nderiv = 0, type.basis1 = type,
                                    nbasis1 = nbasis, type.basis2 = type, nbasis2 = nbasis)#L2 distance
    result      <- adp1(distance = distance_L2, clusters = cluster)
  }

  if( class(fdata) == "list" & method == "FADP1" ){
    basis    <- fdata[[1]]$basis
    nbasis   <- fdata[[1]]$basis$nbasis
    rangeval <- fdata[[1]]$basis$rangeval
    p        <- length(fdata)
    n.curve  <- ncol(fdata[[1]]$coefs)
    dot      <- nbasis*10
    t        <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - rangeval[1])/(dot - 1))
    data     <- array(0, dim = c(n.curve, dot, p))
    for (i in 1:p) {
      data[, , i] <- t(eval.fd(evalarg = t, fdobj = fdata[[i]]))
    }

    ### Data standardization
    for (i in 1:p) {
      m1          <- min(data[, , i])
      m2          <- max(data[, , i])
      data.temp   <- (data[, , i] - m1)/(m2 - m1)
      data[, , i] <- data.temp
    }

    ### Step 1: calculate the eigenfunctions
    lam <- matrix(0, p, dot)  #eigenvalues
    a   <- array(data = 0, dim = c(p, dot, p))  #eigenvectors
    for (i in 1:dot) {
        xir      <- matrix(0, p, n.curve)
      for (r in 1:p) {
        xir[r, ] <- data[, i, r]
      }
      sigma    <- xir %*% t(xir)
      out      <- eigen(sigma)
      lam_i    <- out[["values"]]
      a_i      <- out[["vectors"]]
      lam[, i] <- lam_i
      for (r in 1:p) {
        a[, i, r] <- a_i[, r]
      }
    }

    #################### Step 2: adjust the signs of eigenvectors
    w    <- 5 * (t[2] - t[1])  #bandwidth
    l    <- c()
    l[1] <- 1
    for (k in 2:dot) {
      if ((t[k] - w) < 1) {
        l[k] <- 1
      } else {
        l[k] <- k - 5
      }
    }
    judge <- function(r, a) {
      # a is the original eigenvectors
      adj <- c()
      for (k in 2:dot) {
        sum1 <- 0
        sum2 <- 0
        for (j in l[k]:k - 1) {
          sum1 <- sum1 + sqrt(sum((a[, j, r] - a[, k, r])^2))
          sum2 <- sum2 + sqrt(sum((a[, j, r] + a[, k, r])^2))
        }
        if (sum1 > sum2) {
          adj[k] <- 1
        }
      }
      return(adj)
    }

    for (r in 1:p) {
      loc         <- which(judge(r, a) == 1)
      a[, loc, r] <- -a[, loc, r]
    }

    #################### Step 3: calculate the variability explained
    v  <- lam[1, ] + lam[2, ]
    phi <- c()
    # calculate by phi_2
    for (k in 1:p) {
      phi_temp <- sum(lam[k, ] * 0.1)/sum(v * 0.1)
      phi[k]   <- phi_temp
    }
    m <- 1
    s <- 0
    repeat {
      s <- s + phi[m]
      m <- m + 1
      # 0.8 is the upper limit of variability explained
      if (s >= pve) {
        num <- m - 1
        break
      }
    }
    # print(num) #the number of principal components

    #################### Step 4:calculate the r-th principal components
    component <- function(r) {
        bind.matrix <- data[, , 1]
      for (m in 2:p) {
        bind.matrix <- cbind(bind.matrix, data[, , m])
      }
      z   <- c()
      zir <- matrix(0, n.curve, dot)
      for (i in 1:n.curve) {
        xi <- matrix(bind.matrix[i, ], p, dot, byrow = TRUE)
        for (j in 1:dot) {
          z[j] <- t(a[, j, r]) %*% xi[, j]
        }
        zir[i, ] <- z
      }
      return(zir)
    }

    #################### ADPclust the distance of principal components
    comp_list      <- list()
    beta           <- phi[1:num]/sum(phi[1:num])
    distance_list2 <- list()
    distance_L2    <- matrix(0, n.curve, n.curve)
    for (r in 1:num) {
      comp                <- component(r)
      compfd              <- Data2fd(argvals = t, y = t(comp), basisobj = basis)
      L2_distance         <- semimetric.basis(fdata1 = compfd, fdata2 = compfd, nderiv = 0,
                                      type.basis1 = "bspline", nbasis1 = nbasis, type.basis2 = "bspline",
                                      nbasis2 = nbasis)
      comp_list[[r]]      <- compfd
      distance_list2[[r]] <- L2_distance
      distance_L2         <- distance_L2 + distance_list2[[r]]^2 * beta
    }
    distance <- sqrt(distance_L2)
    result   <- adp1(distance = distance, clusters = cluster)
  }

  if( class(fdata) == "fd" & method == "FADP2" ){
    nbasis   <- fdata$basis$nbasis
    rangeval <- fdata$basis$rangeval
    dot      <- nbasis*10
    t        <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - rangeval[1])/(dot - 1))
    fundata  <- funData(argvals = t, X = t(eval.fd(evalarg = t, fdobj = fdata)))
    pca  <- PACE(fundata, pve = pve)
    score       <- pca$scores
    distance_L2 <- dist(score, method = "euclidean", upper = TRUE)

    result      <- adp2(distance = distance_L2, clusters = cluster, score = score)
  }

  if( class(fdata) == "list" & method == "FADP2" ){
    nbasis   <- fdata[[1]]$basis$nbasis
    rangeval <- fdata[[1]]$basis$rangeval
    dot      <- nbasis*10
    t        <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - rangeval[1])/(dot - 1))
    fundata  <- list()
    type     <- list()
    for(i in 1:length(fdata)){
      fundata[[i]]  <- funData(argvals = t, X = t(eval.fd(evalarg = t, fdobj = fdata[[i]])) )
      type[[i]]     <- list(type = "uFPCA")
    }
    comp_num   <- PACE(fundata[[1]], pve = pve)$npc
    fundata <- multiFunData(fundata)
    pca  <- MFPCA(mFData = fundata, M = comp_num, uniExpansions = type)
    score       <- pca$scores
    distance_L2 <- dist(score, method = "euclidean", upper = TRUE)

    result      <- adp2(distance = distance_L2, clusters = cluster, score = score)
  }

  return(result)
}

