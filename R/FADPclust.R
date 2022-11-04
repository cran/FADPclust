##' @title Functional Data Clustering Using Adaptive Density Peak Detection
##' @description Clustering of univariate or multivariate functional data by finding cluster centers from estimated density peaks. FADPclust is a non-iterative procedure that incorporates KNN density estimation algorithm. The number of clusters can also be selected by the user or selected automatically through an internal clustering criterion.
##' @param fdata for univariate functional data clustering: a functional data object produced by fd() function of fda package; for multivariate functional data clustering: a list of functional data objects produced by fd() function of fda package.
##' @param cluster integer, or a vector of integers specifying the pool of the number of clusters in automatic variation. The default is 2:10.
##' @param method character string specifying the method used to calculate the pseudo functional k-nearest neighbor density. Valid options of are 'FADP1' and 'FADP2' (see details in references). The default is 'FADP1'.
##' @param proportion numeric, a number or numeric vector of numbers within the range [0,1], specifying to automatically select the smoothing parameter k in density estimation (see details). The default is 0.1, 0.2, ... ,1.
##' @param f.cut numeric, a number within the range [0,1], specified to automatically select cluster centroids from the decision plot. The default is 0.15.
##' @param pve numeric, a number within the range [0,1], the proportion of variance explained: used to choose the number of functional principal components. The default is 0.9. When the method is chosen to be 'FADP1', there is no need to specify parameter 'pve' for univariate functional data clustering.
##' @param stats character string specifying the distance based statistics for cluster validation and determining the number of clusters. Valid options are 'silhouette', 'Dunn', and 'CH' (See the description document of the cluster.stats function in the fpc R package for more details about these statistics). The default is "silhouette".
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
##' \item Liu Y, Ma Z, and Yu F. (2017), "Adaptive density peak clustering based on K-nearest neighbors with aggregating strategy," Knowledge-Based Systems, 133(oct.1), 208-220.
##' }
##' @seealso \code{\link{FADPsummary}}, \code{\link{FADPplot}}.
##' @examples
##' ###univariate functional data
##' data("simData1")
##' plot(simData1, xlab = "x", ylab = "y")
##' FADP1.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP1",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        stats = "silhouette")
##' FADPsummary(FADP1.ans); FADPplot(FADP1.ans)
##' \donttest{
##' FADP2.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP2",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        pve = 0.9, stats = "silhouette")
##' FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
##'
##' ###multivariate functional data
##' data("simData2")
##' FADP1.ans <- FADPclust(fdata = simData2, cluster = 2:5, method = "FADP1",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        pve = 0.9, stats = "silhouette")
##' FADPsummary(FADP1.ans); FADPplot(FADP1.ans)
##'
##' FADP2.ans <- FADPclust(fdata = simData2, cluster = 2:5, method = "FADP2",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        pve = 0.9, stats = "silhouette")
##' FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
##' }
##' @import fda
##' @import fpc
##' @import cluster
##' @import MFPCA
##' @import funData
##' @importFrom stats dist
##' @importFrom fda.usc metric.lp
##' @importFrom fda.usc semimetric.basis
##' @export

FADPclust <- function(fdata, cluster = 2:10, method = "FADP1",
                      proportion = NULL, f.cut = 0.15,
                      pve = 0.90, stats = "silhouette"){
  if( !inherits(fdata, "fd") & !inherits(fdata, "list") ){
    stop( "Error in fdata! For univariate FADP: fdata should be a functional data object produced by fd() function of fda package, for multivariate FADP: a list of functional data objects.", sep="\n" )
  }

  if( is.null(proportion) & method == "FADP1" ){ proportion <- seq(0.1,1,0.1) }

  if( is.null(proportion) & method == "FADP2" ){  proportion <- seq(0.1,1,0.1) }

  ### Related functions ###
  knn_density1 <- function(distance, k){
    distmat    <- as.matrix(distance)
    n          <- nrow(distmat)
    v.d        <- pi^(1/2) /gamma(1/2+1)
    r.k        <- apply(distmat, 1, sort)[k+1,]
    den        <- k / (n * v.d * r.k)
    return(den)
  }

  knn_density2 <- function(score, k) {
    if (is.null(ncol(score)) == FALSE) {
      n <- nrow(score); p <- ncol(score)
    } else {n <- length(score); p <- 1}
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

  assignment <- function(center, distance, den){
    n.curve   <- nrow(distance)
    assign <- rep(0, n.curve)
    assign[center] <- 1:length(center)
    assigned.index <- center
    temp.index <- center
    den.max <- which.max(den)
    while (length(which(assign == 0)) != 0) {
      loc <- rep(0, n.curve)
      for (j in setdiff(1:n.curve, assigned.index)) {
        if (j %in% den.max) {
          nearest.neighbor <- center[which.min(distance[center, j])]
        } else {
          neighbors.index <- which(den[j] < den)
          nearest.neighbor <- neighbors.index[which.min(distance[j, neighbors.index])]
        }
        loc[j] <- nearest.neighbor
      }
      for (l in 1:length(temp.index)) {
        loc.index <- which(loc == temp.index[l])
        assign[loc.index] = rep(assign[temp.index[l]], length(loc.index))
      }
      assigned.index <- which(assign != 0)
      temp.index <- setdiff(assigned.index, temp.index)
    }
    return(assign)
  }

  adp1 <- function(distance, clusters){
    distance  <- as.matrix(distance)
    n.curve   <- nrow(distance)
    k.list    <- unique(ceiling(proportion * n.curve^(4/5)))

    lapply.adp1 <- function(para.index){
      k <- paralist[para.index,1]
      m <- paralist[para.index,2]
      den     <- knn_density1(distance = distance, k = k)
      del     <- c()
      den.max <- which(den == max(den), arr.ind = TRUE)
      del <- mapply(function(j) if (j %in% den.max) { del[j] <- max(distance[j, ]) } else { del[j] <- min(distance[j, which(den[j] < den)]) }, 1:n.curve)
      if(round(n.curve*f.cut) < m){
        alpha <- m/n.curve
      }else{ alpha = f.cut }
      den.index              <- which(den %in% sort(den, decreasing = T)[1:round(n.curve*alpha)])
      center.temp            <- which(del %in% sort(del[den.index], decreasing = T)[1:m])
      result.temp            <- assignment(center = center.temp, distance = distance, den = den)
      stats.temp             <- cluster.stats(d = distance, clustering = result.temp)
      if(stats == "silhouette"){ evaluation <- stats.temp$avg.silwidth }
      if(stats == "Dunn"){ evaluation <- stats.temp$dunn2 }
      if(stats == "CH"){ evaluation <- stats.temp$ch }
      return(list(evaluation = evaluation, clustering = result.temp, density = den, delta = del, center = center.temp))
    }

    paralist <- expand.grid(k.list, clusters)
    task.len <- nrow(paralist)
    result.temp <- lapply(1:task.len, lapply.adp1)
    evaluation.list <- unlist(lapply(1:task.len, function(l) result.temp[[l]]$evaluation))
    selected.index <- which.max(evaluation.list)
    k.selected <- paralist[selected.index, 1]
    nclust.selected <- paralist[selected.index, 2]
    clustering.selected <- result.temp[[selected.index]]$clustering
    density.selected <- result.temp[[selected.index]]$density
    delta.selected <- result.temp[[selected.index]]$delta
    center.selected <- result.temp[[selected.index]]$center
    result         <- list(nclust.selected, proportion[which(k.list == k.selected)],
                           method, clustering.selected, density.selected, delta.selected, center.selected, max(evaluation.list))
    #result         <- list(nclust.selected, round(k.selected),
    #                       method, clustering.selected, density.selected, delta.selected, center.selected, max(evaluation.list))
    names(result)  <- c("nclust", "para","method", "clust", "density", "delta", "center", "stats")
    return(result)
  }

  adp2 <- function(distance, clusters, score){
    distance  <- as.matrix(distance)
    n.curve   <- nrow(distance)
    k.list    <- unique(ceiling(proportion * n.curve^(4/5)))

    lapply.adp2 <- function(para.index){
      k <- paralist[para.index,1]
      m <- paralist[para.index,2]

      den     <- knn_density2(score = score, k = k)
      del     <- c()
      den.max <- which(den == max(den), arr.ind = TRUE)
      del <- mapply(function(j) if (j %in% den.max) { del[j] <- max(distance[j, ]) } else { del[j] <- min(distance[j, which(den[j] < den)]) }, 1:n.curve)
      if(round(n.curve*f.cut) < m){
        alpha <- m/n.curve
      }else{ alpha = f.cut }
      den.index              <- which(den %in% sort(den, decreasing = T)[1:round(n.curve*alpha)])
      center.temp            <- which(del %in% sort(del[den.index], decreasing = T)[1:m])
      result.temp            <- assignment(center = center.temp, distance = distance, den = den)
      stats.temp             <- cluster.stats(d = distance, clustering = result.temp)
      if(stats == "silhouette"){ evaluation <- stats.temp$avg.silwidth }
      if(stats == "Dunn"){ evaluation <- stats.temp$dunn2 }
      if(stats == "CH"){ evaluation <- stats.temp$ch }
      return(list(evaluation = evaluation, clustering = result.temp, density = den, delta = del, center = center.temp))
    }
    paralist <- expand.grid(k.list, clusters)
    task.len <- nrow(paralist)
    result.temp <- lapply(1:task.len, lapply.adp2)

    evaluation.list <- unlist(lapply(1:task.len, function(l) result.temp[[l]]$evaluation))
    selected.index <- which.max(evaluation.list)
    k.selected <- paralist[selected.index, 1]
    nclust.selected <- paralist[selected.index, 2]
    clustering.selected <- result.temp[[selected.index]]$clustering
    density.selected <- result.temp[[selected.index]]$density
    delta.selected <- result.temp[[selected.index]]$delta
    center.selected <- result.temp[[selected.index]]$center
    eta.k.opt      <- apply(distance, 1, sort)[k.selected + 1, ]
    mu.k.opt       <- mean(eta.k.opt)
    h.opt          <- mu.k.opt + sqrt(sum((eta.k.opt - mu.k.opt)^2)/(n.curve - 1))

    result         <- list(nclust.selected, proportion[which(k.list == k.selected)],
                           method, clustering.selected, density.selected, delta.selected, center.selected, max(evaluation.list))
    #result         <- list(nclust.selected, h.opt,
    #                       method, clustering.selected, density.selected, delta.selected, center.selected, max(evaluation.list))
    names(result)  <- c("nclust", "para","method", "clust", "density", "delta", "center", "stats")
    return(result)
  }

  if( inherits(fdata, "fd") & method == "FADP1" ){
    nbasis      <- fdata$basis$nbasis
    type        <- fdata$basis$type
    distance_L2 <- semimetric.basis(fdata1 = fdata, fdata2 = fdata, nderiv = 0, type.basis1 = type,
                                    nbasis1 = nbasis, type.basis2 = type, nbasis2 = nbasis)#L2 distance
    result      <- adp1(distance = distance_L2, clusters = cluster)
  }

  if( inherits(fdata, "list") & method == "FADP1" ){
    basis    <- fdata[[1]]$basis
    nbasis   <- fdata[[1]]$basis$nbasis
    rangeval <- fdata[[1]]$basis$rangeval
    p        <- length(fdata)
    n.curve  <- ncol(fdata[[1]]$coefs)
    dot      <- nbasis*10
    t        <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - rangeval[1])/(dot - 1))
    # dot=nbasis+1
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
    v  <- apply(lam, 2, sum)#lam[1, ] + lam[2, ]
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

  if( inherits(fdata, "fd") & method == "FADP2" ){
    nbasis   <- fdata$basis$nbasis
    rangeval <- fdata$basis$rangeval
    dot      <- nbasis*10
    t        <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - rangeval[1])/(dot - 1))
    fundata  <- funData(argvals = t, X = t(eval.fd(evalarg = t, fdobj = fdata)))
    pca  <- PACE(fundata, pve=pve)
    score       <- pca$scores
    distance_L2 <- dist(score, method = "euclidean", upper = TRUE)

    result      <- adp2(distance = distance_L2, clusters = cluster, score = score)
  }

  if( inherits(fdata, "list") & method == "FADP2" ){
    nbasis   <- fdata[[1]]$basis$nbasis
    rangeval <- fdata[[1]]$basis$rangeval
    dot      <- nbasis*10
    t        <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - rangeval[1])/(dot - 1))
    fundata  <- list()
    type     <- list()
    comp_num <- c()
    for(i in 1:length(fdata)){
      fundata[[i]]  <- funData(argvals = t, X = t(eval.fd(evalarg = t, fdobj = fdata[[i]])) )
      type[[i]]     <- list(type = "uFPCA")
      comp_num[i]   <- PACE(fundata[[i]], pve=pve)$npc
    }
    fundata <- multiFunData(fundata)
    pca  <- MFPCA(mFData = fundata, M = max(comp_num), uniExpansions = type)
    score       <- pca$scores
    distance_L2 <- dist(score, method = "euclidean", upper = TRUE)

    result      <- adp2(distance = distance_L2, clusters = cluster, score = score)
  }

  return(result)
}

