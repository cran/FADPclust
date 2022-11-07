##' @title Visualize the result of FADPclust
##' @description Plot the f vs delta plot with selected centroids.
##' @param object object of class 'FADPclust' that is returned from FADPclust().
##' @param cols vector of colors used to distinguish different clusters. Ten default colors are given.
##' @return NULL
##' @seealso \code{\link{FADPclust}}, \code{\link{FADPsummary}}.
##' @importFrom graphics points
##' @importFrom graphics text
##' @importFrom stats dist
##' @export
##'
##' @examples
##' ###univariate functional data
##' data("simData1")
##' plot(simData1, xlab = "x", ylab = "y")
##' FADP1.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP1",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        stats = "Avg.silhouette")
##' FADPsummary(FADP1.ans); FADPplot(FADP1.ans)
##' \donttest{
##' FADP2.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP2",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        pve = 0.9, stats = "Avg.silhouette")
##' FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
##'
##' ###multivariate functional data
##' data("simData2")
##' FADP1.ans <- FADPclust(fdata = simData2, cluster = 2:5, method = "FADP1",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        pve = 0.9, stats = "Avg.silhouette")
##' FADPsummary(FADP1.ans); FADPplot(FADP1.ans)
##'
##' FADP2.ans <- FADPclust(fdata = simData2, cluster = 2:5, method = "FADP2",
##'                        proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
##'                        pve = 0.9, stats = "Avg.silhouette")
##' FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
##' }


FADPplot <- function(object, cols = "default") {
  nclusters <- sils <- NULL  # Null out to remove 'no visible binding for global variable' note from R check.
  defCol    <- function() {
    mycols  <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                 "#A65628", "#F781BF", "#999999", "blue")
    return(mycols)
  }
  if (cols == "default")
    cols <- defCol()

  # Recycle colors
  if ((temp <- ceiling(object$nclust/length(cols))) > 1)
    cols    <- rep(cols, temp)[1:object$nclust]

  f       <- object[["density"]]
  delta   <- object[["delta"]]
  centers <- object[["center"]]

  ##--------------------
  ## f vs delta
  ##--------------------
  plot(f, delta, xlab = "f(x)", ylab = "delta(x)", main = "f(x) vs delta(x)")
  f.range     <- range(f)
  delta.range <- range(delta)
  points(f[centers], delta[centers], col = cols, pch = 19, cex = 1.2)
  text(f[centers], delta[centers], labels = centers, cex = 0.6, pos = 1)

}

