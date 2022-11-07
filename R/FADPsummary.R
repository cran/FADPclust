##' @title Summary of FADPclust
##' @description Summarize the result obetained from the FADPclust() function.
##' @param object object of class 'FADPclust' that is returned from FADPclust().
##' @return NULL
##' @seealso \code{\link{FADPclust}}, \code{\link{FADPplot}}.
##' @export
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

FADPsummary = function(object)
{
  p <- length(object)
  if(names(object)[p] == "Avg.silhouette"){Name <- "Average silhouette"}
  if(names(object)[p] == "Dunn"){Name <- "Dunn statistics"}
  if(names(object)[p] == "CH"){Name <- "CH statistics"}
  cat("----------------------\n")
  cat("-- FADPclust Result -- \n")
  cat("----------------------\n \n")
  cat("Number of observations: ", length(object[['clust']]), "\n")
  cat("Parameter selection: \t", object[['method']], "\n")
  cat("Number of clusters: \t", object[['nclust']], "\n")
  cat(Name, ": \t", object[[names(object)[p]]], "\n")
}
