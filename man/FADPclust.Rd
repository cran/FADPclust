% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FADPclust.R
\name{FADPclust}
\alias{FADPclust}
\title{Functional Data Clustering Using Adaptive Density Peak Detection}
\usage{
FADPclust(fdata, cluster = 2:10, method = "FADP1", proportion = NULL,
  f.cut = 0.15, pve = 0.9, stats = "Avg.silhouette")
}
\arguments{
\item{fdata}{for univariate functional data clustering: a functional data object produced by fd() function of fda package; for multivariate functional data clustering: a list of functional data objects produced by fd() function of fda package.}

\item{cluster}{integer, or a vector of integers specifying the pool of the number of clusters in automatic variation. The default is 2:10.}

\item{method}{character string specifying the method used to calculate the pseudo functional k-nearest neighbor density. Valid options of are 'FADP1' and 'FADP2' (see details in references). The default is 'FADP1'.}

\item{proportion}{numeric, a number or numeric vector of numbers within the range [0,1], specifying to automatically select the smoothing parameter k in density estimation (see details). The default is 0.1, 0.2, ... ,1.}

\item{f.cut}{numeric, a number within the range [0,1], specified to automatically select cluster centroids from the decision plot. The default is 0.15.}

\item{pve}{numeric, a number within the range [0,1], the proportion of variance explained: used to choose the number of functional principal components. The default is 0.9. When the method is chosen to be 'FADP1', there is no need to specify parameter 'pve' for univariate functional data clustering.}

\item{stats}{character string specifying the distance based statistics for cluster validation and determining the number of clusters. Valid options are 'Avg.silhouette', 'Dunn', and 'CH' (See the description document of the cluster.stats function in the fpc R package for more details about these statistics). The default is "Avg.silhouette".}
}
\value{
An 'FADPclust' object that contains the list of the following items.
\itemize{
\item{nclust:}{ number of clusters. }
\item{para:}{ smoothing parameter k selected automatically by KNN estimation.}
\item{method:}{ character string introducing the method used to calculate the smoothing parameter. }
\item{clust:}{ cluster assignments. A vector of the same length as the number of observations. }
\item{density:}{ final density vector f(x). }
\item{delta:}{ final delta vector delta(x). }
\item{center:}{ indices of the clustering centers. }
\item{Avg.silhouette:}{ average silhouette score from the final clustering result. }
\item{Dunn:}{ Dunn statistics from the final clustering result. }
\item{CH:}{ CH statistics from the final clustering result. }
}
}
\description{
Clustering of univariate or multivariate functional data by finding cluster centers from estimated density peaks. FADPclust is a non-iterative procedure that incorporates KNN density estimation algorithm. The number of clusters can also be selected by the user or selected automatically through an internal clustering criterion.
}
\details{
Given n functional objects or curves, FADPclust() calculates f(x) and delta(x) for each object based on the semi-metric distance (see details in references), where f(x) is the local density calculated by the functional k-nearest neighbor density estimator of curve x, and delta(x) is the shortest semi-metric distance between sample curve x and y for all samples y such that f(x) <= f(y). Functional objects or curves with large f and large delta values are labeled class centroids. In other words, they appear as isolated points in the upper right corner of the f vs delta plot (the decision plot, see details in FADPplot). After cluster centroids are determined, other obejects are clustered according to their semi-metric distances to the closes centroids.

The smoothing parameter k in functional k-nearest neighbor density estimation must be explicitly provided. Following Lauter (1988)'s idea, suggest that the optimal size of k satisfies a certain proportion, k = a*n^(4/5), where a is a parameter about the optimal proportion to be determined. Here, users enters variable 'proportion' to specify the parameter a.
}
\examples{
###univariate functional data
data("simData1")
plot(simData1, xlab = "x", ylab = "y")
FADP1.sil.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP1",
                           proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                           stats = "Avg.silhouette")
FADPsummary(FADP1.sil.ans); FADPplot(FADP1.sil.ans)
\donttest{
FADP1.dunn.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP1",
                            proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                            stats = "Dunn")
FADPsummary(FADP1.dunn.ans); FADPplot(FADP1.dunn.ans)

FADP1.ch.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP1",
                          proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                          stats = "CH")
FADPsummary(FADP1.ch.ans); FADPplot(FADP1.ch.ans)

FADP2.ans <- FADPclust(fdata = simData1, cluster = 2:5, method = "FADP2",
                       proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                       pve = 0.9, stats = "Avg.silhouette")
FADPsummary(FADP2.ans); FADPplot(FADP2.ans)

###multivariate functional data
data("simData2")
FADP1.ans <- FADPclust(fdata = simData2, cluster = 2:5, method = "FADP1",
                       proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                       pve = 0.9, stats = "Avg.silhouette")
FADPsummary(FADP1.ans); FADPplot(FADP1.ans)

FADP2.ans <- FADPclust(fdata = simData2, cluster = 2:5, method = "FADP2",
                       proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                       pve = 0.9, stats = "Avg.silhouette")
FADPsummary(FADP2.ans); FADPplot(FADP2.ans)
}
}
\references{
\itemize{
\item Lauter, H. (1988), "Silverman, B. W.: "Density Estimation for Statistics and Data Analysis.," Biometrical Journal, 30(7), 876-877.
\item Wang, X. F., and Xu, Y. (2016), "Fast Clustering Using Adaptive Density Peak Detection," Statistical Methods in Medical Research.
\item Rodriguez, A., and Laio, A. (2014), "Machine learning. Clustering by fast search and find of density peaks," Science, 344(6191), 1492.
\item Liu Y, Ma Z, and Yu F. (2017), "Adaptive density peak clustering based on K-nearest neighbors with aggregating strategy," Knowledge-Based Systems, 133(oct.1), 208-220.
}
}
\seealso{
\code{\link{FADPsummary}}, \code{\link{FADPplot}}.
}
