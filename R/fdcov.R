#' Analysis of Covariance Operators.
#'
#' \code{fdcov} a variety of tools for the analysis of covariance 
#' operators including k-sample tests for equality and 
#' classification and clustering methods.
#'
#' This package contains a collection of tools for performing
#' statistical inference on functional data specifically through
#' an analysis of the covariance structure of the data.  It includes
#' two methods for performing a k-sample test for equality
#' of covariance in \code{ksample.perm} and \code{ksample.com}
#' and two methods for 2-sample tests for equality assuming 
#' Gaussian data in \code{ksample.gauss} and \code{ksample.vstab}.
#' For supervised and unsupervised learning,
#' it contains a method to classify functional data with
#' respect to each category's covariance operator in
#' \code{classif.com},
#' and it contains a method to cluster functional data, \code{cluster.com},
#' again based on the covariance structure of the data.
#'
#' The current version of this package assumes that all functional
#' data is sampled on the same grid at the same intervals.  Future
#' updates are planned to allow for the below methods to interface
#' with the \code{fda} package and its functional basis representations
#' of the data.
#'
#' @author
#' Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it},
#' Adam B Kashlak \email{kashlak@ualberta.ca}
#' @references
#'   Cabassi, A., Pigoli, D., Secchi, P., Carter, P. A. (2017). 
#'   Permutation tests for the equality of covariance operators of 
#'   functional data with applications to evolutionary biology. 
#'   Electron. J. Statist. 11(2), pp.3815--3840.
#'
#'   Kashlak, Adam B, John AD Aston, and Richard Nickl (2016).
#'   "Inference on covariance operators via concentration
#'   inequalities: k-sample tests, classification, and clustering via
#'   Rademacher complexities", April, 2016 (in review)
#'
#'   Pigoli, Davide, John AD Aston, Ian L Dryden, and Piercesare Secchi.
#'   "Distances and inference for covariance operators."
#'   Biometrika (2014): 101(2):409–422.
#'
#'   Panaretos, Victor M., David Kraus, and John H. Maddocks.
#'   "Second-order comparison of Gaussian random functions and the
#'   geometry of DNA minicircles." Journal of the American Statistical
#'   Association 105.490 (2010): 670-682.
#' @importFrom grDevices dev.off postscript setEPS
#' @importFrom stats cov qnorm rbinom rgamma sd pchisq
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
"_PACKAGE"
#> [1] "_PACKAGE"

