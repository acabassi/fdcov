#' k-sample test for equality of covariance operators
#'
#' \code{ksample.gauss} performs a k-sample test for equality
#' of covariance operators under the assumption that the
#' data arises from a Gaussian process.
#'
#' \code{ksample.vstab} applies a similar method that has
#' been modified to stabilize the variance.  See the reference
#' paper for more details on the mathematics of these methods.
#'
#' These two methods use the Karhunen-Loeve expansion
#' (eigen expansion for functional data) to represent
#' the data in terms of K eigen-functions.  Then a
#' test statistic with asymptotic chi-squared distribution
#' is computed in order to test for the equality of the
#' covariance operators based on the two samples.
#' If K is set to be 0, then the methods determine the
#' number of eigen-functions to retain.
#'
#' @param  dat1 the first set of data with one entry per row
#' @param  dat2 the second set of data with one entry per row
#' @param  K the number of basis vectors to use, Default is 5.
#' @return p-value testing whether or not the two samples have
#' differing covariance operators.
#' @references
#' Panaretos, Victor M., David Kraus, and John H. Maddocks.
#' "Second-order comparison of Gaussian random functions and the
#' geometry of DNA minicircles." Journal of the American Statistical
#' Association 105.490 (2010): 670-682.
#' @author Adam B Kashlak \email{kashlak@ualberta.ca}
#' @examples
#' # Load in phoneme data
#' library(fds)
#' # Set up test data
#' dat1 = t(aa$y)[1:20,];
#' dat2 = t(sh$y)[1:20,];
#' dat3 = t(aa$y)[21:40,];
#' # Compare two disimilar phonemes
#' # Resulting in a small p-value
#' ksample.gauss(dat1,dat2,K=5);
#' ksample.vstab(dat1,dat2,K=5);
#' # Compare two sets of the same phonemes
#' # Resulting in a large p-value
#' ksample.gauss(dat1,dat3,K=5);
#' ksample.vstab(dat1,dat3,K=5);
#' @export
#'
ksample.gauss <- function(
  dat1, dat2, K=5
){
  # Centre data
  dat1 = t(t(dat1)-apply(dat1,2,mean));
  dat2 = t(t(dat2)-apply(dat2,2,mean));
  # group sizes
  n1 = nrow(dat1);
  n2 = nrow(dat2);
  N  = n1+n2;
  # sample covariances
  cov1= (n1-1)/n1*cov(dat1);
  cov2= (n2-1)/n2*cov(dat2);
  covP= n1/N*cov1 + n2/N*cov2;#cov(rbind(dat1,dat2));
  # Eigen decomp
  eigP= eigen(covP);
  # Choose value of K if not provided
  if(K==0){
    eig1 = eigen(cov1);
    eig2 = eigen(cov2);
    pfc_scores = rep(0,min(n1,n2));
    for(K in 2:min(n1,n2)){
      vec  = eigP$vectors[,1:K]; # take first K eigenvectors
      prj1 = ksg_projData( dat1, vec );
      prj2 = ksg_projData( dat2, vec );
      gof  = sum((prj1-dat1)^2) + sum((prj2-dat2)^2);
      ip1  = (prj1 %*% eig1$vectors)^2/eig1$values;
      ip2  = (prj2 %*% eig2$vectors)^2/eig2$values;
      pen  = 2*(
        sum(eig1$values)*sum(ip1)/n1 +
        sum(eig2$values)*sum(ip2)/n2
      );
      pfc_scores[K] = gof + pen;
    }
    pfc_scores[1] = max(pfc_scores);
    K = which.min(pfc_scores);
    #cumEig = cumsum(eigP$values)/sum(eigP$values);
    #K      = which(cumEig>0.9)[1];
  }
  # Get Eigenvectors
  vec = eigP$vectors[,1:K]; # take first K eigenvectors
  # Fourier coefficients
  ev1 = t(vec) %*% cov1 %*% vec/n1;
  ev2 = t(vec) %*% cov2 %*% vec/n2;
  # Compute test statistic
  denom = n1*diag( ev1 ) + n2*diag( ev2 );
  mat   = t((ev1-ev2)^2/denom)/denom;
  tstat = (n1*n2*N)/2*sum(mat);
  dof   = K*(K+1)/2;
  return( pchisq(tstat,dof,lower.tail=FALSE)  );
}

# Project data matrix (N x d) onto vectors (d x K)
ksg_projData <- function( dat, vec ){
  return( (dat %*% vec) %*% t(vec) );
}

#' k-sample test for equality of covariance operators
#'
#' \code{ksample.gauss} performs a k-sample test for equality
#' of covariance operators under the assumption that the
#' data arises from a Gaussian process.
#'
#' \code{ksample.vstab} applies a similar method that has
#' been modified to stabilize the variance.  See the reference
#' paper for more details on the mathematics of these methods.
#'
#' These two methods use the Karhunen-Loeve expansion
#' (eigen expansion for functional data) to represent
#' the data in terms of K eigen-functions.  Then a
#' test statistic with asymptotic chi-squared distribution
#' is computed in order to test for the equality of the
#' covariance operators based on the two samples.
#' If K is set to be 0, then the methods determine the
#' number of eigen-functions to retain.
#'
#' @param  dat1 the first set of data with one entry per row
#' @param  dat2 the second set of data with one entry per row
#' @param  K the number of basis vectors to use, Default is 5.
#' @return p-value testing whether or not the two samples have
#' differing covariance operators.
#' @references
#' Panaretos, Victor M., David Kraus, and John H. Maddocks.
#' "Second-order comparison of Gaussian random functions and the
#' geometry of DNA minicircles." Journal of the American Statistical
#' Association 105.490 (2010): 670-682.
#' @author Adam B Kashlak \email{kashlak@ualberta.ca}
#' @examples
#' # Load in phoneme data
#' library(fds)
#' # Set up test data
#' dat1 = t(aa$y)[1:20,];
#' dat2 = t(sh$y)[1:20,];
#' dat3 = t(aa$y)[21:40,];
#' # Compare two disimilar phonemes
#' # Resulting in a small p-value
#' ksample.gauss(dat1,dat2,K=5);
#' ksample.vstab(dat1,dat2,K=5);
#' # Compare two sets of the same phonemes
#' # Resulting in a large p-value
#' ksample.gauss(dat1,dat3,K=5);
#' ksample.vstab(dat1,dat3,K=5);
#' @export
#'
ksample.vstab <- function(
  dat1, dat2, K=5
){
  # Centre data
  dat1 = t(t(dat1)-apply(dat1,2,mean));
  dat2 = t(t(dat2)-apply(dat2,2,mean));
  # group sizes
  n1 = nrow(dat1);
  n2 = nrow(dat2);
  N  = n1+n2;
  # sample covariances
  cov1= (n1-1)/n1*cov(dat1);
  cov2= (n2-1)/n2*cov(dat2);
  covP= n1/N*cov1 + n2/N*cov2;#cov(rbind(dat1,dat2));
  # Eigen decomp
  eigP= eigen(covP);
  # Choose value of K if not provided
  if(K==0){
    cumEig = cumsum(eigP$values)/sum(eigP$values);
    K      = which(cumEig>0.9)[1];
  }
  # Get Eigenvectors
  vec = eigP$vectors[,1:K]; # take first K eigenvectors
  # Fourier coefficients
  ev1 = t(vec) %*% cov1 %*% vec/n1;
  ev2 = t(vec) %*% cov2 %*% vec/n2;
  # Compute test statistic
  tstat1= sum( (log(diag(ev1))-log(diag(ev2)))^2,na.rm=TRUE )/2
  dprod1= sqrt(outer(diag(ev1),diag(ev1)));
  dprod2= sqrt(outer(diag(ev2),diag(ev2)));
  tstat2= ((
    log((dprod1+ev1)/(dprod1-ev1)) -
    log((dprod2+ev2)/(dprod2-ev2))
  )/2)^2;
  tstat2= sum( tstat2*upper.tri(tstat2),na.rm=TRUE );
  tstat = n1*n2/N*(tstat1+tstat2);
  dof   = K*(K+1)/2;
  return( pchisq(tstat,dof,lower.tail=FALSE) );
}

