# fdcov

This is a read-only mirror of the CRAN R package repository.  fdcov — Analysis of Covariance Operators  
https://cran.r-project.org/web/packages/fdcov/index.html

This package contains a collection of tools for performing statistical inference on functional data
specifically through an analysis of the covariance structure of the data. It includes two methods
for performing a k-sample test for equality of covariance in ksample.perm and ksample.com. For
supervised and unsupervised learning, it contains a method to classify functional data with respect to
each category’s covariance operator in classif.com, and it contains a method to cluster functional
data, cluster.com, again based on the covariance structure of the data.
The current version of this package assumes that all functional data is sampled on the same grid at
the same intervals. Future updates are planned to allow for the below methods to interface with the
fda package and its functional basis representations of the data.

Author(s): 
Alessandra Cabassi <ac2051@cam.ac.uk>, Adam B Kashlak <ak852@cam.ac.uk>

References:
Kashlak, Adam B, John AD Aston, and Richard Nickl (2016). "Inference on covariance operators
via concentration inequalities: k-sample tests, classification, and clustering via Rademacher
complexities", April, 2016 (in review)
Pigoli, Davide, John AD Aston, Ian L Dryden, and Piercesare Secchi. "Distances and inference for
covariance operators." Biometrika (2014): asu008.
