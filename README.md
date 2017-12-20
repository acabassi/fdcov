# fdcov

R package for the Analysis of Covariance Operators.
v 1.0 is available on CRAN:  https://cran.r-project.org/web/packages/fdcov/index.html

This package contains a collection of tools for performing statistical inference on functional data specifically through an analysis of the covariance structure of the data. It includes two methods for performing a k-sample test for equality of covariance in ksample.perm and ksample.com. For supervised and unsupervised learning, it contains a method to classify functional data with respect to each category's covariance operator in classif.com, and it contains a method to cluster functional data, cluster.com, again based on the covariance structure of the data.
The current version of this package assumes that all functional data is sampled on the same grid at the same intervals. Future updates are planned to allow for the below methods to interface with the fda package and its functional basis representations of the data.

Authors:

Alessandra Cabassi <ac2051@cam.ac.uk>, Adam B Kashlak <kashlak@ualberta.ca>

Contributors:

Davide Pigoli <davide.pigoli@kcl.ac.uk>

References:

Cabassi, A., Pigoli, D., Secchi, P. and Carter, P.A., 2017. Permutation tests for the equality of covariance operators of functional data with applications to evolutionary biology. Electronic Journal of Statistics, 11(2), pp.3815-3840.

Kashlak, A.B., Aston, J.A. and Nickl, R., 2016. Inference on covariance operators via concentration inequalities: k-sample tests, classification, and clustering via Rademacher complexities. arXiv preprint arXiv:1604.06310.

Pigoli, D., Aston, J.A., Dryden, I.L. and Secchi, P., 2014. Distances and inference for covariance operators. Biometrika, 101(2), pp.409-422.
