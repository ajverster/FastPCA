
This library contains functions for doing PCA quickly using an approximation.

FastPCA() behaves like prcomp() except that instead of using svd() it uses FastSVD()
FastSVD() does an approximation by randomly projecting your data into a subspace and then doing SVD on the orthonormal subspace of this projection. It is based on Halko et al. (http://arxiv.org/pdf/0909.4061v2.pdf). The iterative option can be used to reduce rounding errors


Testing/Basic.R shows the basic usage.

Testing/Example.R contains some example code testing the library out on the classic didgit recognition dataset, which is used in the Kaggle training competition.
Testing/train.csv Is the didgit classification dataset.
