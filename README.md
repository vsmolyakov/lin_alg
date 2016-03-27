# lin alg
Numerical Linear Algebra Algorithms.

### Description

**Non-negative Matrix Factorization**

Non-negative Matrix Factorization (NMF) decomposes matrix A into non-negative factors W and H that minimize the Frobenius norm ||A-WH||. The non-negativity of W and H is useful in application to probability distrubitions. In particular, we apply NMF to topic modeling by decomposing a term-document matrix of word counts of size V x D into a V x K topic matrix W and a K x D topic proportions matrix H.

<p align="center">
<img src="https://github.com/vsmolyakov/lin_alg/blob/master/figures/nmf_merged.png" />
</p>

The figure above shows the term-document matrix (left) decomposed into the topic matrix (middle) and the topic proportions matrix (right) for K = 4 topics. We can see the topic distributions in the middle where each word is assigned to one of K topics. The plot on the right shows the topic proportions in each document and can be used to measure document similarity for document matching and retrieval.

References:  
*D. Lee and S. Seung, "Algorithms for Non-negative Matrix Factorization", NIPS 2000*  

**Rayleigh Quotient Iteration**

Rayleigh Quotient Iteration is an eigenvalue estimation algorithm that iteratively refines eigenvalue estimates via Rayleigh quotient and eigenvector estimates via inverse iteration.

<p align="center">
<img src="https://github.com/vsmolyakov/lin_alg/blob/master/figures/rayleigh_quotient.png" width = "400" />
</p>

Figure above shows fast convergence of the algorithm to the maximum eigenvalue estimates for a Gaussian random matrix.

References:  
*G. Golub and C. Van Loan, "Matrix Computations", 1983*  


**Modified Gram Schmidt**

The Modified Gram Schmidt (MGS) algorithm is used to find QR matrix decomposition. The figure below compares numerical stability of the classical Gram-Schmidt algorithm with MGS.

<p align="center">
<img src="https://github.com/vsmolyakov/lin_alg/blob/master/figures/gram_schmidt.png" width = "400" />
</p>

We can see from the plot that the MGS algorithm is able to maintain higher precision accuracy for the diagonal entries of R compared to the classic algorithm.

References:  
*L. Trefethen and D. Bau III, "Numerical Linear Algebra", 1997*

**Legendre Polynomials**

Legendre polynomials are orthogonal basis polynomials and can be computed via QR decomposition of the Vandermonde matrix.

<p align="center">
<img src="https://github.com/vsmolyakov/lin_alg/blob/master/figures/legendre_poly.png" width = "400" />
</p>

The figure above shows the first four Legendre polynomials.

References:  
*R. Horn and C. Johnson, "Matrix Analysis", 2012*  

**Numerical Stability**

A stable algorithm gives nearly the right answer to nearly the right question. Here, numerical stability is examined for polynomials, householder triangularization, SVD decomposition, least squares, and matrix inversion.

<p align="center">
<img src="https://github.com/vsmolyakov/lin_alg/blob/master/figures/numerical_stability.png" width = "400" />
</p>

The figure above shows the sensitivity of polynomial approximation near its roots showing multiple zero crossings.

References:  
*G. Strang, "Introduction to Linear Algebra", 2009*
 
### Dependencies

Matlab 2014a
