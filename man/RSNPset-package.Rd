\name{RSNPset-package}
\alias{RSNPset-package}
\alias{RSNPset}
\docType{package}
\title{The RSNPset Package}
\description{
    RSNPset is a software implementation of the use of efficient score statistics in genome-wide SNP set analysis with complex traits. 
}
\details{
\tabular{ll}{
Package: \tab RSNPset\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2014-06-10\cr
License: \tab GPL-3\cr
}

This package is designed for the analysis of sets of related SNPs, using genes as the loci of interest, but the methodology can naturally be applied to other genomic loci, including bands and pathways. The core function, \code{rsnpset()}, provides options for three efficient score statistics, Binomial, Gaussian, and Cox, for the analysis of binary, quantitative, and time-to-event outcomes, but is easily extensible to include others. Code implementing the inferential procedure is primarily written in \emph{C++} and utilizes parallelization to reduce runtime. A supporting function, \code{rsnpset.pvalue()}, offers easy computation of observed, permutation, FWER-adjusted, and FDR-adjusted p-values, and summary functions provide diagnostic measures and results metadata.
}
\author{
Chanhee Yi, Alexander Sibley, and Kouros Owzar

Maintainer: Chanhee Yi <chanhee.yi@dm.duke.edu>
}
\note{
The \emph{C++} code implementing the inferential procedures utilizes linear algebra routines from the \pkg{Eigen} library. This implementation is facilitated using the templates provided by the \pkg{Rcpp} and \pkg{RcppEigen} extension packages. Parallelization and reproducibility of the analysis is enabled by using the \pkg{doRNG} package in concert with the \pkg{foreach} package. The FDR-adjusted p-values are obtained using the \pkg{qvalue} package.
}
\seealso{
    Functions available in this package: \code{\link{rsnpset}}, \code{\link{rsnpset.pvalue}}, and supporting summary functions \code{\link{summary.RSNPset}}, \code{\link{summary.RSNPset.pvalue}}
    
    For more information on supporting packages, see: \code{Eigen}, \code{Rcpp}, \code{RcppEigen}, \code{\link{doRNG}}, \code{\link{foreach}}, \code{\link{qvalue}}
    
    The \code{snplist} package can be used to generate sets of SNPs for analysis with this package.
}
\examples{
n <- 200    # Number of patients
m <- 1000   # Number of SNPs

set.seed(123)

G <- matrix(rnorm(n*m), n, m)   # Normalized SNP expression levels
rsids <- paste0("rs", 1:m)      # SNP rsIDs 
colnames(G) <- rsids
 
K <- 10 # Number of SNP sets
genes <- paste0("XYZ", 1:K)     # Gene names 
gsets <- lapply(sample(3:50, size=K, replace=TRUE), sample, x=rsids)
names(gsets) <- genes

# Survival outcome
time <- rexp(n, 1/10)           # Survival time
event <- rbinom(n, 1, 0.9)      # Event indicator

res <- rsnpset(Y=time, delta=event, G=G, snp.sets=gsets, score="cox")
head(res)
summary(res)
rsnpset.pvalue(res)

\dontrun{
# Optional parallel backend
library(doParallel)
registerDoParallel(cores=8)

res <- rsnpset(Y=time, delta=event, G=G, snp.sets=gsets, 
               score="cox", B=1000)
rsnpset.pvalue(res) }

# Binary outcome
set.seed(123)
Y <- rbinom(n, 1, 0.5) 

head(rsnpset(Y=Y, G=G, snp.sets=gsets, score="binomial", v.method="empirical"))
head(rsnpset(Y=Y, G=G, snp.sets=gsets, score="binomial", v.method="asymptotic"))
}
\keyword{ 
htest
}

