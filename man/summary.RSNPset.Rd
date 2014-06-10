\name{summary.RSNPset}
\alias{summary.RSNPset}
\title{RSNPset Analysis Summary Function}
\description{
    Summary function to display execution information from \code{rsnpset()}.
}
\usage{
\method{summary}{RSNPset}(object, verbose = TRUE, ...) 
}

\arguments{
	\item{object}{Result from \code{rsnpset()}, an "RSNPset" S3 class object. Required.}
	\item{verbose}{Boolean inidcating if additional information about the resuts should be reported. Default is \code{TRUE}.}
	\item{...}{Additional arguments affecting the summary produced.}
}

\details{
If \code{verbose = TRUE}, prints a summary of the execution conditions of \code{rsnpset()}. The default report includes:
\itemize{
  \item The number of samples (the length of the \code{rsnpset()} argument \code{Y}).
  \item The range in sizes of the analyzed SNP sets.
  \item The number of SNP sets not analyzed (due to not containing any valid SNPs or other issues).
  \item The number of SNP sets containing SNPs not used in the analysis (due to being missing from the data or other issues).
}
If permutation replicates were generated, the report will also include:
\itemize{
  \item The number of permutation replicates computed (value of the \code{rsnpset()} argument \code{B}).
  \item Whether or not the ranks of the variance matrices of the permutation replicates are included in the results (value of the \code{rsnpset()} argument \code{ret.rank}).
  \item Whether or not the variance matrices were recomputed for each permutation replicate (value of the \code{rsnpset()} argument \code{v.permute}).
}
If \code{rsnpset()} was run with \code{pinv.check = TRUE}, the value of the argument \code{pinv.tol} will also be reported.
}

\value{
    If \code{rsnpset()} was run with \code{pinv.check = TRUE}, a list of \code{data.frame} objects is returned, each containing the following diagnositc measures of the calculated Penrose-Moore inverses for the observed and permutation data.
\tabular{ll}{
\bold{Column} \tab \bold{Absolute largest element of:} \cr
\code{d0} \tab \eqn{\Sigma - \bm{QDQ}}{\Sigma - QDQ} \cr
\code{d1} \tab \eqn{\Sigma \Sigma^{+} \Sigma - \Sigma}{V\%*\%\Sigma\%*\%V-V} \cr
\code{d2} \tab \eqn{\Sigma^{+} \Sigma \Sigma^{+} - \Sigma^{+}}{\Sigma\%*\%V\%*\%\Sigma-\Sigma} \cr
\code{d2} \tab \eqn{(\Sigma \Sigma^{+})^T - \Sigma \Sigma^{+}}{t(V\%*\%\Sigma)-V\%*\%\Sigma} \cr
\code{d4} \tab \eqn{(\Sigma^{+} \Sigma)^T - \Sigma^{+} \Sigma }{t(\Sigma\%*\%V)-\Sigma\%*\%V} \cr
}
where \eqn{\bm{QDQ}}{QDQ} is the spectral decomposition of \eqn{\Sigma}{V}. Departure of these values from zero indicates poor performance of the Penrose-Moore inverse. If \code{rsnpset()} was run with \code{pinv.check = FALSE}, function returns \code{NA}.
}

\note{
If \code{pinv.check = TRUE} and the number of permutations is large, the user may wish to capture the resulting diagnositc measures in an object for examination, as in the example below, as opposed to having them output directly.
}

\seealso{
    The function \code{\link{rsnpset}} provides a discription of the meaning of these reported values, as well as an explanation as to how they influence the results.
}

\examples{
n <- 200    # Number of patients
m <- 1000   # Number of SNPs

set.seed(123)
G <- matrix(rnorm(n*m), n, m)   # Normalized SNP expression levels
rsids <- paste0("rs", 1:m)      # SNP rsIDs 
colnames(G) <- rsids
 
K <- 15 # Number of SNP sets
genes <- paste0("XYZ", 1:K)     # Gene names 
gsets <- lapply(sample(3:50, size=K, replace=TRUE), sample, x=rsids)
names(gsets) <- genes

# Survival outcome
time <- rexp(n, 1/10)           # Survival time
event <- rbinom(n, 1, 0.9)      # Event indicator

\dontrun{
# Optional parallel backend
library(doParallel)
registerDoParallel(cores=8) }

# B >= 1000 is typically recommended
res <- rsnpset(Y=time, delta=event, G=G, snp.sets=gsets, 
               score="cox", B=50, ret.rank=TRUE, pinv.check=TRUE)

pinvcheck <- summary(res)
pinvcheck[["Observed"]]
}
