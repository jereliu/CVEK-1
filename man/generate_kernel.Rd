% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_kernel.R
\name{generate_kernel}
\alias{generate_kernel}
\title{Generating A Single Kernel}
\usage{
generate_kernel(method = "rbf", l = 1, d = 2)
}
\arguments{
\item{method}{(character) A character string indicating which kernel 
is to be computed.}

\item{l}{(numeric) A numeric number indicating the hyperparameter 
(flexibility) of a specific kernel.}

\item{d}{(integer) For polynomial, d is the power; for matern, v = d + 1 / 2; for
rational, alpha = d.}
}
\value{
\item{kern}{(function) A function indicating the generated kernel.}
}
\description{
Generate kernels for the kernel library.
}
\details{
There are seven kinds of kernel available here. For convenience, we define
\eqn{r=\mid x-x'\mid}.

\bold{Gaussian RBF Kernels} \deqn{k_{SE}(r)=exp\Big(-\frac{r^2}{2l^2}\Big)}

\bold{Matern Kernels}
\deqn{k_{Matern}(r)=\frac{2^{1-\nu}}{\Gamma(\nu)}\Big(\frac{\sqrt{2\nu
r}}{l}\Big)^\nu K_\nu \Big(\frac{\sqrt{2\nu r}}{l}\Big)}

\bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
l^2}\Big)^{-\alpha}}

\bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^d} We have intercept
kernel when \eqn{d=0}, and linear kernel when \eqn{d=1}.

\bold{Neural Network Kernels} \deqn{k_{NN}(x,
x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T
\tilde{x}'}{\sqrt{(1+2\tilde{x}^T \tilde{x})(1+2\tilde{x}'^T
\tilde{x}')}}\Big)}
}
\examples{


kern_func_list <- list()
for (j in 1:nrow(kern_par)) {
  kern_func_list[[j]] <- generate_kernel(kern_par[j, ]$method,
                                         kern_par[j, ]$l,
                                         kern_par[j, ]$d)
}


}
\references{
The MIT Press. Gaussian Processes for Machine Learning, 2006.
}
\author{
Wenying Deng
}
