\name{GA.test}
\alias{GA.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the association between an quantitative outcome variable and a region/gene by a generalized score type test.
}
\description{
Once the preliminary work is done using "GA.prelim()", this function tests a specifc region/gene. Single SNP analyses are also incorporated.
}
\usage{
GA.test(result.prelim,G,Gsub.id=NULL,weights='beta',B=5000, B.coef=NULL,
impute.method='fixed')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result.prelim}{
The output of function "GEI.prelim()"
}
  \item{G}{
Genetic variants in the target region/gene, an m*q matrix where m is the subject ID and q is the total number of genetic variables. Note that the number of rows in Z should be same as the number of subjects.
}
  \item{Gsub.id}{
The subject id corresponding to the genotype matrix, an m dimensional vector. This is in order to match the phenotype and genotype matrix. The default is NULL, where the order is assumed to be matched with Y, X and time.
}
  \item{weights}{
Can be a numeric vector of weights for genetic variants (The length should be same as the number of genetic variants in the set.), or pre-determined weights: "beta" (beta weights as in SKAT paper), "rare" (restrited to MAF<0.01), "common" (restrited to MAF>0.01). The default is NULL, where the flat weights are applied.
}
  \item{B}{
Number of Bootstrap replicates. The default is 5000.
}
  \item{B.coef}{
Direct import of Bootstrap coefficients, an m by B matrix. This is in order to efficiently implement the Bootstrap step. The default is NULL.
}
  \item{impute.method}{
Choose the imputation method when there is missing genotype. Can be "random", "fixed" or "bestguess". Given the estimated allele frequency, "random" simulates the genotype from binomial distribution; "fixed" uses the genotype expectation; "Best guess" uses the genotype with highest probability.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{p.value}{
P-value of the set based generalized score type test.
}
  \item{p.single}{
P-values of the incorporated single SNP analyses
}
  \item{n.marker}{
number of heterozygous SNPs in the SNP set.
}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## GA.prelim does the preliminary data management.
# Input: Y, time, X (covariates)
## GA.test tests a region.
# Input: G (genetic variants) and result of GEI.prelim

library(LGEWIS)

# Load data example
# Y: outcomes, n by 1 matrix where n is the total number of observations
# X: covariates, n by p matrix
# time: describe longitudinal structure, n by 2 matrix
# G: genotype matrix, m by q matrix where m is the total number of subjects

data(LGEWIS.example)
Y<-LGEWIS.example$Y
X<-LGEWIS.example$X
time<-LGEWIS.example$time
G<-LGEWIS.example$G

# Preliminary data management
result.prelim<-GA.prelim(Y,time,X=X)

# test with 5000 bootstrap replicates
result<-GA.test(result.prelim,G,B=5000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{test}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
