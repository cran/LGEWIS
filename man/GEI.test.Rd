\name{GEI.test}
\alias{GEI.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the interaction between an environemental exposure and a region/gene by a generalized score type test.
}
\description{
Once the preliminary work is done using "GEI.prelim()", this function tests a specifc region/gene. Single SNP analyses are also incorporated.
}
\usage{
GEI.test(result.prelim,G,Gsub.id=NULL,G.method='wPCA',G.df=floor(sqrt(nrow(G))),
bootstrap=NULL,MinP.adjust=NULL,impute.method='fixed')
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
  \item{G.method}{
The dimension reduction method for main effect adjustment of G. The following are permitted: "wPCA" for weighted principal component analysis; "PCA" for principal component analysis; "R2" for ordering the principal components by their R-squares. The dimension reduction method is in order to analyze large regions, i.e., the number of variants is close to or larger than the number of subjects. The default is "wPCA".
}
  \item{G.df}{
Number of components selected by the dimension reduction method. The default is sqrt(m).
}
  \item{bootstrap}{
Whether to use bootstrap for small sample size adjustement. This is recommended when the number of subjects is small, or the set contains rare variants. The default is NULL, but a suggested number is 10000 when it is needed.
}
  \item{MinP.adjust}{
If the users would like to compare with the MinP test, this parameter specify the adjustment thereshold as in Gao, et al. (2008) "A multiple testing correction method for genetic association studies using correlated single nucleotide polymorphisms". Values from 0 to 1 are permitted. The default is NULL, i.e., no comparison. The value suggested by Gao, et al. (2008) is 0.95.
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
  \item{p.MinP}{
P-value of the MinP test.
}
  \item{n.marker}{
number of heterozygous SNPs in the SNP set.
}
  \item{E.df}{
number of tested SNPs in the SNP set.
}
  \item{G.df}{
number of tested SNPs in the SNP set.
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
## GEI.prelim does the preliminary data management.
# Input: Y, time, E, X (covariates)
## GEI.test tests a region.
# Input: G (genetic variants) and result of GEI.prelim

library(LGEWIS)

# Load data example
# Y: outcomes, n by 1 matrix where n is the total number of observations
# X: covariates, n by p matrix
# E: environmental exposure, n by 1 matrix
# time: describe longitudinal structure, n by 2 matrix
# G: genotype matrix, m by q matrix where m is the total number of subjects

data(LGEWIS.example)
Y<-LGEWIS.example$Y
X<-LGEWIS.example$X
E<-LGEWIS.example$E
time<-LGEWIS.example$time
G<-LGEWIS.example$G

# Preliminary data management
result.prelim<-GEI.prelim(Y,time,E,X=X)

# test without the MinP test
result<-GEI.test(result.prelim,G,MinP.adjust=NULL)

# test with the MinP test
result<-GEI.test(result.prelim,G,MinP.adjust=0.95)

# test with the MinP test and the small sample adjustment
result<-GEI.test(result.prelim,G,MinP.adjust=0.95,bootstrap=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{test}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
