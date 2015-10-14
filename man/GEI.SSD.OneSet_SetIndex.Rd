\name{GEI.SSD.OneSet_SetIndex}
\alias{GEI.SSD.OneSet_SetIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene-environment interaction tests for a single region/gene using SSD format files
}
\description{
Test the interaction between an environmental exposure and one region/gene on a quantitative outcome using SSD format files.
}
\usage{
GEI.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, result.prelim, Gsub.id=NULL,
MinP.adjust=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SSD.INFO}{
SSD format information file, output of function ``Open_SSD". The sets are defined by this file.
}
  \item{SetIndex}{
Set index. From 1 to the total number of sets.
}
  \item{result.prelim}{
Output of function "GEI.prelim()".
}
  \item{Gsub.id}{
The subject id corresponding to the genotype matrix, an m dimensional vector. This is in order to match the phenotype and genotype matrix. The default is NULL, where the order is assumed to be matched with Y, X, E and time.
}
  \item{MinP.adjust}{
If the users would like to compare with the MinP test, this parameter specify the adjustment thereshold as in Gao, et al. (2008) "A multiple testing correction method for genetic association studies using correlated single nucleotide polymorphisms". Values from 0 to 1 are permitted. The default is NULL, i.e., no comparison. The value suggested by Gao, et al. (2008) is 0.95.
}
  \item{...}{
Other options of the generalized score type test. Defined same as in function "GEI.test()".
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
number of tested SNPs in the SNP set.
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

# * Since the Plink data files used here are hard to be included in a R package,
# The usage is marked by "#" to pass the package check.

# library(LGEWIS)

##############################################

# Plink data files: File.Bed, File.Bim, File.Fam
# Files defining the sets: File.SetID, File.SSD, File.Info
# For longitudinal data, outcome and covariates are saved in a separate file: Y, time, X, E.
# Preliminary work was done using function null.LGRF.

# Create the MW File
# File.Bed<-"./example.bed"
# File.Bim<-"./example.bim"
# File.Fam<-"./example.fam"
# File.SetID<-"./example.SetID"
# File.SSD<-"./example.SSD"
# File.Info<-"./example.SSD.info"

# Generate SSD file
# To use binary ped files, you have to generate SSD file first.
# If you already have a SSD file, you do not need to call this function.
# Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

# SSD.INFO<-Open_SSD(File.SSD, File.Info)
# Number of samples
# SSD.INFO$nSample
# Number of Sets
# SSD.INFO$nSets

## Fit the null model
# Y: outcomes, n by 1 matrix where n is the total number of observations
# X: covariates, n by p matrix
# E: covariates, n by 1 matrix
# time: describe longitudinal structure, n by 2 matrix
# result.prelim<-GEI.prelim(Y,time,E,X=X)

## Test a single region
# out_single<-GEI.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO,
# SetIndex=3, result.prelim=result.prelim,MinP.adjust=0.95)

# Example result
# out_single
# $n.marker
# [1] 108

# $E.df
# [1] 22

# $G.df
# [1] 22

# $p.single
#           MAF      p.value
# SNP0254 0.098 5.920758e-01
# SNP0273 0.091 9.468959e-01
# SNP0199 0.118 4.131540e-01
# ...

# $p.MinP
# [1] 3.876731e-05

# $p.value
# [1] 0.1167991

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plink_test_single}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
