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
\details{
Please see SKAT vignettes for using SSD files.
%%  ~~ If necessary, more details than the description above ~~
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plink_test_single}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
