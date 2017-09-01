\name{GA.SSD.OneSet_SetIndex}
\alias{GA.SSD.OneSet_SetIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic association tests for a single region/gene using SSD format files
}
\description{
Test the genetic association between an quantitative outcome and one region/gene using SSD format files.
}
\usage{
GA.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, result.prelim, Gsub.id=NULL,
B=5000, ...)
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
Output of function "GA.prelim()".
}
  \item{Gsub.id}{
The subject id corresponding to the genotype matrix, an m dimensional vector. This is in order to match the phenotype and genotype matrix. The default is NULL, where the order is assumed to be matched with Y, X, E and time.
}
  \item{B}{
Number of Bootstrap replicates. The default is 5000.
}
  \item{...}{
Other options of the generalized score type test. Defined same as in function "GA.test()".
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
  \item{n.marker}{
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
