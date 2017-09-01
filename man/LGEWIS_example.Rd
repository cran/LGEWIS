\name{LGEWIS.example}
\alias{LGEWIS.example}
\docType{data}
\title{
Data example for LGEWIS (tests for genetic association or gene-environment interaction)
}
\description{
Example data for LGEWIS.
}
\format{
  LGEWIS.example contains the following objects:
  \describe{
    \item{G}{a numeric genotype matrix of 10 individuals and 20 SNPs.
	Each row represents a different individual, and each column represents a different SNP marker.}
    \item{Y}{a numeric vector of continuous phenotypes of 10 individuals with 4 repeated measurements.}
    \item{time}{a numeric matrix. The first column is the subject ID and the second column is the measured exam.}
    \item{X}{a numeric matrix of 1 covariates.}
    \item{E}{a numeric vector of environmental exposure.}
  }
}
\usage{data(LGEWIS.example)}
%\format{
%}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
%\examples{
%data(GEI.example)
%## maybe str(LGRF.example) ; plot(LGRF.example) ...
%}
\keyword{datasets}
