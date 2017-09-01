\name{GEI.prelim}
\alias{GEI.prelim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The preliminary data management for GEI (tests for gene-environment interaction)
}
\description{
Before testing a specific region using a generalized score type test, this function does the preliminary data management, such as pareparing spline basis functions for E etc..
}
\usage{
GEI.prelim(Y,time,E,X=NULL,E.method='ns',E.df=floor(sqrt(length(unique(time[,1])))),
corstr="exchangeable")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The outcome variable, an n*1 matrix where n is the total number of observations
}
  \item{time}{
An n*2 matrix describing how the observations are measured. The first column is the subject id. The second column is the measured exam (1,2,3,etc.).
}
  \item{E}{
An n*1 environmental exposure.
}
  \item{X}{
An n*p covariates matrix where p is the total number of covariates.
}
  \item{E.method}{
The method of sieves for the main effect of E. It can be "ns" for natural cubic spline sieves; "bs" for B-spline sieves; "ps" for polynomial sieves. The default is "ns".
}
  \item{E.df}{
Model complexity for the method of sieves, i.e., number of basis functions. The default is sqrt(m).
}
  \item{corstr}{
The working correlation as specified in 'geeglm'. The following are permitted: "independence", "exchangeable", "ar1".
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  It returns a list used for function GEI.test().
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
}

\keyword{preliminary work}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
