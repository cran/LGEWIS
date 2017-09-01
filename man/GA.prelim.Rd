\name{GA.prelim}
\alias{GA.prelim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The preliminary data management for GA (tests for genetic association)
}
\description{
Before testing a specific region using a generalized score type test, this function does the preliminary data management, such as fitting the model under the null hypothesis.
}
\usage{
GA.prelim(Y,time,X=NULL,corstr="exchangeable")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The outcome variable, an n*1 matrix where n is the total number of observations
}
  \item{time}{
An n*2 matrix describing how the observations are measured. The first column is the subject id. The second column is the measured exam (1,2,3,etc.).
}
  \item{X}{
An n*p covariates matrix where p is the total number of covariates.
}
  \item{corstr}{
The working correlation as specified in 'geeglm'. The following are permitted: "independence", "exchangeable", "ar1".
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  It returns a list used for function GA.test().
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
# time: describe longitudinal structure, n by 2 matrix
# G: genotype matrix, m by q matrix where m is the total number of subjects

data(LGEWIS.example)
Y<-LGEWIS.example$Y
X<-LGEWIS.example$X
time<-LGEWIS.example$time
G<-LGEWIS.example$G

# Preliminary data management
result.prelim<-GA.prelim(Y,time,X=X)
}

\keyword{preliminary work}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
