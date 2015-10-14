\name{GA.SSD.All}
\alias{GA.SSD.All}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic association tests for multiple regions/genes using SSD format files
}
\description{
Test the association between an quantitative outcome and multiple region/genes using SSD format files.
}
\usage{
GA.SSD.All(SSD.INFO, result.prelim, Gsub.id=NULL, MinP.adjust=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SSD.INFO}{
SSD format information file, output of function ``Open_SSD". The sets are defined by this file.
}
  \item{result.prelim}{
Output of function "GA.prelim()".
}
  \item{Gsub.id}{
The subject id corresponding to the genotype matrix, an m dimensional vector. This is in order to match the phenotype and genotype matrix. The default is NULL, where the order is assumed to be matched with Y, X, E and time.
}
  \item{MinP.adjust}{
If the users would like to compare with the MinP test, this parameter specify the adjustment thereshold as in Gao, et al. (2008) "A multiple testing correction method for genetic association studies using correlated single nucleotide polymorphisms". Values from 0 to 1 are permitted. The default is NULL, i.e., no comparison. The value suggested by Gao, et al. (2008) is 0.95.
}
  \item{...}{
Other options of the generalized score type test. Defined same as in function "GA.test()".
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{results}{
Results of the set based analysis. First column contains the set ID; Second column (second and third columns when the MinP test is compared) contains the p-values; Last column contains the number of tested SNPs.
}
  \item{results.single}{
Results of the single variant analysis for all variants in the sets. First column contains the regions' names; Second column is the variants' names; Third column contains the minor allele frequencies; Last column contains the p.values.
}
}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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

#  library(LGEWIS)

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
# time: describe longitudinal structure, n by 2 matrix
# result.prelim<-GA.prelim(Y,time,X=X)

## Test all regions
# out_all<-GA.SSD.All(SSD.INFO, result.prelim, MinP.adjust=0.95)

# Example result
# out.all$results
#      SetID    P.value P.value_MinP N.Marker
# 1  GENE_01 0.91592151   1.00000000       94
# 2  GENE_02 0.31103681   0.06451609       84
# 3  GENE_03 0.05976685   0.49923701      108
# 4  GENE_04 0.09389408   1.00000000      101
# 5  GENE_05 0.10339403   0.67139314      103
# 6  GENE_06 0.94666614   1.00000000       94
# 7  GENE_07 0.47955756   1.00000000      104
# 8  GENE_08 0.99415148   1.00000000       96
# 9  GENE_09 0.29271727   0.78406742      100
# 10 GENE_10 0.29061419   0.67045802      100

# out.all$results.single
#     Region.name SNP.name                MAF              p.value
# 1       GENE_01  SNP0056              0.097    0.343641626676761
# 2       GENE_01  SNP0083               0.11     0.75884859949668
# 3       GENE_01  SNP0035              0.097    0.796572155483814
# ...


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plink_test_all}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
