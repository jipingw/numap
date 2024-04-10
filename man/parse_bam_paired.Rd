\name{parse_bam_paired}
\alias{parse_bam_paired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to parse paired-end chemical mapping bam file}
\description{
function to parse paired-end chemical mapping bam file}
\usage{
parse_bam_paired(bam_file,result_path,chrom_to_include)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file}{bam file from paired-end reads}
  \item{result_path}{directory where results are stored}
\item{chrom_to_include}{a vector listing the chromosome names of which the data to output}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
output cleavage frequency at each genomic location on watson and crick strands separately for each chromosome.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
