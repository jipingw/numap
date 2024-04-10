\name{NCP_wrapper4}
\alias{NCP_wrapper4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A wrapper function to compute the pipeline including parseing bam, 4-template deconvolution, unique and redundant map calling and occupancy calculation.
}
\description{
A wrapper function to compute the pipeline including parseing bam, 4-template deconvolution, unique and redundant map calling and occupancy calculation.}
\usage{
NCP_wrapper4(bam_file,paired=TRUE, result_path, chromosome_path,
                      temp4=NULL,center=TRUE,unique_map_wsize=120, 
wsize=73,parse_skip=FALSE,chrom_to_include =
                 NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file}{name of a single bam file}
  \item{paired}{paired or single-end}
  \item{result_path}{the folder where the results files are saved}
  \item{chromosome_path}{the folder where the chromosome files are}
  \item{temp4}{four template file}
  \item{center}{center or uniformly weighted occupancy to be calculate}
  \item{unique_map_wsize}{the window size to define the local maximum for the unique map}
  \item{wsize}{nucleosome occupancy window size for weights}
    \item{parse_skip}{skip parse bam file if TRUE}
\item{chrom_to_include}{a vector listing chromosome names of which the mapping data to output}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
