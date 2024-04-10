\name{occupancy_Mnase_paired_par}
\alias{occupancy_Mnase_paired_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to calculate Mnase map occupancy in  parallel}
\description{
This function compute the nucleosome occupancy based on MNase maps in parallel.}
\usage{
occupancy_Mnase_paired_par(reads_file_list, result_path, chrom_length=NULL,center=TRUE,wsize=73, insert_length_low, insert_length_high,cores=NULL, maxhits=Inf,nameExt=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads_file_list}{list of parsed bam file from `parse_bam_paired_Mnase function}
  \item{result_path}{directory to save the results}
  \item{chrom_length}{vector to contain the chromosome length. if `NULL`, the maximum end position + wsize will be used}
  \item{center}{whether to calcualte center weighted or uniformly weighted occupancy}
  \item{wsize}{the widown size to calcualte the center weighted occupany. Ignored if `center=FALSE`}
    \item{insert_length_low}{
the lower limit of read length to keep for calculation of occupancy
}
  \item{insert_length_high}{
the upper limit of read length to keep for calculation of occupancy
}
  \item{cores}{the number of cores to use in parallel}
  \item{maxhits}{the maximum of hits that a read can have to be included in occupancy calculation}
  \item{nameExt}{the extension to be added to the occupancy output file name}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Output occupancy value by chromosomes.
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
