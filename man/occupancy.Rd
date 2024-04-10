\name{occupancy}
\alias{occupancy}
\title{
function to compute nucleosome occupancy}
\description{
function to compute nucleosome occupancy}

\usage{
occupancy(redundant_map_list,result_path, chrom_length=NULL,T4=FALSE,
center=TRUE,wsize=73)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{redundant_map_list}{a list of redundant maps by chromosomes}
  \item{result_path}{directory where the results to be saved}
  \item{chrom_length}{the length of each chromosomes. If not provided, the occupancy is calculated based on the maximum coordinate of the redundant map + wsize.}
  \item{T4}{T4=false to compute single-template, true for 4-template}
  \item{center}{center-weighted or uniformly weighted occupancy}
  \item{wsize}{range of weights, +/- wsize}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function outputs calculated occupancy per chromosome in files.
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
