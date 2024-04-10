\name{plotAATT}
\alias{plotAATT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot AATT frequency based on unique map
}
\description{
Function to plot AA/TT/TA/AT fequency based on a nucleosome map}
\usage{
plotAATT(chromosome_path,center_list,result_path,NCP_thresh=0, wsize=73)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chromosome_path}{path to chromosome fasta files}
  \item{result_path}{directory where the results to be saved}
  \item{center_list}{list of nucleosome maps by chromosomes}
  \item{NCP_thresh}{NCP_thresh to define which nucleosome to be used from the unique map}
  \item{wsize}{window size to plot the AATT frequency}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It returns AA/TT/TA/AT frequency plot in-screen around nucleosome centers.
The function also returns a vector of AA/TT/TA/AT frequency in the defined region.
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
