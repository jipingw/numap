\name{TSS_occu}
\alias{TSS_occu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate occupancy around TSS
}
\description{
This funciton calculates occupancy around TSS.}
\usage{
TSS_occu(TSS_file, occu_list, wsize=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TSS_file}{TSS_file must tab or space_delimited file and contain three columns, chr, strand and start (start is TSS start), 1 for watson and 2 for crick strand. chr must match the chromosome names in Occupancy file.
 }
  \item{occu_list}{occupancy file outputted by `occupancy` function}
  \item{wsize}{window size. default is +/-100 bp.}
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
