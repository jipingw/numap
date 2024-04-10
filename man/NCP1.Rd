\name{NCP1}
\alias{NCP1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to perform deconvolution algorithm using single template.
}
\description{
Main function to perform deconvolution algorithm using 4 templates.
}
\usage{
NCP1(result_path, chromosome_path, temp1=NULL,unique_map_wsize=120)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result_path}{directory where the results are stored}
  \item{chromosome_path}{path to where the chromosomes are stored, one file for each chromosome}
  \item{temp1}{file contains the template}
  \item{unique_map_wsize}{distance to define the local maximum for unique map}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NCP1 generates the NCP scores and the unique and redudnant maps
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
