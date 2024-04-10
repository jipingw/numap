\name{trainTEMP1}
\alias{trainTEMP1}
\title{
This function trains the single-template model based on cleavages.
}
\description{
This function trains the single-template model based on cleavages.
}
\usage{
trainTEMP1(w_files,c_files,result_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result_path}{directory where the results to be saved}
  \item{w_files}{a list of raw cleavage files on Watson strands by chromosomes}
  \item{c_files}{a list of raw cleavage files on Crick strands by chromosomes}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Tt returns a vector of length 8 to define the relative cleavage frequency at 8 positions including -2,-1,0,1, 4,5,6,7.
The single template is needed for calling deconvolution algorithm NCP1.
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
