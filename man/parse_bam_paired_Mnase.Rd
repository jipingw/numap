\name{parse_bam_paired_Mnase}
\alias{parse_bam_paired_Mnase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to parse bam file of paired MNase map.
}
\description{
This function parses bam file of paired MNase map.
}
\usage{
parse_bam_paired_Mnase(bam_file, result_path,chrom_to_include)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file}{
the bam file from paired end alignment
}
  \item{result_path}{
directory path to save the parsed bam file
}
\item{chrom_to_include}{a vector listing chromosome names of which the mapping data to be 
processed}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
output parsed bam file, each row contains start and end positions.
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
function (bam_file, result_path, insert_length_low, insert_length_high)
{
    setwd(result_path)
    bf = BamFile(bam_file, asMates = TRUE, qnameSuffixStart = ".")
    gr = as(seqinfo(bf), "GRanges")
    chrom_length = end(gr)
    chrom = as.vector(unlist(runValue(seqnames(gr))))
    if (!file.exists(paste(bam_file, ".bai", sep = ""))) {
        cat("Index", bam_file, "\n")
        indexBam(bam_file)
    }
    what = c("rname", "pos", "mpos", "isize")
    for (i in 1:length(chrom)) {
        cat("\t", chrom[i])
        ch = chrom[i]
        if (grepl("M", ch) == FALSE) {
            which <- gr[i]
            param = ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                isSecondaryAlignment = FALSE), what = what, mapqFilter = 255,
                which = which)
            bam <- scanBam(bam_file, param = param)
            bamout = data.table(chr = bam[[1]]$rname, pos1 = bam[[1]]$pos,
                mate_pos = bam[[1]]$mpos, frag_size = bam[[1]]$isize)
            rm(bam)
            bamout1 = bamout[frag_size >= insert_length_low &
                frag_size <= insert_length_high, ]
            rm(bamout)
            bamout2 = data.table(chr = bamout1$chr, start = bamout1$pos1,
                end = bamout1$pos1 + bamout1$frag_size - 1)
            rm(bamout1)
            fwrite(bamout2, file = paste(ch, ".reads.txt", sep = ""),
                sep = "\t")
        }
    }
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
