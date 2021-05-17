\name{network}
\alias{network}
\docType{data}
\title{
dataframe of interacting proteins of C.elegans from STRING database
}
\description{
The first two columns corresponds to the interacting proteins and other columns represent scores on different qualifiers for the pair of interaction and the combined score.
}
\usage{network}
\format{
  A data frame with 405668 observations on the following 11 variables.
  \describe{
    \item{\code{protein1}}{a character vector-- transcript id of protein1}
    \item{\code{protein2}}{a character vector-- transcript id of protein2}
    \item{\code{neighborhood}}{a numeric vector}
    \item{\code{fusion}}{a numeric vector}
    \item{\code{cooccurence}}{a numeric vector}
    \item{\code{coexpression}}{a numeric vector}
    \item{\code{experimental}}{a numeric vector}
    \item{\code{database}}{a numeric vector}
    \item{\code{textmining}}{a numeric vector}
    \item{\code{combined_score}}{a numeric vector}
    \item{\code{combined_score_weight}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
STRING DB
}
\references{\url{http://string-db.org/}}
\examples{
}
\keyword{datasets}
