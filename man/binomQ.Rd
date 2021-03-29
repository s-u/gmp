\name{binomQ}
\title{Exact Rational Binomial Probabilities}
\alias{dbinomQ}
%NOT YET \alias{pbinomQ}
\description{
  Compute \emph{exact} binomial probabilities using (big integer and) big
  rational arithmetic.
}
\usage{
dbinomQ(x, size, prob, log = FALSE)
}
%% pbinomQ(q, size, prob, lower.tail = TRUE, log.p = FALSE)
%
\arguments{
  \item{x, size}{integer or big integer (\code{"\link{bigz}"}), will be
    passed to \code{\link{chooseZ}()}.}
  \item{prob}{the probability; should be big rational
    (\code{"\link{bigq}"}); if not it is coerced with a warning.}
  \item{log}{logical; must be \code{FALSE} on purpose.}
}
\value{a big rational (\code{"\link{bigq}"}) of the \code{\link{length}} of
  (recycled) \code{x+size+prob}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{chooseZ}}; \R's (\pkg{stats} package) \code{\link{dbinom}()}.
}
%% \examples{
%% TODO
%% }
