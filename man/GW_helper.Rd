\name{GW_helper}
\alias{GW_helper}
\title{
Function to compute weighted change statistics for geometrically weighted terms.}
\description{
Computes geometric weights for a pre-specified subgraph count given a known decay parameter. Intended primarily for use with \code{ergm.MSE} to simplify specifications for values of lower-order terms.}
\usage{
GW_helper(term_count,
            decay,
            lower_bound=0)
}
\arguments{
  \item{term_count}{is the count of subgraphs to which a researcher wants to assign geometric weights. Should be provided as a scalar (single numeric value).
  }
  \item{decay}{is the decay parameter used to assign geometric weights provided by the researcher.

  }
  \item{lower_bound}{is the baseline value for the subgraph count when calculating the geometric weights. Default value is 0.

  }
}
\details{
Computes the weighted value to assign to change statistics for a desired subgraph count. \code{term_count} corresponds to the subgraph count to assign geometric weights. For example, consider an ERGM specification where triadic effects are specified with a \code{gwesp} statistic. \code{GW_helper} would be used to determine the change statistic for a given triad count. If the triad count is 2, this would be specified by setting the \code{term_count} argument to equal 2 and providing the decay parameter used to estimate the \code{gwesp} term in an \code{ergm} call. \code{GW_helper} would return the change statistic for the \code{gwesp} term when a potential tie increases the triad count to a value of 2, weighted by \code{decay}.

The \code{lower_bound} argument assigns the baseline value for the subgraph count in question. By default, \code{lower_bound} is set to 0. This means that, in the above example, if \code{term_count} is set to 2, the default behavior of \code{GW_helper} is to calculate the weighted change statistic for a \code{gwesp} term that increases the number of triangles from 0 to 2. If \code{lower_bound} is set to 1 and \code{term_count} is set to 3, then \code{GW_helper} will calculate the change statistic when increasing the number of triangles from 1 to 3. Since geometrically weighted terms assign diminishing weights for each statistic added, selection of \code{lower_bound} and \code{term_count} will affect the value of the resulting change statistic in meaningful ways.

 If the geometrically weighted term is estimated with fixed decay, \code{decay} is the value provided by the researcher within the \code{ergm} formula. If the term is estimated without fixed decay, \code{decay} is the decay parameter estimate from the data corresponding to a single entry in the ergm summary table parameter vector.

}
\value{
Returns a scalar value reflecting the change statistic for a potential tie that increases the subgraph count captured by a geometrically weighted term by the values assigned to \code{lower_bound} and \code{term_count}. Intended primarily to facilitate use of \code{ergm.MSE} and \code{ergm.msma}.
}
\references{
Duxbury, Scott W. 2025. "Marginal Structural Effects for ERGM" \emph{Working paper.}

Hunter, David. 2007. "Curved Exponential Family Models for Social Networks" \emph{Social Networks.}
}
\author{
Scott Duxbury, Associate Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}


\examples{
\dontshow{
require(ergm)

GW_helper(1,.5)

}
\donttest{

library(ergm)
set.seed(21093)
data("faux.dixon.high")

# fit a model
my.ergm<-ergm(faux.dixon.high~edges+
                nodeicov("grade")+
                nodeocov("grade")+
                nodeifactor("sex")+
                nodeofactor("sex")+
                absdiff("grade")+
                nodematch("sex")+
                mutual+
                gwidegree(.5,fixed=TRUE))

#for a fixed decay of .5, what is the change statistic for a tie that increases
  # degree from 0 to 2?

GW_helper(2,.5)

#for a fixed decay of .5, what is the change statistic for a tie that increases
  # degree from 1 to 2?
GW_helper(2,.5,lower_bound=1)


##within a call to ergm.MSE
  #creating a single mutual tie requires, at a minimum, increasing indegree
    #from 0 to 1 for the recipient node

ergm.MSE(my.ergm,
          substructural_effect="mutual",
          lower_order_term="gwideg.fixed.0.5",
          at.lower_order_term=GW_helper(1,.5),
          estimate="MSEm")


#What happens if the tie moves indegree from 4 to 5 for the recipient node?

ergm.MSE(my.ergm,
          substructural_effect="mutual",
          lower_order_term="gwideg.fixed.0.5",
          at.lower_order_term=GW_helper(5,.5,lower_bound=4),
          estimate="MSEm")
  }

}

\keyword{ ~ergm }
\keyword{ ~networks }
\keyword{ ~marginal effects}
