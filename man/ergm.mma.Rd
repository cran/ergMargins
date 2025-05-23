\name{ergm.mma}
\alias{ergm.mma}
\title{
Function to conduct marginal effects mediation analysis in ERGM.
}
\description{
ergm.mma is a function to conduct mediation analysis in ERGM using marginal effects. It computes the total, direct, and indirect average marginal effects for two nested ERGMs fit to the same data. Can be used with discrete or continuous mediators and for multiple mediator analyses. Currently accepts \code{ergm}, \code{btergm}, \code{mtergm}, \code{ergm.count}, \code{mlergm}, and \code{ergm.multi} objects.}
\usage{
ergm.mma(restricted.model,
          full.model,
          direct.effect,
          mediator,
          at.controls=NULL,
          control_vals=NULL,
          ME="AME")
}
\arguments{
  \item{restricted.model}{the ergm without the mediator variable.
}
  \item{full.model}{the ergm including the mediator variable.}

  \item{mediator}{a character string containing the name of the mediator variable. When examining joint mediation with multiple mediators, this is a vector with each entry containing the character strings for each mediator variable.}
  \item{direct.effect}{a character string containing the name of the direct effect.}
  \item{at.controls}{a vector of character strings containing the names of control variables to hold at prespecified values.}
  \item{control_vals}{a vector of numeric values to hold at.controls at during estimation.}
  \item{ME}{the type of marginal effect to calculate. Currently accepts either "AME" or "MEM". Defaul it "AME"}
}
\details{
ergm.mma conducts marginal effects mediation analysis for nested ERGMs fit to the same network by computing the change in average marginal effects between models. The AMEs are estimated on the scale of tie probabilities, not the linear component (log-odds). Standard errors are computed using the Delta method.

A few notes on using ergm.mma:
\itemize{

\item ergm.mma default specification uses the AME, which can be computationally intensive in large networks or high dimensional models. If run times are an issue, try specifying ME="MEM".

\item ergm.mma can be applied to discrete mediators or multiple mediator analysis.


}

}
\value{
Returns a table including the point estimates and standard errors for total, direct, and indirect effects. A single indirect effect is reported for joint mediation analyses. All marginal effects are computed on the scale of tie probabilities, not the linear component (log-odds).

}
\references{
Duxbury, Scott W. 2023. "The Problem of Scaling in Exponential Random Graph Models." \emph{Sociological Methods & Research}.

}
\author{
Scott Duxbury, Associate Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}
\note{
Note that ergm.mma may have restrictive run times in large networks or high dimensional models.

}


\seealso{
\code{\link{ergm.AME}}
\code{\link{ergm.mod.mma}}

}
\examples{

\dontshow{
require(ergm)


set.seed(21093)
a1<-network::as.network(matrix(c(rbinom(10, 1,.3),
            rbinom(10, 1,.3),
             rbinom(10, 1,.3),
              rbinom(10, 1,.3),
               rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3)),
          nrow=10,ncol=10))

network::set.vertex.attribute(a1,"var.1",rbinom(10,1,.3))

a<-ergm(a1~edges+nodeifactor("var.1"))
ab<-ergm(a1~edges+nodeifactor("var.1")+nodeofactor("var.1"))
ergm.mma(a,ab,mediator="nodeofactor.var.1.1",
             direct.effect="nodeifactor.var.1.1")
}
\donttest{
library(ergm)


data("faux.mesa.high")

faux.mesa.high

set.seed(21093)

my.ergm1<-ergm(faux.mesa.high~edges+
            nodecov("Grade")+
            nodefactor("Race")+
            nodefactor("Sex")+
            nodematch("Grade")+
            nodematch("Sex")+
            nodematch("Race"),
            estimate = "MLE")


set.seed(21093)

my.ergm2<-ergm(faux.mesa.high~edges+
                nodecov("Grade")+
                nodefactor("Race")+
                nodefactor("Sex")+
                nodematch("Grade")+
                nodematch("Sex")+
                nodematch("Race")+
               gwesp(.5,fixed=TRUE))


###gwesp as single mediator
ergm.mma(my.ergm1,my.ergm2,mediator="gwesp.fixed.0.5",
             direct.effect="nodecov.Grade")


###multiple mediators

set.seed(21093)

my.ergm3<-ergm(faux.mesa.high~edges+
                 nodecov("Grade")+
                 nodefactor("Race")+
                 nodefactor("Sex")+
                 nodematch("Grade")+
                 nodematch("Sex")+
                 nodematch("Race")+
                 gwdsp(1,fixed=TRUE)+
                 gwesp(.5,fixed=TRUE),
                 #increase MC sample to facilitate convergence
               control = control.ergm(MCMC.samplesize=5000)
)

##joint mediation for gwdsp and gwesp
ergm.mma(my.ergm2,my.ergm3,mediator=c("gwdsp.fixed.1","gwesp.fixed.0.5"),
         direct.effect=names(my.ergm3$coef[8]))
}

}
\keyword{ ~ergm }
\keyword{ ~networks }
\keyword{ ~mediation }
\keyword{ ~marginal effects }
