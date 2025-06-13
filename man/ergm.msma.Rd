\name{ergm.msma}
\alias{ergm.msma}
\title{
Function to compare marginal substructrual effects between nested ERGMs.
}
\description{
\code{ergm.msma} tests the difference in a marginal substructural effect between two ERGMs fit to the same network that differ only in their specification. It computes the total, direct, and indirect marginal substructural effects. Can be used with discrete or continuous mediators and for multiple mediator analyses. Currently accepts \code{ergm}, \code{btergm}, \code{mtergm}, \code{ergm.count}, \code{mlergm}, and \code{ergm.multi} objects.}
\usage{
ergm.msma(restricted.model,
                    full.model,
                    direct_substructural_effect,
                    higher_order_term=NULL,
                    lower_order_term=NULL,
                    at.lower_order_term=NULL,
                    mediator,
                    at.controls=NULL,
                    control_vals=NULL,
                    estimate="aMSE")
}
\arguments{
  \item{restricted.model}{the ergm without the mediator variable.
}
  \item{full.model}{the ergm including the mediator variable.}
  \item{direct_substructural_effect}{a character string identifying to the marginal substructural effect to be compared between models}
  \item{higher_order_term}{an optional parameter identifying any higher-order substructures to be held at 0 when calculating marginal substructural effects. \code{higher_order_term} must appear in both models when specified}
  \item{lower_order_term}{an optional parameter identifying lower-order substructures nested in \code{direct_substructural_effect}. When provided, the lower-order terms will be held at fixed values calculating direct and indirect marginal substructural effects}
  \item{at.lower_order_term}{optional parameter identifying the values to assign to \code{lower_order_term}. Should be provided as a vector of numeric values with each entry corresponding to a single entry in \code{lower_order_term}. Default is 1}

  \item{mediator}{a character string containing the name of the mediator variable. Can be specified as a vector with each entry containing the character strings for several mediator variables when there are multiple mediators.}
  \item{at.controls}{a vector of character strings containing the names of control variables to hold at prespecified values.}
  \item{control_vals}{a vector of numeric values to hold at.controls at during estimation.}
  \item{estimate}{the type of marginal substructural effect to calculate. Must be one of \code{"aMSE"}, \code{"MSEm"},\code{"tMSE"}, or \code{"tMSEm"}}
}
\details{
\code{ergm.msma} calculates the change in a marginal substructural effect between two models fit to the same data. Formally, define two models \eqn{Mod^1} and \eqn{Mod^2} that differ by the inclusion of one or more mediating variables. For a substructural network term of interest given by the change statistic \eqn{\delta ^+_{ij}(y)}, the \code{ergm.msma} provides the following quantities:
\deqn{MSE^{Mod^1}(\delta ^{+}_{ij}(y))}
\deqn{MSE^{Mod^2}(\delta ^{+}_{ij}(y))}
\deqn{MSE^{Mod^1}(\delta ^{+}_{ij}(y))-MSE^{Mod^2}(\delta ^{+}_{ij}(y))}

Where \eqn{MSE^{Mod^1}(\delta ^{+}_{ij}(y))} is the total effect corresponding to the marginal substructural effect of \eqn{\delta ^{+}_{ij}(y)} operating through all indirect pathways, \eqn{MSE^{Mod^2}(\delta ^{+}_{ij}(y))} is the direct effect of \eqn{\delta ^+_{ij}(y)} on tie formation, and \eqn{MSE^{Mod^1}(\delta ^{+}_{ij}(y))-MSE^{Mod^2}(\delta ^{+}_{ij}(y))} is the portion of the total effect explained by controlling for a confounding variable.

}
\value{
Returns a table including the point estimates and standard errors for total, direct, and indirect effects. A single indirect effect is reported for joint mediation analyses. All marginal substructural effects are computed on the scale of tie probabilities, not the linear component (log-odds).

}
\references{
Duxbury, Scott W. 2023. "The Problem of Scaling in Exponential Random Graph Models." \emph{Sociological Methods & Research}.

Duxbury, Scott W., and Jenna Wertsching. 2023. "Scaling Bias in Pooled Exponential Random Graph Models." Social Networks.

Duxbury, Scott W. 2025. "Marginal Substructural Effects for Exponential Random Graph Models." Working Paper.
}
\author{
Scott Duxbury, Associate Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}
\note{
Note that ergm.msma may have restrictive run times in large networks or high dimensional models.

}


\seealso{
\code{\link{ergm.AME}}
\code{\link{ergm.mod.mma}}
\code{\link{ergm.mma}}
\code{\link{ergm.MSE}}
\code{\link{ergm.MEM}}

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
ergm.msma(a,ab,
             direct_substructural_effect="nodeifactor.var.1.1",
             mediator="nodeofactor.var.1.1")
}
\donttest{
library(ergm)


data("faux.mesa.high")

faux.mesa.high

#does sociality explain triadic effects?

  #model without sociality effects
my.ergm1<-ergm(faux.mesa.high~edges+
                nodematch("Grade")+
                nodematch("Sex")+
                nodematch("Race")+
                gwesp(.5,fixed=TRUE),
               control = control.ergm(MCMLE.termination="Hummel",
                                      MCMLE.effectiveSize = NULL,
                                      MCMC.samplesize = 5000,
                                      MCMLE.maxit=120))


  #model with sociality effects
my.ergm2<-ergm(faux.mesa.high~edges+
                nodecov("Grade")+
                nodefactor("Race")+
                nodefactor("Sex")+
                nodematch("Grade")+
                nodematch("Sex")+
                nodematch("Race")+
                 gwesp(.5,fixed=TRUE),
                control = control.ergm(MCMLE.termination="Hummel",
                                      MCMLE.effectiveSize = NULL,
                                      MCMC.samplesize = 5000,
                                      MCMLE.maxit=120))


###sociality effects as mediators
ergm.msma(my.ergm1,
            my.ergm2,
            direct_substructural_effect="gwesp.fixed.0.5",
            mediator=c("nodecov.Grade","nodefactor.Sex.M","nodefactor.Race.White",
            "nodefactor.Race.Other","nodefactor.Race.Hisp","nodefactor.Race.NatAm"),
            estimate="aMSE")


}

}
\keyword{ ~ergm }
\keyword{ ~networks }
\keyword{ ~mediation }
\keyword{ ~marginal effects }
