\name{ergm.MSE}
\alias{ergm.MSE}
\title{
Function to compute marginal substructural effects in ERGM.
}
\description{
ERGM function to compute marginal substructural effects in ERGM, as described by Duxbury (2025). Intended for use with endogenous substructural terms in ERGMs. Also provides utilities to test interactions between substructural terms and exogenous covariates. See \code{\link{ergm.AME}} for average marginal effects and \code{\link{ergm.MEM}} for marginal effects at means.}

\usage{
ergm.MSE(model,
         substructural_effect,
         higher_order_term=NULL,
         lower_order_term=NULL,
         at.lower_order_term=NULL,
         estimate="aMSE",
         var2=NULL,
         inter=NULL,
         at.2=NULL,
         at.controls=NULL,
         control_vals=NULL,
         return_Jac=FALSE)
}
\arguments{
  \item{model}{the model to obtain MSE estimates.Currently accepts \code{ergm}, \code{btergm}, \code{mtergm}, \code{ergm.count}, \code{mlergm}, and \code{ergm.multi} objects.
}
\item{substructural_effect}{the character string containing the name of the endogenous susbtructural term to use to estimate the MSE.}
\item{higher_order_term}{character string or vector of character strings containing the name(s) of all higher-order endogenous substructures. Optional value. Default value is NULL.}
\item{lower_order_term}{character string or vector of character strings containing the name(s) of all lower-order endogenous substructures. Terms provided will be held at a fixed value as provided with the \code{at.lower_order_term} argument.}
\item{at.lower_order_term}{a vector of numeric values to assign to variables provided in \code{lower_order_term} argument. \code{lower_order_term} will be held at those values when calculating MSE. Should be a 1:1 correspondence between the terms in \code{at.lower_order_term} and \code{lower_order_term}. Defaults to 1 if left \code{NULL}}
\item{estimate}{the type of MSE summary statistic to compute. Should be one of "aMSE", "MSEm","tMSE", or "tMSEm". See details below for information on each estimand. Default is "aMSE"}
\item{var2}{the character string containing the name of an exogenous variable that interactions with \code{substructural_term}. Used when there the researcher provides an interaction between an endogenous graph statistic and an exogenous attribute. Default value is \code{NULL}. Optional; can be ignored if no interactions are specified.}
\item{inter}{the character string containing the interaction between \code{substructural_term} and the exogenous covariate assigned to \code{var2}. Default value is \code{NULL}. Optional; can be ignored if no interactions are specified.}
\item{at.2}{the levels of \code{var2} at which to compute the MSE in cases where there is an interaction between \code{substructural_term} and \code{var2}. Provided as a vector of values. If var2=NULL, this parameter is ignored. If var2!=NULL and at.2=NULL, computed holding \code{var2} at \code{c(0,1)}.}
\item{at.controls}{A character string containing the names of all control variables to be specified at fixed levels when calculating the MSE. Should match ERGM output exactly.}
\item{control_vals}{A vector specifying the values to hold at.controls variables. Ignored if at.controls = NULL}
\item{return_Jac}{logical parameter indicating whether to return the Jacobian matrix of parital derivatives used to compute the MSE covariance matrix. Defaults to FALSE.}
}





\details{

Computes the marginal substructural effect for a specified substructural term. Formally, the MSE for a unit is given by

\deqn{MSE_{ij}=Pr(Y_{ij}=1|-y_{ij},\delta ^+_{ij}(y)=1,x_{ij})-Pr(Y_{ij}=0|-y_{ij},\delta ^+_{ij}(y)=0,x_{ij})}

where \eqn{\delta ^+_{ij}(y)} is the substructural effect of interest. The \code{estimate} argument determines whether the MSE should be taken holding controls at observed values (when set to \code{"aMSE"}), holding controls at means (\code{"MSEm"}), or by combining the estimates from multiple substructures (when set to \code{"aTSE"}).

The terms provided in the \code{lower_order_terms} argument will be held constant at a value of 1 when estimating the MSE, unless alternative values are provided with the \code{at.lower_order_term} argument. If the \code{higher_order_term} argument is provided, the specified term will be held at a value of 0 when estimating the marginal substructural effect

When an endogenous network term is specified as an interaction with an exogenous term, the MSE for each level of the interaction can be obtained by providing the \code{var2} and \code{inter} argumnets. This returns the MSE at different levels of \code{var2} along with a test of the differences between the MSEs. If at.controls and control_vals are provided, the MEM is calculated while holding specified control variables at fixed values provided in control_vals, rather than at their means
}
\value{
A table containing the MSE estimate and variance estimates. If \code{inter} and \code{var2} are provided, a list is returned containing:

\item{change_in_MSE`}{is the change in MSE when the value of \code{var2} changes in value, usually from 0 to 1 unless specified with the \code{at.2} argument}
\item{MSE}{is a table containing MSE estimates, standard errors, and Wald tests for the MSE at each level of \code{var2}}

}
\references{
Duxbury, Scott W. 2023. "The Problem of Scaling in Exponential Random Graph Models." \emph{Sociological Methods & Research}.

Duxbury, Scott W., and Jenna Wertsching. 2023. "Scaling Bias in Pooled Exponential Random Graph Models." Social Networks.

Duxbury, Scott W. 2025. "Marginal Substructural Effects for Exponential Random Graph Models." Working Paper.

Long, J. Scott, and Sarah Mustillo. 2018. "Using Predictions and Marginal Effects to Compare Groups in Regression Models for Binary Outcomes." \emph{Sociological Methods and Research}

Mize, Trenton D. 2019. "Best Practices for Estimating, Interpreting, and Presenting Nonlinear Interaction Effects." \emph{Sociological Science}
}
\author{
Scott Duxbury, Associate Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}
\note{
Syntax is loosely based on the margins package by Thomas Leeper and Stata's margins commands. See ergm.AME to compute average marginal effects. Currently only considers two-ways interactions.}


\seealso{
\code{\link{ergm.AME}}
\code{\link{ergm.mma}}
\code{\link{ergm.msma}}
\code{\link{ergm.MEM}}
\code{\link[margins]{margins}}
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

a<-ergm(a1~edges+nodeifactor("var.1")+nodeofactor("var.1"))
ergm.MSE(a,"nodeifactor.var.1.1",estimate="MSEm")

}
\donttest{

library(ergm)


data("faux.dixon.high")


set.seed(21093)

my.ergm<-ergm(faux.dixon.high~edges+
                nodeicov("grade")+
                nodeocov("grade")+
                nodeifactor("sex")+
                nodeofactor("sex")+
                absdiff("grade")+
                nodematch("sex")+
                mutual+
                gwidegree(.5,fixed=TRUE))

#MSE at means
ergm.MSE(my.ergm,
          substructural_effect="mutual",
          lower_order_term="gwideg.fixed.0.5",
          estimate="MSEm")

#total effect of both endogenous terms
ergm.MSE(my.ergm,
          substructural_effect="mutual",
          lower_order_term="gwideg.fixed.0.5",
          estimate="tMSEm")


}

}
\keyword{ ~ergm }
\keyword{ ~networks }
\keyword{ ~marginal effects }
\keyword{ ~moderation }
