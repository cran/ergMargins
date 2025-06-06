\name{ergm.MEM}
\alias{ergm.MEM}
\title{
Function to compute marginal effects at means in ERGM.
}
\description{
ERGM function to compute marginal effects at means for main effects or at various levels of a moderator. Also provides tests of significance for second differences for interaction terms. Currently accepts \code{ergm}, \code{btergm}, \code{mtergm}, \code{ergm.count}, \code{mlergm}, and \code{ergm.multi} objects. See ergm.AME to compute average marginal effect.
}
\usage{
ergm.MEM(model,
         var1,
         var2=NULL,
         inter=NULL,
         at.2=NULL,
         at.controls=NULL,
         control_vals=NULL,
         return.dydx=FALSE)
}
\arguments{
  \item{model}{the ergm object to be analyzed.
}
\item{var1}{the character string containing the name of the focal variable on which to compute the MEM. A vector of character strings can be provided matching model output to compute the marginal effect of a linear combination of covariates.}
\item{var2}{the character string containing the name of the moderator variable. Default value is NULL.}
\item{inter}{the character string containing the name of the interaction between var1 and var2. Default value is NULL.}
\item{at.2}{the levels of var2 at which to compute the MEMs. If var2=NULL, this parameter is ignored. If var2!=NULL and at.2=NULL, the function computes the AME at each unique value of var2. at.2 can also be set to meaningful values, such as the mean of var2 or other discrete values of interest. If var2 is discrete, leaving at.2 NULL will compute the AME for each discrete category of the moderator variable. It is usually safe to leave at.2=NULL unless examining a large network where var2 has a large number of unique values, in which case computation times can be restrictive.}
\item{at.controls}{A character string containing the names of all control variables to be specified at fixed levels when calculating the marginal effect. Should match ERGM output exactly.}
\item{control_vals}{A vector specifying the values to hold at.controls variables. Ignored if at.controls = NULL}
\item{return.dydx}{logical parameter indicating whether to return the marginal effects for each dyad observation. Defaults to FALSE.}
}
\details{
Largely a wrapper for ergm.AME. The only difference is that marginal effects are computed at variable means.

Computes the MEM for variables specified. If var2 and inter are left at NULL, this is the simple MEM. All MEMs are computed as tie probabilities (not on the linear scale). Standard errors are computed using the Delta method. If a moderator is specified using var2 and inter, a table of MEMs is returned at each level of at.2. If at.2=NULL, this is computed for every unique level of var2.

When at.2 and var2 are specified, a table of second differences is also returned. The table of second differences assesses whether the difference in marginal effects is significant between levels of a moderator (values of at.2). If the moderator is binary, the second difference tests the overall significance of an interaction. The tests statistic is a Wald test.

If more than 2 values of at.2 are specified, the function also computes the average absolute Wald test statistic and average second difference. The average second difference measures the average change in the MEM of var1 between levels of a moderator. The average absolute Wald test is the mean absolute value of the Wald tests for all second differences. It assesses whether, on average, there is a significant difference between the second differences and zero. A significant result means the interaction is significant.

If at.controls and control_vals are provided, the MEM is calculated while holding specified control variables at fixed values provided in control_vals, rather than at their means
}
\value{
If var2 and inter=NULL, then a table is returned with the MEM for var1, its standard error, and test statistic. Note again that all MEMs are in the metric of tie probabilities, not log-odds. If a moderator is specified, a list is returned with the following elements:
\item{`Marginal effects at means`}{is the table of MEMs at each level of at.2. If at.2=NULL, the table contains unique entries for every unique value of at.2.}
\item{`Second differences`}{is a table containing point estimates, standard errors, and Wald tests for the difference in MEMs between each level of at.2. It has a number of rows equal to length(at.2)-1}
\item{`Aggregate output`}{is a table of summary statistics when at.2 has length > 2. Reports the mean second difference and mean absolute Wald statistic for the second differences. If at.2 has length<3, it is not returned}
\item{`Marginal effects`}{if return.dydx=TRUE, returns the marginal effects (partial derivatives with respect to var1) for each dyad in the ERGM sample space.}

}
\references{
Duxbury, Scott W. 2023. "The Problem of Scaling in Exponential Random Graph Models." \emph{Sociological Methods & Research}.

Duxbury, Scott W., and Jenna Wertsching. 2023. "Scaling Bias in Pooled Exponential Random Graph Models." \emph{Social Networks}.

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
ergm.MEM(a,var1="nodeifactor.var.1.1")

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
                mutual)

#main effect
main.eff<-ergm.MEM(my.ergm,"mutual")

#main effect of "mutual" when sender and receiver effects for sex == 1
ergm.MEM(my.ergm,"mutual",
              at.controls=c("nodeifactor.sex.2","nodeofactor.sex.2"),
              control_vals=c(1,1))


##with moderator
sex.inter<-ergm.MEM(my.ergm,"nodeifactor.sex.2","nodeofactor.sex.2","nodematch.sex")

##with moderator at specific levels
grade.inter<-ergm.MEM(my.ergm,"nodeicov.grade","nodeocov.grade","absdiff.grade",
          at.2=c(9,12))

##when moderator changes from mean to mean +sd
tie.probs<-edge.prob2(my.ergm)
m<-mean(tie.probs$nodeocov.grade)
m.sd<-m+sd(tie.probs$nodeocov.grade)

grade.inter.at.means<-ergm.MEM(my.ergm,"nodeicov.grade","nodeocov.grade",
                            "absdiff.grade",at.2=c(m,m.sd))

}

}
\keyword{ ~ergm }
\keyword{ ~networks }
\keyword{ ~marginal effects }
\keyword{ ~moderation }
