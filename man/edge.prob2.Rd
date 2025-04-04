\name{edge.prob2}
\alias{edge.prob2}

\title{Function to compute tie probabilities for each dyad.
}
\description{Computes tie probabilities for each dyad. Largely a wrapper for edgeprob() in btergm.
}
\usage{
edge.prob2(model,verbose=FALSE)
}
\arguments{
  \item{model}{ model is the ergm object.
}
  \item{verbose}{Boolean operator. If set to T, updates user on progress of edge.prob2. May be useful for large networks.}
}
\details{
Computes tie probabilities for each dyad in the ERGM sample space.This function is largely a wrapper for the edgeprob function in btergm, but can be applied to curved ERGMs.
}
\value{
Returns a dataframe where weach row is a dyad and each column is a term in the model, an identification variable for the vertices/dyads, or the estimated probability of observing a tie. See edgeprob() in btergm for more details.
}

\author{
Scott Duxbury, Associate Professor, University of North Carolina--Chapel Hill, Department of Sociology.
original author is Philip Leifeld at University of Essex.
}


\seealso{
\code{\link[btergm]{edgeprob}}
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
edge.prob2(a)

}
\donttest{
library(ergm)


data("faux.dixon.high")


my.ergm<-ergm(faux.dixon.high~edges+
                nodeifactor("sex")+
                nodeofactor("sex")+
                nodematch("sex")+
                mutual)
edge.prob2(my.ergm)
}
}
\keyword{ ~ERGM }
\keyword{ ~tie probabilities }
