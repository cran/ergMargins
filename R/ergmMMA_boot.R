#' Function to conduct marginal mediation analysis
#'
#'computes indirect effect using average marginal effects with bootstrapped standard errors
#' @param restricted.model is the model without the mediator variable(s)
#' @param full.model is the full model including all mediating variables
#' @param mediator is a mediating varaible or a character vector containing the names of the mediator variables if joint mediation is of interest
#' @param direct.effect is the direct effect of interest


ergm.mma_boot<-function(restricted.model,full.model,direct.effect,mediator,
                   at.controls=NULL,
                   control_vals=NULL,
                   ME="AME"){

  if(!ME%in%c("AME","MEM")){
    warning("ME must be specified as AME or MEM. Returning the AME.")
    ME<-"AME"
  }


      theta1<-btergm::coef(restricted.model)
      theta2<-btergm::coef(full.model)



  ##check at.controls appear in both modelsz
  if(!is.null(at.controls)){
    if(!all(at.controls%in%names(theta1))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
    if(!all(at.controls%in%names(theta2))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
  }


  if(ME%in%"AME"){

    tot.AME<-ergm.AME(restricted.model,direct.effect,return.dydx=TRUE,
                    at.controls=at.controls,control_vals=control_vals)
    tot.bootstraps<-tot.AME$bootstraps
    tot.AME<-tot.AME$AME

    p.AME<-ergm.AME(full.model,direct.effect,return.dydx=TRUE,
                  at.controls=at.controls,control_vals=control_vals)
    p.bootstraps<-p.AME$bootstraps
    p.AME<-p.AME$AME
  }else{
    tot.AME<-ergm.MEM(restricted.model,direct.effect,return.dydx=TRUE,
                      at.controls=at.controls,control_vals=control_vals)
    tot.bootstraps<-tot.AME$bootstraps
    tot.AME<-tot.AME$MEM

    p.AME<-ergm.MEM(full.model,direct.effect,return.dydx=TRUE,
                    at.controls=at.controls,control_vals=control_vals)
    p.bootstraps<-p.AME$bootstraps
    p.AME<-p.AME$MEM

  }

    #ensure equal dimensions. Typically only necessary in large networks.
 # if(length(tot.dydx)!=length(p.dydx)){

  #  if(length(tot.dydx)<length(p.dydx)){
 #     p.dydx<-p.dydx[1:length(tot.dydx)]
  #  }else{
#      tot.dydx<-tot.dydx[1:length(p.dydx)]
  #  }

 # }



  rownames(tot.AME)<-paste("total effect:",rownames(tot.AME))
  rownames(p.AME)<-paste("partial effect:",rownames(p.AME))

    ###indirect effect

    mma.me<-tot.AME[1,1]-p.AME[1,1]
    mma.se<-sd(tot.bootstraps-p.bootstraps)
    mma.z<-mma.me/mma.se
    p.mma<-2*stats::pnorm(-abs(mma.z))

  ind<-matrix(signif(c(mma.me,mma.se,mma.z,p.mma),digits=5),nrow=1,ncol=4)
  colnames(ind)<-colnames(tot.AME)
   if(length(mediator)>1) {
    mediator<-paste(mediator,collapse=", ")}
  rownames(ind)<-paste("indirect effect:",direct.effect,"->",mediator)

  out<-rbind(tot.AME,p.AME,ind)
  proportion.mediated<-1-(p.AME[1]/tot.AME[1])
  attr(out,"description")<-paste("proportion of",direct.effect,"mediated by",mediator," = ",round(proportion.mediated,digits=3))


  return(out)
}
