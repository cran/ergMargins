#' Function to conduct marginal mediation analysis
#'
#'computes indirect effect using average marginal effects with bootstrapped standard errors
#' @param restricted.model is the model without the mediator variable(s)
#' @param full.model is the full model including all mediating variables
#' @param mediator is a mediating varaible or a character vector containing the names of the mediator variables if joint mediation is of interest
#' @param direct_substructural_effect is the direct effect of interest


ergm.msma_boot<-function(restricted.model,
                    full.model,
                    direct_substructural_effect,
                    higher_order_term=NULL,
                    lower_order_term=NULL,
                    at.lower_order_term=NULL,
                    mediator,
                    at.controls=NULL,
                    control_vals=NULL,
                    estimate="aMSE"){

  
  theta1<-btergm::coef(restricted.model)
  theta2<-btergm::coef(full.model)
  
  ##check at.controls appear in both models
  if(!is.null(at.controls)){
    if(!all(at.controls%in%names(theta1))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
    if(!all(at.controls%in%names(theta2))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
  }
  
  

  
  ##check at.controls appear in both models
  if(!is.null(at.controls)){
    if(!all(at.controls%in%names(theta1))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
    if(!all(at.controls%in%names(theta2))){
      stop("Variables specified in any.controls must appear in both models to be used.")
    }
  }

  

    tot.MSE<-ergm.MSE(restricted.model,
                      direct_substructural_effect,
                      higher_order_term=higher_order_term,
                      lower_order_term=lower_order_term,
                      at.lower_order_term=at.lower_order_term,
                      var2=NULL,
                      inter=NULL,
                      at.2=NULL,
                      at.controls=at.controls,
                      control_vals=control_vals,
                      estimate=estimate,
                      return_Jac=TRUE)
    tot.BS<-tot.MSE$bootstraps
    tot.MSE<-tot.MSE$MSE
    
    p.MSE<-ergm.MSE(full.model,
                    direct_substructural_effect,
                    higher_order_term=higher_order_term,
                    lower_order_term=lower_order_term,
                    at.lower_order_term=at.lower_order_term,
                    var2=NULL,
                    inter=NULL,
                    at.2=NULL,
                    at.controls=at.controls,
                    control_vals=control_vals,
                    estimate=estimate,
                    return_Jac=TRUE)
    p.BS<-p.MSE$bootstraps
    p.MSE<-p.MSE$MSE

  


  rownames(tot.MSE)<-paste("total effect:",rownames(tot.MSE))
  rownames(p.MSE)<-paste("partial effect:",rownames(p.MSE))
  
  ###indirect effect
  
  msma.me<-tot.MSE[1,1]-p.MSE[1,1]
  msma.se<-sd(tot.BS-p.BS,na.rm=T)
  msma.z<-msma.me/msma.se
  p.msma<-2*stats::pnorm(-abs(msma.z))
  
  ind<-matrix(signif(c(msma.me,msma.se,msma.z,p.msma),digits=5),nrow=1,ncol=4)
  colnames(ind)<-colnames(tot.MSE)
  if(length(mediator)>1) {
    mediator<-paste(mediator,collapse=", ")}
  rownames(ind)<-paste("indirect effect:",direct_substructural_effect,"->",mediator)
  
  out<-rbind(tot.MSE,p.MSE,ind)
  proportion.mediated<-1-(p.MSE[1]/tot.MSE[1])
  attr(out,"description")<-paste("proportion of",direct_substructural_effect,"mediated by",mediator," = ",round(proportion.mediated,digits=3))
  
  
  return(out)
}
