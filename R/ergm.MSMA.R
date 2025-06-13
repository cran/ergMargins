#' Function to conduct marginal mediation analysis
#'
#'computes indirect effect using average marginal effects with bootstrapped standard errors
#' @param restricted.model is the model without the mediator variable(s)
#' @param full.model is the full model including all mediating variables
#' @param mediator is a mediating varaible or a character vector containing the names of the mediator variables if joint mediation is of interest
#' @param direct_substructural_effect is the direct effect of interest


ergm.msma<-function(restricted.model,
                    full.model,
                    direct_substructural_effect,
                    higher_order_term=NULL,
                    lower_order_term=NULL,
                    at.lower_order_term=NULL,
                    mediator,
                    at.controls=NULL,
                    control_vals=NULL,
                    estimate="aMSE"){


  if(class(restricted.model)%in%"btergm"){
    out<-ergm.msma_boot(restricted.model=restricted.model,
                        full.model=full.model,
                        direct_substructural_effect=direct_substructural_effect,
                        higher_order_term=higher_order_term,
                        lower_order_term=lower_order_term,
                        at.lower_order_term=at.lower_order_term,
                        mediator=mediator,
                        at.controls=at.controls,
                        control_vals=control_vals,
                        estimate=estimate)
    return(out)
  }


  if(class(restricted.model)%in%"mlergm"){
    theta1<-restricted.model$theta
    theta2<-full.model$theta

  }else{
    theta1<-btergm::coef(restricted.model)
    theta2<-btergm::coef(full.model)

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


   # dse<-direct_substructural_effect
    tot.MSE<-ergm.MSE(model=restricted.model,
                      substructural_effect=direct_substructural_effect,
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
    tot.Jac<-tot.MSE$Jac
    tot.MSE<-tot.MSE$MSE

    p.MSE<-ergm.MSE(model=full.model,
                    substructural_effect=direct_substructural_effect,
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
    p.Jac<-p.MSE$Jac
    p.MSE<-p.MSE$MSE




  ###estimate cross-model covariance
  names(p.Jac)<-names(theta2)
  names(tot.Jac)<-names(theta1)
  tot.Jac_expanded<-rep(0,length(p.Jac))
  names(tot.Jac_expanded)<-names(p.Jac)
  tot.Jac_expanded[names(tot.Jac_expanded)%in%names(tot.Jac)]<-tot.Jac

  #get vcovs
  if(class(restricted.model)%in%"mtergm"){
    tot.vc <- stats::vcov(restricted.model@ergm)
    tot.vc <-tot.vc[!rownames(tot.vc)%in%c("edgecov.offsmat","offset(edgecov.offsmat)"),
                    !colnames(tot.vc)%in%c("edgecov.offsmat","offset(edgecov.offsmat)")]
    p.vc <- stats::vcov(full.model@ergm)
    p.vc <-p.vc[!rownames(p.vc)%in%c("edgecov.offsmat","offset(edgecov.offsmat)"),
                !colnames(p.vc)%in%c("edgecov.offsmat","offset(edgecov.offsmat)")]

  }else{

    tot.vc <- stats::vcov(restricted.model)
    p.vc <- stats::vcov(full.model)

  }

  if(class(restricted.model)%in%"mlergm"){
    tot.vc<-solve(tot.vc)
    p.vc<-solve(p.vc)
  }


  tot_vc_expanded<-matrix(0,nrow=nrow(p.vc),ncol=ncol(p.vc))
  rownames(tot_vc_expanded)<-colnames(tot_vc_expanded)<-names(theta2)
  tot_vc_expanded[rownames(tot_vc_expanded)%in%names(theta1),
                  colnames(tot_vc_expanded)%in%names(theta1)]<-tot.vc


  #calculate gradient and cross_mod covariance
  cross_mod_grad<-t(p.Jac)%*%t(as.matrix(tot.Jac_expanded))

  cross_mod_cov<-p.vc%*%cross_mod_grad%*%tot_vc_expanded


  rownames(tot.MSE)<-paste("total effect:",rownames(tot.MSE))
  rownames(p.MSE)<-paste("partial effect:",rownames(p.MSE))

  ###indirect effect

  msma.me<-tot.MSE[1,1]-p.MSE[1,1]
  cov.mse<-cross_mod_cov[direct_substructural_effect,direct_substructural_effect]
  msma.se<-sqrt(tot.MSE[,2]^2+p.MSE[,2]^2-cov.mse)
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
