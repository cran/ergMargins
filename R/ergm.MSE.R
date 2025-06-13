#' Function to compute marginal substructural effects in ERGM
#'
#' @param model is an ergm object
#' @param substructural_effect is the name of the substructural effect of interest, character string
#' @param higher_order_term is a vector listing names of all higher order structural dependencies, character string
#' @param lower_order_term is a vector listing names of all lower order structural dependencies, character string
#' @param at.lower_order_term is a vector specifying the values to assign lower_order_term. If left NULL, it defaults to 1
#' @param estimate defines the MSE summary statistic. Must be equal to "aMSE", "MSEm", "tMSE", or "tMSEm". Default is "aMSE"
#' @param var2 optional parameter when examining interactions with exogenous attributes. Corresponds to name of moderating term
#' @param inter optional parameter when examining interactions with exogenous attributes. Corresponds to name of interaction
#' @param at.2 optional parameter to assign values to variables in var2. Default is c(0,1)
#' @param at.controls optional parameter to hold controls at pre-specified levels. Expects a vector of character strings corresponding to the names for each control variable
#' @param control_vals the values to assign to variables in at.controls. Optional parameter



ergm.MSE<-function(model,
                   substructural_effect,
                   higher_order_term=NULL,
                   lower_order_term=NULL,
                   at.lower_order_term=NULL,
                   estimate="aMSE", #must be one of c("aMSE","MSEm","tMSE","tMSEm")
                   var2=NULL,
                   inter=NULL,
                   at.2=NULL,
                   at.controls=NULL,
                   control_vals=NULL,
                   return_Jac=FALSE){
  if(!estimate%in%c("aMSE","MSEm","tMSEm","tMSE")){
    message("Value of estimate argument not recognized. Returning the MSEm")
    estimate<-"aMSE"
  }

  if(!is.null(lower_order_term) & is.null(at.lower_order_term)){
    at.lower_order_term<-rep(1,length(lower_order_term))
  }

  if(is.null(at.2)&!is.null(inter)){
    at.2<-c(0,1)
  }

  ####start prepackage

  ##need to incorporate ergm.MSE_boot and ergm.MSE_count
  if(class(model)[1]%in%"btergm"){
    out<-ergm.MSE_boot(model,
                       substructural_effect,
                       higher_order_term=higher_order_term,
                       lower_order_term=lower_order_term,
                       at.lower_order_term=at.lower_order_term,
                       estimate=estimate,
                       var2=var2,
                       inter=inter,
                       at.2=at.2,
                       at.controls=at.controls,
                       control_vals=control_vals,
                       return_Jac = return_Jac)
    return(out)
  }

  if(!class(model)[1]%in%"mtergm"){
    if(is.valued(model)){
      out<-ergm.MSE_count(model,
                          substructural_effect,
                          higher_order_term=higher_order_term,
                          lower_order_term=lower_order_term,
                          at.lower_order_term=at.lower_order_term,
                          estimate=estimate,
                          var2=var2,
                          inter=inter,
                          at.2=at.2,
                          at.controls=at.controls,
                          control_vals=control_vals,
                          return_Jac=return_Jac)

      return(out)
    }
  }

  if(class(model)[1]%in%"mtergm"){
    dyad.mat<-ergmMPLE(model@ergm$formula,output="dyadlist",basis=model@ergm$network)
    dyad.mat<-dyad.mat$predictor[,-c(1:2,ncol(dyad.mat$predictor))]
    vc <- stats::vcov(model@ergm)
    vc<-vc[!rownames(vc)%in%"offset(edgecov.offsmat)",
           !colnames(vc)%in%"offset(edgecov.offsmat)"]


  }else{

    dyad.mat<-ergmMPLE(model$formula,output="dyadlist",basis=model$network,
                       control=control.ergm(
                         term.options = list(interact.dependent="message")))
    dyad.mat<-dyad.mat$predictor[,-c(1:2)]
    vc <- stats::vcov(model)

  }


  if(class(model)[1]%in%"mlergm"){
    theta<-model$theta
    vc<-solve(vc) #invert the fisher matrix
  }else{
    theta<-btergm::coef(model)
  }
  #handle mlergm objects
  if("mlergm"%in%class(model)){
    class(model)<-"ergm"
  }

  ##handle curved ergms by removing decay parameter
  #note that the micro-level change statistics are already properly weighted,
  #so decay term is not needed for predictions

  if(class(model)[1]%in%"mtergm"){

    if(ergm::is.curved(model@ergm)){
      curved.term<-curved.term_main<-vector(length=length(model@ergm$etamap$curved))
      curved_loc_from<-curved_loc_to<-list()
      for(i in 1:length(model@ergm$etamap$curved)){
        curved.term_main[i]<-model@ergm$etamap$curved[[i]]$to[1]
        curved.term[i]<-model@ergm$etamap$curved[[i]]$from[2]
        curved_loc_from[[i]]<-model@ergm$etamap$curved[[i]]$from
        curved_loc_to[[i]]<-model@ergm$etamap$curved[[i]]$to

      }

      curved_loc_from_list<-unlist(curved_loc_from)
      theta<-theta[-c(curved.term)]

      all_drops<-unlist(curved_loc_to)
      all_drops<-all_drops[!all_drops%in%c(curved.term_main)]
      dyad.mat<-dyad.mat[,-c(all_drops)]
      vc<-vc[-c(curved.term),-c(curved.term)]


    }

  }else{
    if(ergm::is.curved(model)){
      curved.term<-curved.term_main<-vector(length=length(model$etamap$curved))
      curved_loc_from<-curved_loc_to<-list()
      for(i in 1:length(model$etamap$curved)){
        curved.term_main[i]<-model$etamap$curved[[i]]$to[1]
        curved.term[i]<-model$etamap$curved[[i]]$from[2]
        curved_loc_from[[i]]<-model$etamap$curved[[i]]$from
        curved_loc_to[[i]]<-model$etamap$curved[[i]]$to

      }

      curved_loc_from_list<-unlist(curved_loc_from)
      theta<-theta[-c(curved.term)]
      all_drops<-unlist(curved_loc_to)
      all_drops<-all_drops[!all_drops%in%c(curved.term_main)]
      dyad.mat<-dyad.mat[,-c(all_drops)]
      vc<-vc[-c(curved.term),-c(curved.term)]


    }
  }




  #handle offset
  offset_ind<-pmatch("offset",names(theta))
  if(!is.na(offset_ind)){
    if(!class(model)[1]%in%"ergm.ego"){
      dyad.mat<-dyad.mat[,-offset_ind]
    }
    theta<-theta[-offset_ind]
    vc<-vc[-offset_ind,-offset_ind]
  }

  ##check for correct names
  if(any(!substructural_effect%in%names(theta))){
    stop("at least one substructural_effect not found in model. Respecify")
  }

  if(!is.null(lower_order_term)&any(!lower_order_term%in%names(theta))){
    stop("at least one substructural_effect not found in model. Respecify")
  }

  if(!is.null(higher_order_term)&any(!higher_order_term%in%names(theta))){
    stop("at least one substructural_effect not found in model. Respecify")
  }


  if(any(names(theta)!=colnames(dyad.mat))){
    colnames(dyad.mat)<-names(theta) #make sure names align
  }


  #end prepackage



  ##start value assignments
    #tMSE is the focal effect + change in all lower_order_term
  if(estimate%in%c("tMSE","tMSEm")){
    substructural_effect<-c(substructural_effect,lower_order_term)
    lower_order_term<-NULL
  }


  ##assign fixed values for controls when specified
  substr_controls<-NULL
  substr_vals<-NULL

  if(is.null(lower_order_term)&!is.null(higher_order_term)){
    substr_controls<-higher_order_term
    substr_vals<-rep(0,length(higher_order_term))
  }

  if(!is.null(lower_order_term)&is.null(higher_order_term)){
    substr_controls<-lower_order_term
    substr_vals<-at.lower_order_term
  }

  if(!is.null(lower_order_term)&!is.null(higher_order_term)){
    substr_controls<-c(higher_order_term,lower_order_term)
    substr_vals<-c(rep(0,length(higher_order_term)),
                   at.lower_order_term)

  }


  if(is.null(at.controls)){
    at.controls<-substr_controls
    control_vals<-substr_vals
  }else{
    at.controls<-c(substr_controls,at.controls)
    control_vals<-c(substr_vals,control_vals)
  }

  if(!is.null(at.controls)){
    if(is.null(control_vals)){
      stop("control_vals must be specified to use at.controls argument.")
    }

    if(length(at.controls)==1){
      dyad.mat[,at.controls]<-control_vals
    }else{
      for(i in 1:length(at.controls)){
        dyad.mat[,at.controls][,i]<-control_vals[i]
      }
    }


  }


  #MSEm holds all values at means
  if(estimate%in%c("MSEm","tMSEm")){
    dyad.mat<-t(as.matrix(colMeans(dyad.mat)))
  }


  #end value assignments


  #start MSE calculation--no moderator
    #accomodate moderator
  loop_length<-ifelse(is.null(inter),1,length(at.2))

  for(i in 1:loop_length){

   dyad_mat0<-dyad.mat
   dyad_mat0[,colnames(dyad_mat0)%in%substructural_effect]<-0
   dyad_mat1<-dyad.mat
   dyad_mat1[,colnames(dyad_mat1)%in%substructural_effect]<-1

   if(loop_length>1){
     dyad_mat0[,colnames(dyad_mat0)%in%c(var2,inter)]<-at.2[i]
     dyad_mat1[,colnames(dyad_mat1)%in%c(var2,inter)]<-at.2[i]
   }



   MSE_fun<-function(theta){
     p_0<-1/(1+exp(-(dyad_mat0%*%theta)))
     p_1<-1/(1+exp(-(dyad_mat1%*%theta)))

     mean(p_1-p_0,na.rm=T)

   }


   #return(list(theta,dyad_mat0))


    #no moderator
   if(i==1){
     MSE<-MSE_fun(theta)
     J<-numDeriv::jacobian(MSE_fun,theta)
     Jac<-as.matrix(J)
    # return(list(Jac,vc))
     variance.MSE<-J%*%vc%*%t(J)
     MSE.se<-sqrt(variance.MSE)

   }else{
     #with moderator
     MSE<-c(MSE,MSE_fun(theta))
     J<-numDeriv::jacobian(MSE_fun,theta)
     Jac<-rbind(Jac,as.matrix(J))

     return(list(Jac,vc))

     variance.MSE<-J%*%vc%*%t(J)

     MSE.se<-c(MSE.se,sqrt(variance.MSE))

   }
  }

  MSE.z<-MSE/MSE.se
  P.MSE<-2*(stats::pnorm(-abs(MSE.z)))

  if(i==1){
    MSE<-matrix(c(MSE,MSE.se,MSE.z,P.MSE),nrow=1,ncol=4)
  }else{
    MSE<-cbind(MSE,MSE.se,MSE.z,P.MSE)
  }

  colnames(MSE)<-c(estimate,"Delta SE","Z","P")

  if(!is.null(inter)){
    rownames(MSE)<-paste(substructural_effect,"with",var2,"held at",at.2,sep=" ")
  }else{
    rownames(MSE)<-ifelse(length(substructural_effect)>1,
                          paste(substructural_effect,collapse="+"),
                          substructural_effect)
  }



  MSE<-signif(MSE,digits=5)

  #stop if no interaction
  if(is.null(inter)){
    if(isTRUE(return_Jac)){
      MSE<-list(MSE=MSE,
                Jac=Jac)
      return(MSE)
    }else{
       return(MSE)
    }
  }

  #yes, interaction
  SD<-diff(MSE[,1])
  covar<-Jac%*%vc%*%t(Jac)
  se<-1:(loop_length-1)
  for(i in 1:length(se)){
    k<-i+1
    se[i]<-sqrt(covar[i,i]+covar[k,k]-(2*covar[i,k]))
  }

  z<-SD/se
  P_SD<-2*(stats::pnorm(-abs(z)))

  SD_dat<-matrix(c(SD,se,z,P_SD),nrow=1,ncol=4)
  colnames(SD_dat)<-c(estimate,"Delta SE","Z","P")
  rownames(SD_dat)<-paste(at.2[-c(length(at.2))],"to",at.2[-c(1)])

  if(isTRUE(return_Jac)){

    MSE<-list(change_in_MSE=SD_dat,
              MSE=MSE,
              Jac=Jac)
  }else{
    MSE<-list(change_in_MSE=SD_dat,
              MSE=MSE)

  }

  return(MSE)



}
