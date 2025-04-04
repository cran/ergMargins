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




ergm.MSE_boot<-function(model,
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
                   return_Jac=TRUE){

  dyad.mat<-edge.prob2(model)[,-c(1)]
  dyad.full.mat<-dyad.mat
  p<-dyad.mat$probability
  start.drops<-ncol(dyad.mat)-5
  dyad.mat<-dyad.mat[,-c(start.drops:ncol(dyad.mat))]
  theta<-btergm::coef(model)
  if(any(names(theta)!=colnames(dyad.mat))){
    colnames(dyad.mat)<-names(theta) #make sure names align
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
  theta_boot<-model@boot$t

  #end value assignments


  #start MSE calculation
  loop_length<-ifelse(is.null(inter),1,length(at.2))

  for(i in 1:loop_length){

   dyad_mat0<-as.matrix(dyad.mat)
   dyad_mat0[,colnames(dyad_mat0)%in%substructural_effect]<-0
   dyad_mat1<-as.matrix(dyad.mat)
   dyad_mat1[,colnames(dyad_mat1)%in%substructural_effect]<-1

   if(loop_length>1){
     dyad_mat0[,colnames(dyad_mat0)%in%c(var2,inter)]<-at.2[i]
     dyad_mat1[,colnames(dyad_mat1)%in%c(var2,inter)]<-at.2[i]
   }

 #  return(dim(theta))
   MSE_fun<-function(theta){
     theta<-as.matrix(theta)
    # return(dim(theta))
     p_0<-1/(1+exp(-(dyad_mat0%*%theta)))
     p_1<-1/(1+exp(-(dyad_mat1%*%theta)))

     mean(p_1-p_0,na.rm=T)

   }


   MSE_boots<-apply(theta_boot,1,FUN=MSE_fun)

    #no moderator
   if(i==1){
     MSE<-MSE_fun(theta)
     MSE.se<-sd(MSE_boots)
     bootMat<-matrix(MSE_boots)


   }else{
     #with moderator
     MSE<-c(MSE,MSE_fun(theta))
     MSE.se<-c(MSE.se,sd(MSE_boots))
     bootMat<-cbind(bootMat,matrix(MSE_boots))

   }
  }

  MSE.z<-MSE/MSE.se
  P.MSE<-2*(stats::pnorm(-abs(MSE.z)))

  if(i==1){
    MSE<-matrix(c(MSE,MSE.se,MSE.z,P.MSE),nrow=1,ncol=4)
  }else{
    MSE<-cbind(MSE,MSE.se,MSE.z,P.MSE)
  }

  colnames(MSE)<-c(estimate,"Bootstrap SE","Z","P")

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
                bootstraps=bootMat)
    }

    return(MSE)

  }

  #yes, interaction
  SD<-diff(MSE[,1])
  se<-1:(loop_length-1)
  for(i in 1:length(se)){
    k<-i+1
    se[i]<-sd(bootMat[,i]-bootMat[,k],na.rm=T)
  }

  z<-SD/se
  P_SD<-2*(stats::pnorm(-abs(z)))

  SD_dat<-matrix(c(SD,se,z,P_SD),nrow=1,ncol=4)
  colnames(SD_dat)<-c(estimate,"Bootstrap SE","Z","P")
  rownames(SD_dat)<-paste(at.2[-c(length(at.2))],"to",at.2[-c(1)])

  if(isTRUE(return_Jac)){
    MSE<-list(change_in_MSE=SD_dat,
              MSE=MSE,
              bootstraps=bootMat)
  }else{
     MSE<-list(change_in_MSE=SD_dat,
            MSE=MSE)
  }
  return(MSE)



}
