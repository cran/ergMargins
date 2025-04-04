#' Function to compute average marginal effects in ERGM
#' If var2 and inter are left NULL, the function returns the average marginal effect for var 1.
#' if var2 and inter are specified, function conducts tests of second differences to assess significance of an interaction
#'
#' @param model is an ergm object
#' @param var1 is the name of the main effect, character string
#' @param var2 is the name of the moderator, character string
#' @param inter is the name of the interaction, character string
#' @param at.2 is a vector specifying the levels of var2 at which to compute the marginal effects. If left NULL, it computes the AME at all unique values of var2. Default is NULL.


#If the moderator is binary or specified at only 2 levels,
#the output object is a 2 dimensional list with one matrix of
#marginal effects and another of the second differences (ADC)

#If the moderator is continuous and specified at 3 or more levels,
#the output object contains a 3 dimensional list, where the Aggregate output is the
#average second difference and t tests based on the average second difference,
#the second differnece matrix is the matrix of second differneces between
#all adjacent levels of var 2, and the marginal matrix is the matrix of average marginal effects

#standard errors are computed using the delta method.


ergm.AME_count<-function(model,
                   var1,
                   var2=NULL,
                   inter=NULL,
                   at.2=NULL,
                   at.controls=NULL,
                   control_vals=NULL,
                   return.dydx=FALSE,
                   return.at.2=FALSE){


  reference<-model$reference

  if(grepl("Binomial", as.character(reference)[2],fixed = TRUE)){
    ref<-"binomial"
  }
  if("Poisson"%in%as.character(reference)){
    ref<-"poisson"
  }
  if("Geometric"%in%as.character(reference)){
    ref<-"geometric"
  }
  if("Uniform"%in%as.character(reference)){
    ref<-"uniform"
  }

  max_y<-max(as.sociomatrix(model$network,model$network$gal$ergm$response[1]))
  #get sufficient statistics
  suff_stats<-ergmCntPrep(formula=model$formula,
                  nw=model$network,
                  response=model$network$gal$ergm$response[1],
                  reference=ref,
                  max.count = max_y)


  out<-lapply(suff_stats$cs,as.data.frame)
  for(i in 1:length(out)){
    out[[i]]$m<-nrow(out[[i]])
  }
  dyad.mat<-dplyr::bind_rows(out)

  #get predicted probs then do m x p to get predicted counts
  theta<-btergm::coef(model)
  m<-dyad.mat$m
  dyad.mat<-dyad.mat[,-ncol(dyad.mat)]
  lp<-as.matrix(dyad.mat)%*%theta
  p<-m*(1/(1+exp(-lp)))

  vc<-stats::vcov(model)



  ##assign fixed values for controls when specified

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

    #get predicted counts
    lp<-as.matrix(dyad.mat)%*%theta
    p<-m*(1/(1+exp(-lp)))

  }


  if(nrow(dyad.mat)>1e06){
    message("There are over 1 million dyads in the ERGM sample space. Variance estimates for marginal effects may take a moment to compute.")
  }

  if(length(theta)>20){
    message("There are more than 20 parameters in the model. Variance estimates for marginal effects may take a moment to compute.")
  }

  ##identify unique values of at.2--not used if var2==NULL
  if(is.null(at.2)){
    at.2<-sort(unique(dyad.mat[,var2]))
  }

  if(length(at.2)>10){
    warning("More than 10 values of at.2 exist for the moderating variable. It may take awhile to compute average marginal effects. Consider specifying fewer values of at.2.")
  }


  ##marginal effects with no interaction
  if(is.null(var2)){

    if(length(var1)>1){
      AME.fun<-function(theta){

        ME.ergm<-sapply(names(theta),function(x)
          p*sum(theta[var1]))

        mean(ME.ergm,na.rm = TRUE)}
    }else{

    AME.fun<-function(theta){

      ME.ergm<-sapply(names(theta),function(x)
        p*theta[var1])

      mean(ME.ergm,na.rm = TRUE)}
    }

    AME<-AME.fun(theta)
    Jac<-numDeriv::jacobian(AME.fun,theta)
    variance.ame<-Jac%*%vc%*%t(Jac)

    AME.se<-sqrt(variance.ame)
    AME.z<-AME/AME.se
    P.AME<-2*(stats::pnorm(-abs(AME.z)))

    AME<-matrix(c(AME,AME.se,AME.z,P.AME),nrow=1,ncol=4)
    colnames(AME)<-c("AME","Delta SE","Z","P")
    rownames(AME)<-ifelse(length(var1)>1,paste(var1,collapse="+"),var1)
    AME<-signif(AME,digits=5)

    if(return.dydx==TRUE){
      dydx<-sapply(names(theta),function(x)
        p*theta[var1])
      AME<-list(AME,dydx[,var1],Jac)
      names(AME)<-c("AME","dydx","Jac")

    }

    return(AME)

  }else{


    ##marginal effects for interaction that does not vary with covariates
    if(!is.na(pmatch("nodematch",inter))){
      if(!is.na(pmatch("nodecov",var1)) | !is.na(pmatch("nodeicov",var1)) |
         !is.na(pmatch("nodeocov",var1))){
        ##matched nodal characteristics are not a product term for node covariates, so compute marginal effects
        #for var 1 and var 2, then use results to compute marignal effect for interaction

        AME.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            p*theta[var1])
          mean(ME.ergm,na.rm = TRUE)}

        AME1<-AME.fun(theta)
        Jac<-numDeriv::jacobian(AME.fun,theta)
        variance.ame1<-Jac%*%vc%*%t(Jac)
        AME1.se<-sqrt(variance.ame1)
        AME1.z<-AME1/AME1.se
        P.AME1<-2*(stats::pnorm(-abs(AME1.z)))

        AME.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            p*theta[var2])
          mean(ME.ergm,na.rm = TRUE)}

        AME2<-AME.fun(theta)
        Jac<-numDeriv::jacobian(AME.fun,theta)
        variance.ame2<-Jac%*%vc%*%t(Jac)
        AME2.se<-sqrt(variance.ame2)
        AME2.z<-AME2/AME2.se
        P.AME2<-2*(stats::pnorm(-abs(AME2.z)))

        AME<-matrix(c(AME1,AME1.se,AME1.z,P.AME1,
                      AME2,AME2.se,AME2.z,P.AME2),nrow=2,ncol=4,byrow=TRUE)
        colnames(AME)<-c("AME","Delta SE","Z","P")
        rownames(AME)<-c(var1,var2)
        marginal.matrix<-AME


        ##compute marginal effect
        AME.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            p*theta[inter])
          mean(ME.ergm,na.rm = TRUE)}

        AME<-AME.fun(theta)
        Jac<-numDeriv::jacobian(AME.fun,theta)
        variance.inter<-Jac%*%vc%*%t(Jac)
        AME.se<-sqrt(variance.inter)
        AME.z<-AME/AME.se
        P.AME<-2*(stats::pnorm(-abs(AME.z)))

        AME<-matrix(c(AME,AME.se,AME.z,P.AME),nrow=1,ncol=4)
        colnames(AME)<-c("AME","Delta SE","Z","P")
        rownames(AME)<-inter
        #  message("NOTE: Nodematch is an interaction, but it is not a product of the main effects (e.g., inter!=var1*var2). Returning the simple AME for the interaction. Consider respecifying ERGM using nodefactor for main effects or absdiff instead of nodematch to measure homophily.")
        marginal.matrix<-signif(marginal.matrix,digits=5)
        AME<-signif(AME,digits=5)
        AME<-list(AME,marginal.matrix)
        names(AME)<-c("Marginal effect for nodematch","Marginal effects for nodal covariates")
        return(AME)
      }

    }

    #for undirected networks, binarize factor variables
    if(!is.na(pmatch("nodefactor",var1))){
      dyad.mat[,var1][which(dyad.mat[,var1]>=2)]<-1
    }

    if(!is.na(pmatch("nodefactor",var2))){
      at.2<-c(0,1)
    }

    if(var1==var2){
      self.int<-TRUE
      var2<-paste(var1,".mod")
      dyad.mat[,var2]<-dyad.mat[,var1]
    }else{
      self.int<-FALSE
    }

    ##marginal effects for interactions
    marginal.matrix<-matrix(0,nrow=length(at.2),ncol=5)
    colnames(marginal.matrix)<-c("AME","Delta SE","Z","P","N")
    rownames(marginal.matrix)<-paste(var2,"==",at.2)
    dydx.list<-list()


    for(i in 1:nrow(marginal.matrix)){
      dyad.submat<-dyad.mat
      dyad.submat[,var2]<-at.2[i]

      if(class(model)%in%"mtergm"){
        if(ergm::is.curved(model@ergm)){
          theta<-ergm::ergm.eta(theta,model@ergm$etamap)
        }
      }else{
        if(ergm::is.curved(model)){
          theta<-ergm::ergm.eta(theta,model$etamap)
        }
      }
      #marginal effects for absolute differences
      if(!is.na(pmatch("absdiff",inter))){
        dyad.submat[,inter]<-abs(dyad.submat[,var1]-dyad.submat[,var2])
        if(i ==1){message(paste("Note that marginal effects for absolute differences are computed holding",var1,"at its mean. The mean for",var1,"is", mean(dyad.mat[,var1])))}
        if(self.int==TRUE){
          dyad.submat<-dyad.submat[,!colnames(dyad.submat)%in%var2]
        }
         p<-m*(1/(1+exp(-(apply(dyad.submat,1,function(x) t(x)%*%theta)))))


        if(class(model)%in%"mlergm"){
          theta<-model$theta
          vc<-solve(vc) #invert the fisher matrix
        }else{
          theta<-btergm::coef(model)
        }

        at.diffs<-abs(at.2[i]-mean(dyad.mat[,var1]))

        AME.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            (p*(theta[var1]+(theta[inter]*at.diffs))))
          mean(ME.ergm,na.rm = TRUE)}

        dydx.list[[i]]<-sapply(names(theta),function(x)
          (p*(theta[var1]+(theta[inter]*at.diffs))))


      }else{

        #marginal effects for product terms
        dyad.submat[,inter]<-dyad.submat[,var1]*dyad.submat[,var2]
        if(self.int==TRUE){
          dyad.submat<-dyad.submat[,!colnames(dyad.submat)%in%var2]
        }
        p<-m*(1/(1+exp(-(apply(dyad.submat,1,function(x) t(x)%*%theta)))))

        if(class(model)%in%"mlergm"){
          theta<-model$theta
          vc<-solve(vc) #invert the fisher matrix
        }else{
          theta<-btergm::coef(model)
        }

        AME.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            (p*(theta[var1]+(theta[inter]*at.2[i]))))
          mean(ME.ergm,na.rm = TRUE)}

        dydx.list[[i]]<-sapply(names(theta),function(x)
          (p*(theta[var1]+(theta[inter]*at.2[i]))))
      }

      AME<-AME.fun(theta)
      Jac<-numDeriv::jacobian(AME.fun,theta)
      marginal.matrix[i,1]<-AME

      if(i==1){
        Jac1<-matrix(Jac)
      }else{
        Jac1<-cbind(Jac1,matrix(Jac))}

    }


    variance.ame<-t(Jac1)%*%vc%*%Jac1
    AME.se<-sqrt(diag(variance.ame))

    marginal.matrix[,2]<-AME.se
    marginal.matrix[,3]<-marginal.matrix[,1]/AME.se
    marginal.matrix[,4]<-2*(stats::pnorm(-abs(marginal.matrix[,3])))
    marginal.matrix[,5]<-length(p)

    if(length(at.2)==1){
      AME<-t(as.matrix(marginal.matrix[,-c(5)]))
      rownames(AME)<-rownames(marginal.matrix)
      if(return.dydx==TRUE){
        AME<-list(AME,dydx.list,Jac1)
        names(AME)<-c("Average Marginal effects","Marginal effects","Jac")
      }
      if(return.at.2==TRUE){
        AME<-list(main.results=AME,
                  at.2=at.2)
      }
      return(AME)
    }

    second.diffs.mat<-as.matrix(diff(marginal.matrix)[,1:4])
    if(length(at.2)==2){
      second.diffs.mat<-t(second.diffs.mat)
    }


    for(j in 1:nrow(second.diffs.mat)){
      k<-j+1

      diff.se<-sqrt((marginal.matrix[j,2]^2)+(marginal.matrix[k,2]^2)-2*variance.ame[j,k])

      df<-marginal.matrix[j,5]-length(theta)
      z.ADC<-(second.diffs.mat[j,1])/diff.se
      P.ADC<-2*stats::pnorm(-abs(z.ADC))

      second.diffs.mat[j,2]<-diff.se
      second.diffs.mat[j,3]<-z.ADC
      second.diffs.mat[j,4]<-P.ADC

    }

    colnames(second.diffs.mat)<-c("Second. diff.","SE","Wald Z","P")
    rownames(second.diffs.mat)<-paste(at.2[-c(length(at.2))],"to",at.2[-c(1)])
    marginal.matrix<-signif(marginal.matrix,digits=5)
    second.diffs.mat<-signif(second.diffs.mat,digits=5)

    if(length(at.2)==2){

      if(return.dydx==FALSE){
        ADC<-list(second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))])
        names(ADC)<-c("Second differences","Average Marginal effects")}else{
          ADC<-list(second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))],dydx.list,Jac1)
          names(ADC)<-c("Second differences","Average Marginal effects","Marginal effects","Jac")
        }
      if(return.at.2==TRUE){
        ADC<-list(ADC,at.2)
        names(ADC)<-c("main.results","at.2")
      }

      return(ADC)
    }else{

      #use absolute t value in case of extreme negatives or positives
      summary.output<-matrix(c(mean(second.diffs.mat[,1]),mean(abs(second.diffs.mat[,3])),NA),nrow=1,ncol=3)
      colnames(summary.output)<-c("Mean Second diff.","Mean |Z|", "P")
      summary.output[1,3]<-2*stats::pnorm(abs(summary.output[1,2]),lower.tail = FALSE)
      summary.output<-signif(summary.output,digits=5)

      if(return.dydx==FALSE){
        ADC<-list(summary.output,second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))])
        names(ADC)<-c("Aggregate output","Second differences","Average Marginal effects")}else{

          ADC<-list(summary.output,second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))],dydx.list,Jac1)
          names(ADC)<-c("Aggregate output","Second differences","Average Marginal effects","Marginal effects","Jac")
        }
      if(return.at.2==TRUE){
        ADC<-list(ADC,at.2)
        names(ADC)<-c("main.results","at.2")
      }

      return(ADC)
    }
  }

}
