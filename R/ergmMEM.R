#' Function to compute marginal effects at means in ERGM
#' If var2 and inter are left NULL, the function returns the marginal effect for var 1.
#' if var2 and inter are specified, function conducts tests of second differences to assess significance of an interaction
#'
#' @param model is an ergm object
#' @param var1 is the name of the main effect, character string
#' @param var2 is the name of the moderator, character string
#' @param inter is the name of the interaction, character string
#' @param at.2 is a vector specifying the levels of var2 at which to compute the marginal effects. If left NULL, it computes the MEM at all unique values of var2. Default is NULL.
#If the moderator is binary or specified at only 2 levels,
#the output object is a 2 dimensional list with one matrix of
#marginal effects and another of the second differences

#If the moderator is continuous and specified at 3 or more levels,
#the output object contains a 3 dimensional list, where the Aggregate output is the
#average second difference and t tests based on the average second difference,
#the second differnece matrix is the matrix of second differneces between
#all adjacent levels of var 2, and the marginal matrix is the matrix of  marginal effects

#standard errors are computed using the delta method.


ergm.MEM<-function(model,
                   var1,
                   var2=NULL,
                   inter=NULL,
                   at.2=NULL,
                   at.controls=NULL,
                   control_vals=NULL,
                   return.dydx=FALSE){
  if(length(var1)>1&&!is.null(var2)){
    stop("Joint tests not allowed for interactions. Please respecify.")
  }
  if(length(var1)>1&&return.dydx==TRUE){
    message("dydx not avaiable for joint tests")
    return.dydx<-FALSE
  }
  if(class(model)[1]%in%"btergm"){
    out<-ergm.MEM_boot(model=model,
                       var1=var1,
                       var2=var2,
                       inter=inter,
                       at.2=at.2,
                       at.controls=at.controls,
                       control_vals=control_vals,
                       return.dydx=return.dydx)
    return(out)
  }

  if(!class(model)[1]%in%"mtergm"){
    if(is.valued(model)){
      out<-ergm.MEM_count(model=model,
                         var1=var1,
                        var2=var2,
                        inter=inter,
                        at.2=at.2,
                        at.controls=at.controls,
                        control_vals=control_vals,
                        return.dydx=return.dydx)

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
    dyad.mat<-ergmMPLE(model$formula,output="dyadlist",basis=model$network)
    dyad.mat<-dyad.mat$predictor[,-c(1:2)]
    theta<-coef(model)
    vc <- stats::vcov(model)

  }

  if(class(model)[1]%in%"mlergm"){
    theta<-model$theta
    vc<-solve(vc)
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
  ##handle decay term in curved ergms
  if(class(model)[1]%in%"mtergm" | class(model)[1]%in%"btergm"){

    if(ergm::is.curved(model@ergm)){
      curved.term<-vector(length=length(model@ergm$etamap$curved))
      for(i in 1:length(model@ergm$etamap$curved)){
        curved.term[i]<-model@ergm$etamap$curved[[i]]$from[2]
      }
    #  theta<-theta[-c(curved.term)]
      theta<-ergm::ergm.eta(theta,model@ergm$etamap)
      theta<-theta[-c(length(theta))]
      names(theta)<-colnames(dyad.mat)

    }

  }else{
    if(ergm::is.curved(model)){
      curved.term<-vector(length=length(model$etamap$curved))
      for(i in 1:length(model$etamap$curved)){
        curved.term[i]<-model$etamap$curved[[i]]$from[2]
      }
      #theta<-theta[-c(curved.term)]
      theta<-ergm::ergm.eta(theta,model$etamap)
      names(theta)<-colnames(dyad.mat)



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

  if(any(names(theta)!=colnames(dyad.mat))){
    colnames(dyad.mat)<-names(theta) #make sure names align
  }


  ###incorporate control vals

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




  #create marginal effects
  dyad.means<-colMeans(dyad.mat,na.rm=TRUE)

  p<-1/(1+exp(-dyad.means%*%theta))

  if (class(model)[1] %in% "mlergm") {
    theta <- model$theta
    vc <- solve(vc)
  }
  else {
    theta <- btergm::coef(model)
    if(!is.na(offset_ind)){
      theta<-theta[-offset_ind]
    }
  }

  ##identify unique values of at.2--not used if var2==NULL
  if(is.null(at.2)){
    at.2<-sort(unique(dyad.mat[,var2]))
  }


  ##marginal effects with no interaction
  if(is.null(var2)){

    if(length(var1)>1){
      MEM.fun<-function(theta){

        ME.ergm<-sapply(names(theta),function(x)
          (p*(1-p)*sum(theta[var1])))
        mean(ME.ergm,na.rm = TRUE)}
    }else{
       MEM.fun<-function(theta){

        ME.ergm<-sapply(names(theta),function(x)
        (p*(1-p)*theta[var1]))
        mean(ME.ergm,na.rm = TRUE)}
    }
    MEM<-MEM.fun(theta)
    Jac<-numDeriv::jacobian(MEM.fun,theta)
    variance.mem<-Jac%*%vc%*%t(Jac)


    MEM.se<-sqrt(variance.mem)
    MEM.z<-MEM/MEM.se
    P.MEM<-2*(stats::pnorm(-abs(MEM.z)))

    MEM<-matrix(c(MEM,MEM.se,MEM.z,P.MEM),nrow=1,ncol=4)
    colnames(MEM)<-c("MEM","Delta SE","Z","P")
    rownames(MEM)<-ifelse(length(var1)==1,var1,paste(var1,collapse="+"))
    MEM<-signif(MEM,digits=5)

    if(return.dydx==TRUE){
      MEM<-list(MEM,Jac)
      names(MEM)<-c("MEM","Jac")

    }


    return(MEM)

  }else{


    ##marginal effects for interaction that does not vary with covariates
    if(!is.na(pmatch("nodematch",inter))){
      if(!is.na(pmatch("nodecov",var1)) | !is.na(pmatch("nodeicov",var1)) |
         !is.na(pmatch("nodeocov",var1))){          ##matched nodal characteristics are not a product term, so compute marginal effects
          #for var 1 and var 2, then use results to compute marignal effect for interaction

          MEM.fun<-function(theta){

            ME.ergm<-sapply(names(theta),function(x)
              (p*(1-p)*theta[var1]))
            mean(ME.ergm,na.rm = TRUE)}

          MEM1<-MEM.fun(theta)
          Jac<-numDeriv::jacobian(MEM.fun,theta)
          variance.ame1<-Jac%*%vc%*%t(Jac)
          MEM1.se<-sqrt(variance.ame1)
          MEM1.z<-MEM1/MEM1.se
          P.MEM1<-2*(stats::pnorm(-abs(MEM1.z)))

          MEM.fun<-function(theta){

            ME.ergm<-sapply(names(theta),function(x)
              (p*(1-p)*theta[var2]))
            mean(ME.ergm,na.rm = TRUE)}

          MEM2<-MEM.fun(theta)
          Jac<-numDeriv::jacobian(MEM.fun,theta)
          variance.ame2<-Jac%*%vc%*%t(Jac)
          MEM2.se<-sqrt(variance.ame2)
          MEM2.z<-MEM2/MEM2.se
          P.MEM2<-2*(stats::pnorm(-abs(MEM2.z)))

          MEM<-matrix(c(MEM1,MEM1.se,MEM1.z,P.MEM1,
                        MEM2,MEM2.se,MEM2.z,P.MEM2),nrow=2,ncol=4,byrow=TRUE)
          colnames(MEM)<-c("MEM","Delta SE","Z","P")
          rownames(MEM)<-c(var1,var2)
          marginal.matrix<-MEM


          ##compute marginal effect
          MEM.fun<-function(theta){

            ME.ergm<-sapply(names(theta),function(x)
              (p*(1-p)*theta[inter]))
            mean(ME.ergm,na.rm = TRUE)}

          MEM<-MEM.fun(theta)
          Jac<-numDeriv::jacobian(MEM.fun,theta)
          variance.inter<-Jac%*%vc%*%t(Jac)
          MEM.se<-sqrt(variance.inter)
          MEM.z<-MEM/MEM.se
          P.MEM<-2*(stats::pnorm(-abs(MEM.z)))

          MEM<-matrix(c(MEM,MEM.se,MEM.z,P.MEM),nrow=1,ncol=4)
          colnames(MEM)<-c("MEM","Delta SE","Z","P")
          rownames(MEM)<-inter
          message("NOTE: Nodematch is an interaction, but it is not a product of the main effects (e.g., inter!=var1*var2). Returning the simple MEM for the interaction. Consider respecifying ERGM using nodefactor for main effects or absdiff instead of nodematch to measure homophily.")
          marginal.matrix<-signif(marginal.matrix,digits=5)
          MEM<-signif(MEM,digits=5)
          MEM<-list(MEM,marginal.matrix)
          names(MEM)<-c("Marginal effect for nodematch","Marginal effects for nodal covariates")

          if(return.dydx==TRUE){
            MEM<-list(MEM,Jac)
            names(MEM)<-c("MEM","Jac")

          }

          return(MEM)
        }

    }

    #for undirected networks, binarize factor variables
    if(!is.na(pmatch("nodefactor",var1))){
      dyad.mat[,var1][which(dyad.mat[,var1]>=2)]<-1
    }

    if(!is.na(pmatch("nodefactor",var2))){
      at.2<-c(0,1)
    }

      #check whether self interaction
    if(var1==var2){
      self.int<-TRUE
      var2<-paste(var1,".mod",sep="")
      dyad.means[var2]<-dyad.means[var1]
    }else{
      self.int<-FALSE
    }

    ##marginal effects for interactions
    marginal.matrix<-matrix(NA,nrow=length(at.2),ncol=5)
    colnames(marginal.matrix)<-c("MEM","Delta SE","Z","P","N")
    rownames(marginal.matrix)<-paste(var2,"==",at.2)


    for(i in 1:nrow(marginal.matrix)){

      dyad.submeans<-dyad.means
      dyad.submeans[var2]<-at.2[i]
      if(class(model)[1]%in%"mtergm"){
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
        dyad.submeans[inter]<-abs(dyad.submeans[var1]-dyad.submeans[var2])

        if(self.int==TRUE){
          dyad.submeans<-dyad.submeans[!names(dyad.submeans)%in%var2]
        }

        p<-1/(1+exp(-(dyad.submeans%*%theta)))

        at.diffs<-abs(at.2[i]-dyad.submeans[var1])
        if(class(model)[1]%in%"mlergm"){
          theta<-model$theta
          vc<-solve(vc) #invert the fisher matrix
        }else{
          theta<-btergm::coef(model)
        }

        MEM.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            (p*(1-p)*(theta[var1]+(theta[inter]*at.diffs))))
          mean(ME.ergm,na.rm = TRUE)}

      }else{

        #marginal effects for product terms
        dyad.submeans[inter]<-dyad.submeans[var1]*dyad.submeans[var2]
        if(self.int==TRUE){
          dyad.submeans<-dyad.submeans[!names(dyad.submeans)%in%var2]
        }

        p<-1/(1+exp(-(dyad.submeans%*%theta)))
        if(class(model)[1]%in%"mlergm"){
          theta<-model$theta
          vc<-solve(vc) #invert the fisher matrix
        }else{
          theta<-btergm::coef(model)
        }

        MEM.fun<-function(theta){

          ME.ergm<-sapply(names(theta),function(x)
            (p*(1-p)*(theta[var1]+(theta[inter]*at.2[[i]]))))
          mean(ME.ergm,na.rm = TRUE)}
      }

      MEM<-MEM.fun(theta)
      Jac<-numDeriv::jacobian(MEM.fun,theta)
      marginal.matrix[i,1]<-MEM

      if(i==1){
        Jac1<-matrix(Jac)
      }else{
        Jac1<-cbind(Jac1,matrix(Jac))}

    }


    variance.mem<-t(Jac1)%*%vc%*%Jac1
    MEM.se<-sqrt(diag(variance.mem))

    marginal.matrix[,2]<-MEM.se
    marginal.matrix[,3]<-marginal.matrix[,1]/MEM.se
    marginal.matrix[,4]<-2*(stats::pnorm(-abs(marginal.matrix[,3])))
    marginal.matrix[,5]<-length(p)


    if(length(at.2)==1){
      MEM<-t(as.matrix(marginal.matrix[,-c(5)]))
      rownames(MEM)<-rownames(marginal.matrix)
      if(return.dydx==TRUE){
        MEM<-list(MEM,Jac1)
        names(MEM)<-c("MEM","Jac")

      }

      return(MEM)
    }


    second.diffs.mat<-as.matrix(diff(marginal.matrix)[,1:4])
    if(length(at.2)==2){
      second.diffs.mat<-t(second.diffs.mat)
    }


    for(j in 1:nrow(second.diffs.mat)){
      k<-j+1

      diff.se<-sqrt((marginal.matrix[j,2]^2)+(marginal.matrix[k,2]^2)-2*variance.mem[j,k])

      df<-marginal.matrix[j,5]-length(theta)
      z.DCR<-(second.diffs.mat[j,1])/diff.se
      P.DCR<-2*stats::pnorm(-abs(z.DCR))

      second.diffs.mat[j,2]<-diff.se
      second.diffs.mat[j,3]<-z.DCR
      second.diffs.mat[j,4]<-P.DCR

    }

    colnames(second.diffs.mat)<-c("Second diff","SE","Wald Z","P")
    rownames(second.diffs.mat)<-paste(at.2[-c(length(at.2))],"to",at.2[-c(1)])
    marginal.matrix<-signif(marginal.matrix,digits=5)
    second.diffs.mat<-signif(second.diffs.mat,digits=5)


    if(length(at.2)==2){
      if(return.dydx==TRUE){

        DCR<-list(second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))],Jac1)
        names(DCR)<-c("Second differences","Marginal effects at means","Jac")

      }else{

         DCR<-list(second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))])
         names(DCR)<-c("Second differences","Marginal effects at means")
      }
      return(DCR)
    }else{

      #use absolute t value in case of extreme negatives or positives
      summary.output<-matrix(c(mean(second.diffs.mat[,1]),mean(abs(second.diffs.mat[,3])),NA),nrow=1,ncol=3)
      colnames(summary.output)<-c("Mean Second diff.","Mean |Z|", "P")
      summary.output[1,3]<-2*stats::pnorm(abs(summary.output[1,2]),lower.tail = FALSE)
      summary.output<-signif(summary.output,digits=5)

      if(return.dydx==TRUE){

        DCR<-list(summary.output,second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))],Jac1)
        names(DCR)<-c("Aggregate output","Second differences","Marginal effects at means","Jac")

      }else{

        DCR<-list(summary.output,second.diffs.mat,marginal.matrix[,-c(ncol(marginal.matrix))])
        names(DCR)<-c("Aggregate output","Second differences","Marginal effects at means")
      }
      return(DCR)
    }
  }

}
