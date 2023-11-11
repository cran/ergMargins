ergmCntPrep<-function(formula,nw,response,reference=c("uniform", "binomial", "geometric", "poisson"), max.count=Inf, max.count.safety=4, max.count.edgewise=TRUE, count.samples=Inf, sample.size=Inf, sample.method=c("weighted", "random"), weight.type=c("TNT", "flatval"), cores=1, seed=NULL, must.count=1,WtSumAsSampSiz=T){
  #Get the network, and set things up
  nw<-ergm.getnetwork(formula)
  isdir<-is.directed(nw)
  n<-network.size(nw)
  nev<-network.dyadcount(nw)              #Number of edge variables
  ss<-min(sample.size,nev)                #Sample size
  lrmr<-switch(match.arg(reference),      #Function to calculate log ref mes rat
               uniform = function(k,y) {rep(0,length(k))},
               binomial = function(k,y) {lchoose(max.count,k)-lchoose(max.count,y)},
               geometric = function(k,y) {rep(0,length(k))},
               poisson = function(k,y) {lfactorial(y)-lfactorial(k)}
  )
  #Set up edge variable sampling - use inverse freq weighting to smooth a bit
  ally<-as.sociomatrix(nw,response)       #All y values
  alli<-row(ally)                         #Senders
  allj<-col(ally)                         #Receivers
  ally<-gvectorize(ally,mode=ifelse(isdir,"digraph","graph"),censor.as.na=FALSE)
  alli<-gvectorize(alli,mode=ifelse(isdir,"digraph","graph"),censor.as.na=FALSE)
  allj<-gvectorize(allj,mode=ifelse(isdir,"digraph","graph"),censor.as.na=FALSE)
  if(ss==nev){                            #If ss==nev, no sampling
    samp<-1:nev                              #Everyone's in the sample
    iw<-rep(1,nev)                           #Inclusion weight is 1
  }else{                                  #Else, sample EVs
    if(match.arg(sample.method)=="random"){
      set.seed(seed)                         #Set seed for sampling
      samp<-sample(1:nev,ss)                 #Choose at random
      iw <- rep(1,ss)                        #Inclusion weight, before regularization
    }else{
      if(match.arg(weight.type)=="TNT"){          #"Tie/No-Tie" style weighting
        taby<-table(ally>0)  #c(FALSE, TRUE)
        wght<-(1/(2*taby))[1+(ally>0)]
        wght<-inclusionprobabilities(wght,ss)
        set.seed(seed)                             #Set seed for sampling
        samp<-which(UPpoisson(wght)>0)             #Draw the sample
        iw<- 1/(wght[samp])                        #Inclusion weight, before regularization
      }else if(match.arg(weight.type)=="flatval"){ #"Flat" value distribution weighting
        if(ss/nev>0.15)
          warning("Target sample size is ",round(ss/nev*100),"% of the total EV count.  Weighted sampling may be unreliable here - you may want to consider random sampling.")
        taby<-table(ally)
        wght<-(1/length(taby)/taby)[match(ally,names(taby))] #Ideal weights
        wght<-inclusionprobabilities(wght,ss)
        set.seed(seed)                           #Set seed for sampling
        samp<-which(UPpoisson(wght)>0)           #Draw the sample
        iw<- 1/(wght[samp])                        #Inclusion weight, before regularization
      }else{
        stop("Unknown weighting method ",weight.type," in ergmCntPrep.  Cannot go on like this.\n")
      }
    }
  }
  ss<-length(samp)                        #Should not have changed, but can...
  #Inclusion weight: regularize it to sample size or total dyads.
  if(WtSumAsSampSiz){
    iw <- iw*ss/sum(iw)
  }else{
    iw <- iw*nev/sum(iw)
  }
  #There may be some ignoreable numerical difference
  if(WtSumAsSampSiz & sum(iw)!=ss) warning("WtSumAsSampSiz=T, but sum of inclusion weight unequal to sample size, diff=",sum(iw)-ss)
  if(WtSumAsSampSiz==F & sum(iw)!=nev) warning("WtSumAsSampSiz=F, but sum of inclusion weight unequal to N of edge variables, diff=",sum(iw)-ss)

  y<-ally[samp]                           #Observed sampled y values
  snd<-alli[samp]                         #Sampled senders
  rec<-allj[samp]                         #Sampled receivers
  #Walk through edge variables and calculate exciting things
  if(max.count.edgewise){        #Use per-edge values
    yub<-pmin(pmax(must.count+1,ceiling(y+max.count.safety*4*sqrt(y))),max.count)  #Upper bounds
    ylb<-pmax(must.count+1,floor(y-max.count.safety*4*sqrt(y)))    #Lower bounds (0:must.count always included)
    yrng<-vector(mode="list",length=ss)        #y values for each edge
    ycwt<-vector(mode="list",length=ss)        #Weights for approximation
    for(i in 1:ss){
      if(count.samples>=yub[i]-ylb[i]){
        yrng[[i]]<-c(0:must.count,ylb[i]:yub[i])
        ycwt[[i]]<-rep(1,length(yrng[[i]]))
      }else{
        yrng[[i]]<-c(0:must.count,round(seq(from=ylb[i],to=yub[i],length=count.samples)))
        ycwt[[i]]<-c(rep(1, must.count+1),1,diff(round(seq(from=ylb[i],to=yub[i],length=count.samples))))
      }
    }
  }else{                         #Use uniform edge value ranges
    maxy<-max(ally)                         #Maximum y value
    yrng<-vector(mode="list",length=ss)     #y values for each edge
    ycwt<-vector(mode="list",length=ss)     #Weights for approximation
    if(max.count<Inf){                      #Range of y values for normalization
      if(count.samples>=max.count+1){
        for(i in 1:ss){
          yrng[[i]]<-0:max.count
          ycwt[[i]]<-rep(1,1+max.count)
        }
      }else{
        for(i in 1:ss){
          yrng[[i]]<-c(0,round(seq(from=1,to=max.count,length=count.samples)))
          ycwt[[i]]<-c(1,1,diff(yrng[[i]]))
        }
      }
    }else{
      if(count.samples>=max.count.safety*maxy+1){
        for(i in 1:ss){
          yrng[[i]]<-0:(max.count.safety*maxy)
          ycwt[[i]]<-rep(1,1+max.count.safety*maxy)
        }
      }else{
        for(i in 1:ss){
          yrng[[i]]<-round(seq(from=0,to=ceiling(max.count.safety*maxy),length=count.samples))
          ycwt[[i]]<-c(diff(yrng[[i]]),1)
        }
      }
    }
  }
  rmr<-parallel::mclapply(1:ss,function(i){         #Ref meas ratio (one entry per EV)
    lrmr(yrng[[i]],y[i])
  },mc.cores=cores)
  cs<-parallel::mclapply(1:ss,function(i){          #Changescore list (one entry per EV)
    ergm.godfather(formula=formula, changes=lapply(yrng[[i]],function(z){matrix(c(snd[i],rec[i],z),nrow=1)}), response=response, changes.only=TRUE)
  },mc.cores=cores)
  #Return a list with all the goodies
  out<-list(snd=snd,rec=rec,y=y,rmr=rmr,cs=cs,ycwt=ycwt,iw=iw)
  class(out)<-"ERGMCntPrep"
  out
}
