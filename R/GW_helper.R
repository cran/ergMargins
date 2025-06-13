#' A helper function for computing change statistics for geometrically weighted terms
#' returns a scalar representing the weighted value assigned to a pre-determined number of subgraph change statistics provided by the researcher
#'
#' @term_count is the subgraph count (e.g., 2 triangles/ edgewise shared partners)to be converted into geometrically weighted change statistics
#' @decay is the decay parameter may be specified apriori by the researcher in cases of fixed decay or estimated from the data in the case of curved ERGM
#' @lower_bound is the subgraph count prior to introducing a focal tie. Default value is 0



GW_helper<-function(term_count,decay,lower_bound=0){

  # start checks
  if(!is.numeric(term_count)){
    stop("term_count must be numeric. Cannot continue")
  }
  if(!is.numeric(decay)){
    stop("decay must be numeric. Cannot continue")
  }

  if(decay<0){
    stop("decay must be greater than 0 Cannot continue")
  }

  if(term_count<0 | lower_bound<0){
    stop("either term_count or lower_bound < 0. Cannot continue. ")
  }

  if(length(term_count)>1){
    term_count<-term_count[1]
    warning("GW_helper currently only accepts scalars. Only the first entry in term_count used for computation.")

  }

  if(length(decay)>1){
    decay<-decay[1]
    warning("Multiple entries provided for decay. Only the first entry being used.")

  }

  if(lower_bound>=term_count){
    stop("lower_bound meets or exceeds term_count. Please respecify")
  }

  #start calculation
  if(term_count==0){
    return(0)
  }

  if(lower_bound==0){
    lower_bound<-1
  }


  out_vec<-matrix(lower_bound:term_count-1,ncol=1)
  out_vec<-unique(out_vec) #handle lower_bound = 0 and upper_bound = 1

  gw_decay<-function(x){
    (1-exp(-decay))^x
  }

  out_vec<-apply(out_vec,1,gw_decay)

  out_val<-sum(out_vec)
  return(out_val)


}



