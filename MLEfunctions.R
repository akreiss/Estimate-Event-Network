## Main Function which computes the MLE.
#### Input ####
## mu_init    - start value for the solver, it is important that the dimension
##              is correct, i.e., its length must equal the number of parameters
## t_0        - Vector of all time points at which the estimator shall be computed
## h          - Bandwidth
## events     - List of events, see pdf for explanation
## covariates - List of covariates, see pdf for explanation
## T          - Observation horizon (observations are on the interval [0,T])
## K          - A function specifying the kernel to be used: Below you find the
##                function K which is the triangular kernel. The kernel must be supported on [-1,1],
##                integrate to one, have int_{-1}^1 uK(u)du=0 and be non-negative.
## KIntegrate - Function of two variables: KIntegrate(a,b) equals the integral over K
##                 from a to b.
## par        - If TRUE the foreach framework is used to perform parallel computations if available,
##                if it is not available then a warning is returned
## iterlim    - Number of interations after which the Newton Method is forced to end (cf. ?nlm)
## gradtol    - Another parameter for the Newton Method, see ?nlm for further information
##
#### Output ####
## A list of three elements:
## estimates - Matrix, each row corresponds to an estimate, row k corresponds to the estimate at
##                time t_0[k].
## hessians  - List containg the hessian matrices at the optimal points, again list element k is
##              the hessian corresponding to time t_0[k] and estimate[k,]
## t_0       - Just equals the input variable t_0 for convenience
## codes     - Vector of exit codes of the corresponding optimisation, entry k corresponds to t_0[k],
##               see ?nlm for detailed explanation, in short: 1=very good, 2=ok, 3=not so good,
##               4=worse (an increase of iterlim could fix this), 5=bad
MLE <- function(mu_init,t_0,h,events,covariates,T,K,KIntegrate,par=TRUE,iterlim=100,gradtol=0.000001) {
  ## This function is just a wrapper.
  ## MLE_raw is doing the actual optimisation
  if(par==TRUE) {
    MLEs <- foreach(t=t_0, .export = c("MLE_raw","nlikelihood","likelihood","int1","summing_function1","expo","summing_function2","summing_function3","integrator","xexpo","x2expo"), .packages = c("Matrix")) %dopar% MLE_raw(t_0=t,mu_init=mu_init,h=h,events=events,covariates=covariates,T=T,K=K,KIntegrate=KIntegrate,iterlim=iterlim,gradtol=gradtol)
  } else {
    MLEs <- lapply(t_0,MLE_raw,mu_init=mu_init,h=h,events=events,covariates=covariates,T=T,K=K,KIntegrate=KIntegrate,iterlim=iterlim,gradtol=gradtol)
  }
  
  ## Rearange Output
  N <- length(MLEs)
  q <- length(MLEs[[1]]$estimate)
  
  estimates <- matrix(0,ncol=q,nrow=N)
  hessians <- lapply(1:N,function(x) return(matrix(0,nrow=q,ncol=q)))
  codes <- 1:N
  
  for(i in 1:N) {
    estimates[i,] <- MLEs[[i]]$estimate
    hessians[[i]] <- MLEs[[i]]$hessian
    codes[i] <- MLEs[[i]]$code
  }
  
  return(list(estimates=estimates,hessians=hessians,t_0=t_0,codes=codes))
}

## Computes the MLE in the locally linear estimation procedure. The syntax of this function
## is identical to that of MLE, see the explanations there. Some options of MLE are not
## available for LL_MLE. Particularly, there is no possibility to chose the kernel function,
## it is always the triangular kernel. Note that in the local linear model the number of
## parameters equals twice the dimension of the covariates. The vector "mu_init" is required
## to have this dimension.
#### Input ####
## Up to "mu_init" which has dimension=2*(dimension of covariates) all input variables have
## the same interpretation as for MLE.
#### Output ####
## A list with elements "estimates", "hessians", "t_0", "codes" and "values".
## The first four have the same meaning as for the function MLE. Note that
## the values in "estimates" have dimension equal to the dimension of "mu_init",
## see also the remarks in the description above. The first half of estimates
## provides the local linear estimators. The element "values" contains
## the estimated optimal value of the likelihood.
LL_MLE <- function(mu_init,t_0,h,events,covariates,T,par=FALSE,iterlim=100,gradtol=0.000001) {
  ## This function is just a wrapper
  ## LL_MLE_raw is doing the optimisation
  if(par==TRUE) {
    MLEs <- foreach(t=t_0, .export = c("LL_MLE_raw","nlikelihood","likelihood","int1","summing_function1","expo","summing_function2","summing_function3","integrator","xexpo","x2expo"), .packages = c("Matrix")) %dopar% LL_MLE_raw(t_0=t,mu_init=mu_init,h=h,events=events,covariates=covariates,T=T,iterlim=iterlim,gradtol=gradtol)
    
  } else {
    MLEs <- lapply(t_0,LL_MLE_raw,mu_init=mu_init,h=h,events=events,covariates=covariates,T=T,iterlim=iterlim,gradtol=gradtol)
  }
  
  N <- length(MLEs)
  q <- length(MLEs[[1]]$estimate)
  
  estimates <- matrix(0,ncol=q,nrow=N)
  hessians <- lapply(1:N,function(x) return(matrix(0,nrow=q,ncol=q)))
  codes <- 1:N
  value <- 1:N
  
  for(i in 1:N) {
    estimates[i,] <- MLEs[[i]]$estimate
    hessians[[i]] <- MLEs[[i]]$hessian
    codes[i] <- MLEs[[i]]$code
    value[i] <- MLEs[[i]]$minimum
  }
  
  return(list(estimates=estimates,hessians=hessians,t_0=t_0,codes=codes,values=value))
}

## Computes a parametric estimator in the model, i.e., it is supposed that the parameter 
## function is a constant and the constant is estimated.
#### Input ####
## The options for this function have the same meaning as for MLE.
#### Output ####
## A list with the same elements as for MLE. But each element has only one
## entry which corresponds to the one estimate.
parametric_MLE <- function(mu_init,events,covariates,T,par=FALSE,iterlim=100,gradtol=0.000001) {
  res <- MLE(mu_init,T/2,T/2,events,covariates,T,K=Kuniform,KIntegrate = KuniformIntegrate,par=par,iterlim=iterlim,gradtol=gradtol)
  return(list(estimate=res$estimates[1,],hessian=res$hessians[[1]],code=res$codes[[1]]))
}

## Function which computes likelihood at t0, the meaning of the input variables
## is the same as for MLE, but t_0 cannot be a vector (it is just one number)
likelihood <- function(mu,t_0,h,events,covariates,T,K,KIntegrate) {
  ################################################
  ## Compute the first integral and derivatives ##
  ################################################
  relevant_indices <- which(events[,1]>=max(0,t_0-h) & events[,1]<=min(T,t_0+h))
  covar_indices <- events[relevant_indices,5]
  
  A <- mapply(int1,lapply(relevant_indices,function(i) events[i,]),covariates$covars[covar_indices],covariates$index[covar_indices],MoreArgs = list(t_0=t_0,h=h,K=K),SIMPLIFY=FALSE)
  
  integral1_derivative <- Reduce('+',A)
  integral1 <- integral1_derivative%*%mu
  integral1_hessian <- 0
  
  #################################################
  ## Compute the second integral and derivatives ##
  #################################################
  covar_timestamps <- sort(unique(covar_indices))
  if(covar_timestamps[1]>1) {covar_timestamps <- c(covar_timestamps[1]-1,covar_timestamps)}
  
  A <- lapply(covariates$covars[covar_timestamps],summing_function1,mu=mu)
  B <- lapply(covariates$covars[covar_timestamps],summing_function2,mu=mu)
  C <- lapply(covariates$covars[covar_timestamps],summing_function3,mu=mu)
  
  ct <- length(covar_timestamps)
  if(ct>1) {
    integral2 <- Reduce('+',mapply(integrator,covariates$time[covar_timestamps[1:(ct-1)]],covariates$time[covar_timestamps[2:ct]],A[1:(ct-1)],MoreArgs = list(KIntegrate=KIntegrate,T=T,t_0=t_0,h=h),SIMPLIFY=FALSE))
    integral2_derivative <- Reduce('+',mapply(integrator,covariates$time[covar_timestamps[1:(ct-1)]],covariates$time[covar_timestamps[2:ct]],B[1:(ct-1)],MoreArgs = list(KIntegrate=KIntegrate,T=T,t_0=t_0,h=h),SIMPLIFY=FALSE))
    integral2_hessian <- Reduce('+',mapply(integrator,covariates$time[covar_timestamps[1:(ct-1)]],covariates$time[covar_timestamps[2:ct]],C[1:(ct-1)],MoreArgs = list(KIntegrate=KIntegrate,T=T,t_0=t_0,h=h),SIMPLIFY=FALSE))
  } else {
    integral2 <- 0
    integral2_derivative <- 0
    integral2_hessian <- 0
  }
  integral2 <- integral2+integrator(covariates$time[covar_timestamps[ct]],min(t_0+h,T),A[[ct]],KIntegrate,T,t_0,h)
  integral2_derivative <- integral2_derivative+integrator(covariates$time[covar_timestamps[ct]],min(t_0+h,T),B[[ct]],KIntegrate,T,t_0,h)
  integral2_hessian <- integral2_hessian+integrator(covariates$time[covar_timestamps[ct]],min(t_0+h,T),C[[ct]],KIntegrate,T,t_0,h)
  
  ## Prepare Output
  output <- integral1-integral2
  attr(output,"gradient") <- integral1_derivative-integral2_derivative
  attr(output,"hessian") <- integral1_hessian-integral2_hessian
  
  return(output)
}


## Examples for kernel functions
K <- function(x) {
  out <- rep(0,length(x))
  indis <- which(abs(x)<=1)
  out[indis] <- 1-abs(x[indis])
  
  return(out)
}
KIntegrate_raw <- function(x) {
  if(x<=-1) {
    return(0)
  } else if (x<=0) {
    return(0.5*x^2+x+0.5)
  } else if (x<=1) {
    return(0.5-0.5*x^2+x)
  } else {
    return(1)
  }
}

KIntegrate <- function(a,b) {
  return(KIntegrate_raw(b)-KIntegrate_raw(a))
}

Kuniform <- function(x) {
  if(abs(x)>0.5) {return(0)} else {return(1)}
}
KuniformIntegrate <- function(a,b) {
  return(min(c(0.5,b))-max(c(-0.5,a)))
}

## This function adds a fifth column to the events matrix which contains the index
## of the covariate index which is active at time of the event.
#### Input ####
## "events" and "covariates" are both as in MLE.
#### Output ####
## A matrix which first four columns are identical to the input "events". Its fifth
## column contains the index of the covariate entry which is active at the time of
## the corresponding event.
add_covariate_index <- function(events,covariates) {
  events <- cbind(events,0)
  
  for(i in 1:(length(covariates$time)-1)) {
    indis <- which(events[,1]>=covariates$time[i] & events[,1]<covariates$time[i+1])
    events[indis,5] <- i
  }
  i <- length(covariates$time)
  indis <- which(events[,1]>=covariates$time[i])
  events[indis,5] <- i
  
  return(events)
}

## Computes estimates for the number of events in the network given an
## estimated model.
#### Input ####
## time_span  - Vector of increasing time points. Estimates are computed for
##              the intervals [time_span[1],time_span[2]], [time_span[3],time_span[4]],...
## parameter  - Estimated model, this must be an output from LL_MLE
## covariates - Covariate list as explained in MLE
## size       - Number of vertices in the network
## par        - Set to TRUE if cluster for foreach parallel computations is available.
#### Output ####
## ee           - List of "size"x"size" matrices: Each matrix in this list corresponds to one interval in time_span:
##                The first element to [time_span[1],time_span[2]], the second element to [time_span[3],time_span[4]],...
##                Each matrix contains the estimated number of events in the corresponding time interval and on
##                the corresponding edge.
## index_matrix - Just a copy of covariates$index
## time_span    - Just the input vector time_span
event_estimates <- function(time_span,parameter,covariates,size,par=FALSE) {
  ## Create list of intervals
  time_list <- lapply(1:(length(time_span)/2),function(i) return(c(time_span[2*(i-1)+1],time_span[2*(i-1)+2])))
  
  if(par==TRUE){
    ee <- foreach(t=time_list,.export=c("event_estimates_internal","single_estimate"),.packages=c("Matrix")) %dopar% event_estimates_internal(t=t,parameter=parameter,covariates=covariates,size=size)
  } else {
    ee <- lapply(time_list,event_estimates_internal,parameter=parameter,covariates=covariates,size=size)
  }
  
  return(list(ee=ee,index_matrices=covariates$index,time_span=time_span))
}

## This function compares the estimated number of events to the actually
## observed number of events. For each graph in "reality_graphs" the squared
## difference between estimated and actual events is computed (that is, sum
## of squares of differences on each edge).
#### Input ####
## ee             - Number of estimated events, must be output from event_estimates
## reality_graphs - List of graphs (igraph objects). Each graph encodes the 
##                  number of observed events per edge (in the edge attribute "count").
##                  This list must have the same length as ee and the k-th element of
##                  ee corresponds the the k-th element of reality_graphs.
#### Output ####
## Average of the squared differences between observed and extimated event numbers.
compare_estimate_reality <- function(ee,reality_graphs) {
  MSEs <- rep(0,length(reality_graphs))
  for(k in 1:length(reality_graphs)) {
    if(length(E(reality_graphs[[k]]))==0) {
      A <- 0
    } else{
      A <- as_adjacency_matrix(reality_graphs[[k]],attr="count",sparse=TRUE)
    }
    MSEs[k] <- sum((ee$ee[[k]]-A)^2)
  }
  return(mean(MSEs))
}

## Generates random events in the network according to the Cox Model with
## parameter theta and covariate vector given in covariates. The simulation ends,
## if reason="time" when the the first simulated event would occur after
## time E, or if reason="events" after occurence of E events. Currently, the
## simulation can only be computed within one covariate segment, i.e., a value
## k has to be specified and then only the k-th segment of the covariates is used.
## Simulation of events on a longer time scale can be achieved by manually calling
## this function on different covariate segments.
##### Input #####
## theta      - Parameter Vector to be used
## n          - Size of the network
## covariates - Covariates which should be used for intensity computation
## k          - Entry of covariates which should be used
## start_time - Start time of the simulation, most likely this equals covariates$time[k]
## E          - Either the number of events which should be simulated or the time point
##              until which the simulation should run, see also "reason" below
## reason     - Either "time" or "events": If "events" then "E" many events will be
##              simulated, if "time" the simulation will run until time "E" is exceeded
##              by an event
## seed       - If different from NULL this value is used as seed value
##              for random number calculation.
##### Output ####
## List of three elements
## graph       - An igraph obejct which contains a graph of size "n" with eddge attribute
##               "weight" which gives the number of simulated events on the corresponding edge.
## events      - An event list of the simulated events in the format ef event lists
## intensities - A vector which contains the intensities of all edges, entry j of this
##               vector corresponds to edge covariates$edgelist[[k]][j,]
simulate_network <- function(theta,n,covariates,k,start_time,E,reason,seed=NULL) {
  set.seed(seed)
  
  ## Compute all intensities
  intensities <- .Call("intensities",as.double(theta),covariates$covars[[k]])
  S <- sum(intensities)
  
  ## Simulate time points of events
  if(reason=="events") {
    ## We know that "E" many events are to be considered
    event_times <- start_time+cumsum(rexp(E,S))
  } else {
    ## Generate events until we reach
    event_times <- start_time+cumsum(rexp(1000,S))
    while(max(event_times)<=E) {
      event_times <- c(event_times,max(event_times)+cumsum(rexp(1000,S)))
    }
    event_times <- event_times[event_times<=E]
  }
  
  ## Simulate locations of events
  helper          <- rmultinom(length(event_times),1,intensities)
  event_locations <- apply(helper==1,2,which)
  
  ## Set initial data including output
  N          <- length(covariates[[1]])
  G          <- graph_from_adjacency_matrix(Matrix(0,nrow=n,ncol=n),mode="undirected")
  events     <- matrix(1,ncol=5,nrow=length(event_times))
  events[,1] <- event_times
  events[,5] <- k
  
  cat("Events have been simulated.\n")
  
  ## Generate event list and graph
  for(i in 1:length(event_times)) {
    ## Extract next event
    event_loc <- event_locations[i]
    event_t   <- event_times[i]
    
    ## Get edge of event
    eoe <- c(covariates$edgelist[[k]][event_loc,])
    
    ## Add event to event list
    events[,c(2,3)] <- eoe
    
    ## Add the event to graph
    if(are_adjacent(G,eoe[1],eoe[2])==FALSE) {
      ## This was the first for the pair => Create the edge
      G <- add_edges(G,eoe,attr=list(weight=1))
    } else {
      ## The edge alread exists => Increase the count
      E(G,eoe)$weight <- E(G,eoe)$weight+1
    }
  }
  
  return(list(graph=G,events=events,intensities=intensities))
}

## Given the events in "events" a graph "G" is returned which contains an edge (i,j)
## for all pairs for which there was at least one event between times t0 and t1.
## The graph "G" has an edge attribute "count" which gives the number of events
#### Input ####
## events - An events matrix, see MLE
## t0,t1  - Start and end point of period of interest, respectively.
## n      - Number of vertices which the resulting graph should have.
#### Output ####
## An igraph object with edge attribut "count" which contains the number of events
## on that edge in the period [t0,t1]. If not all "n" vertices are mentioned in the
## time window in events, the missing number of singletons is just added.
graph_from_events <- function(events,t0,t1,n) {
  ## Create Graph from Event List
  indis <- which(events[,1]>=t0 & events[,1]<=t1)
  G     <- graph_from_edgelist(events[indis,c(2,3)],directed=FALSE)
  E(G)$count <- 1
  G     <- simplify(G,edge.attr.comb = "sum")
  
  ## Add vertices if necessary
  if(length(V(G))<n) {
    G <- add_vertices(G,n-length(V(G)))
  }
  
  return(G)
}

## This function computes the test static together with the scaling and shifting. It is
## necessary to run the functions MLE and parametric_MLE first and provide their outputs
## to this function. In Priciple it does not matter for which time points t0 was called.
## However, the quality of the estimation of the interals is probably better if MLE was
## called on a finer vector t0. The output of MLE is the non-parametric estimator and
## the result of parametric_MLE is the parametric estimator to which it is compared.
##### Input ####
## est        - Result of the function MLE => Nonparametric Estimator
## param_est  - Result of the function parametric_MLE => Parametric Estimator
## covariates - List of covariates (see MLE for details)
## events     - Observed events (see MLE for details)
## h          - Bandwidth to be used (has to be the same as the one used in MLE)
## Kf         - Kernel function (see MLE for details), default is a triangular kernel.
## wf         - Weight function, must be supported on [0,1] and be non-negative, the function
##                is scaled to the interval [0,T] without adjusting for the integral. Default
##                is cutting-off the first and last percent of [0,T].
## Kfour      - Kernel constant, default is the constant for the triangular kernel
## Delta      - The grid-mesh for the grid which is used for approximating the integrals. The
##                choice of the value influences the run time.
#### Output ####
## List of five elements
## Tn        - Value of the test-statistic (includes the factor r_n)
## An        - Centering
## B         - Asympotic variance of centred test statistic
## h         - Bandwidth used
## Tn_scaled - Centred and scaled test statistic, under H0 this is N(0,1) distributed
L2_test_statistic <- function(est,param_est,covariates,events,h,Kf=K,wf=standard_weight,Kfour=151/630,Delta=1) {
  ## Some useful variables
  ## Setup Grid for Integral Approximations
  grid <- seq(from=est$t_0[1],to=ceiling(T/Delta)*Delta,by=Delta)
  
  ## Create weight vector
  w <- sapply(grid/T,wf)
  
  ## Create vector of covariate and estimator indices
  ci <- sapply(grid,covar_index,covar_times=covariates$time)
  ei <- sapply(grid,covar_index,covar_times=est$t_0)
  
  ## Create a vector which contains number of active edges at a time.
  nqpn <- rep(0,length(grid))
  for(i in 1:length(grid)) {
    nqpn[i] <- dim(covariates$edgelist[[ci[i]]])[1]
  }
  
  ## Compute smoothed sparsity \bar{p}_n
  snqpn <- rep(0,length(grid))
  for(i in 1:length(grid)) {
    kweights <- 1/h*Kf((grid-grid[i])/h)
    snqpn[i] <- sum(kweights*nqpn)*Delta
  }
  
  ## Compute test statistic
  Tn <- 0
  for(k in 1:length(est$t_0)) {
    i <- which(ei==k)
    Tn <- Tn+sum((est$estimates[k,]-param_est$estimate)^2)*sum(snqpn[i]*w[i])*Delta
  }
  
  ## Compute Sigma Inverses
  SigmaInv <- list()
  for(k in 1:length(covariates$time)) {
    SigmaInv[[k]] <- solve(Reduce('+',lapply(covariates$covars[[k]],sig_help,theta=param_est$estimate))/dim(covariates$edgelist[[k]])[1])
  }
  
  ## Compute B
  B <- 0
  for(k in 1:length(covariates$time)) {
    i <- which(ci==k)
    B <- B+sum(diag(SigmaInv[[k]]%*%SigmaInv[[k]]))*sum(w[i]^2)*Delta
  }
  B <- 4*Kfour*B
  
  ## Compute A_n
  N <- dim(events)[1]
  An <- 0
  for(k in 1:N) {
    it <- inner_integral_testA(events[k,1],Kf,h,grid,Delta,w,snqpn,SigmaInv,covariates$time)
    j <- covariates$index[[events[k,5]]][events[k,2],events[k,3]]
    X <- matrix(covariates$covars[[events[k,5]]][[j]],ncol=1)
    An <- An+t(X)%*%it%*%X
    if(k %% 1000 ==0) {cat("Step ",k," of ",N,"\n")}
  }
  
  ##Output
  return(list(Tn=Tn,An=An,B=B,h=h,Tn_scaled=(sqrt(h)*Tn-An/sqrt(h))/sqrt(B)))
}

## Example weight function which cuts-off one percent of boundary
standard_weight <- function(x) {
  if(x>=0.01 & x <=0.99) {
    return(1)
  } else {
    return(0)
  }
}


######### Internal Functions ####################################################
## The following functions are needed for the execution of the above functions ##
## but in most scenarios they are probably not directly called by the user.    ##
#################################################################################

## Computes MLE
MLE_raw <- function(t_0,mu_init,h,events,covariates,T,K,KIntegrate,iterlim,gradtol) {
  cat("Begin Maximizing Likelihood\n")
  output <- nlm(nlikelihood,mu_init,hessian=FALSE,check.analyticals = FALSE,iterlim=iterlim,gradtol=gradtol,print.level=2,t_0=t_0,h=h,events=events,covariates=covariates,T=T,K=K,KIntegrate=KIntegrate)
  cat("Done Maximizing Likelihood\n")
  
  output$hessian <- attr(nlikelihood(output$estimate,t_0,h,events,covariates,T,K,KIntegrate),"hessian")
  return(output)
}


## This function computes the negative likelihood, the meaning of the
## input variables is the same as for MLE, however t_0 cannot be a vector
nlikelihood <- function(mu,t_0,h,events,covariates,T,K,KIntegrate) {
  ## Compute Likelihood
  oraw <- likelihood(mu,t_0,h,events,covariates,T,K,KIntegrate)
  
  ## Change the sign of everything
  o <- -oraw
  attr(o,"gradient") <- -attr(oraw,"gradient")
  attr(o,"hessian") <- -attr(oraw,"hessian")
  
  return(o)
}





## Collection of some supportive Functions
int1 <- function(ev,cov,ind,t_0,h,K) {
  return(1/h*K((ev[1]-t_0)/h)*cov[[ind[ev[2],ev[3]]]])
}

expo <- function(x,mu) {return(exp(sum(x*mu)))}
xexpo <- function(x,mu) {return(x*exp(sum(x*mu)))}
x2expo <- function(x,mu) {return(as.matrix(x)%*%t(as.matrix(x))*exp(sum(x*mu)))}

summing_function1 <- function(co,mu) {
  emux <- lapply(co,expo,mu=mu)
  return(Reduce('+',emux))
}

summing_function2 <- function(co,mu) {
  emux <- lapply(co,xexpo,mu=mu)
  return(Reduce('+',emux))
}

summing_function3 <- function(co,mu) {
  emux <- lapply(co,x2expo,mu=mu)
  return(Reduce('+',emux))
}

integrator <- function(t_low,t_high,x,KIntegrate,T,t_0,h) {
  return(KIntegrate((max(t_low,t_0-h,0)-t_0)/h,(min(t_high,t_0+h,T)-t_0)/h)*x)
}

## Computes the maximum likelihood estimator in the locally
## linear model. This is the internal function which is usually
## not called by the user. See LL_MLE
LL_MLE_raw <- function(t_0,mu_init,h,events,covariates,T,iterlim,gradtol) {
  cat("Begin Maximizing Likelihood\n")
  output <- nlm(likelihoodLL,mu_init,hessian=FALSE,check.analyticals = FALSE,iterlim=iterlim,gradtol=gradtol,print.level=2,t_0=t_0,h=h,events=events,covariates=covariates,T=T)
  cat("Done Maximizing Likelihood\n")
  output$hessian <- attr(likelihoodLL(output$estimate,t_0,h,events,covariates,T),"hessian")
  return(output)
}

## Compute the "likelihood" for the locally linear model
## This function is for internal use and usually does
## not have to be called by the user.
likelihoodLL <- function(mu,t_0,h,events,covariates,T) {
  ## Split Parameter
  q <- length(mu)
  mu0 <- mu[1:(q/2)]
  mu1 <- mu[(q/2+1):q]
  
  K=Khalf
  
  cat("Call with mu=",mu,"\n")
  
  ################################################
  ## Compute the first integral and derivatives ##
  ################################################
  relevant_indices <- which(events[,1]>=max(0,t_0-h) & events[,1]<min(T,t_0))
  covar_indices <- events[relevant_indices,5]
  
  A <- mapply(int1,lapply(relevant_indices,function(i) events[i,]),covariates$covars[covar_indices],covariates$index[covar_indices],MoreArgs = list(t_0=t_0,h=h,K=Khalf),SIMPLIFY=FALSE)
  B <- Map('*',A,events[relevant_indices,1]-t_0)
  
  integral1_derivative <- c(Reduce('+',A),Reduce('+',B))
  integral1 <- sum(integral1_derivative[1:(q/2)]*mu0)+sum(integral1_derivative[(q/2+1):q]*mu1)
  integral1_hessian <- 0
  
  #################################################
  ## Compute the second integral and derivatives ##
  #################################################
  covar_timestamps <- sort(unique(covar_indices))
  if(covar_timestamps[1]>1) {covar_timestamps <- c(covar_timestamps[1]-1,covar_timestamps)}
  ibounds <- (c(covariates$time[covar_timestamps],min(T,t_0))-t_0)/h
  
  A1 <- mapply(sf,covariates$covars[covar_timestamps],ibounds[1:(length(ibounds)-1)],ibounds[2:length(ibounds)],MoreArgs=list(mu0=mu0,mu1=mu1,h=h,mode=1),SIMPLIFY = FALSE)
  A2 <- mapply(sf,covariates$covars[covar_timestamps],ibounds[1:(length(ibounds)-1)],ibounds[2:length(ibounds)],MoreArgs=list(mu0=mu0,mu1=mu1,h=h,mode=2),SIMPLIFY = FALSE)
  A3 <- mapply(sf,covariates$covars[covar_timestamps],ibounds[1:(length(ibounds)-1)],ibounds[2:length(ibounds)],MoreArgs=list(mu0=mu0,mu1=mu1,h=h,mode=3),SIMPLIFY = FALSE)
  A4 <- mapply(sf,covariates$covars[covar_timestamps],ibounds[1:(length(ibounds)-1)],ibounds[2:length(ibounds)],MoreArgs=list(mu0=mu0,mu1=mu1,h=h,mode=4),SIMPLIFY = FALSE)
  A5 <- mapply(sf,covariates$covars[covar_timestamps],ibounds[1:(length(ibounds)-1)],ibounds[2:length(ibounds)],MoreArgs=list(mu0=mu0,mu1=mu1,h=h,mode=5),SIMPLIFY = FALSE)
  A6 <- mapply(sf,covariates$covars[covar_timestamps],ibounds[1:(length(ibounds)-1)],ibounds[2:length(ibounds)],MoreArgs=list(mu0=mu0,mu1=mu1,h=h,mode=6),SIMPLIFY = FALSE)
  
  
  integral2 <- Reduce('+',A1)
  integral2_derivative <- c(Reduce('+',A2),Reduce('+',A3))
  integral2_hessian <- matrix(0,ncol=q,nrow=q)
  integral2_hessian[1:(q/2),1:(q/2)] <- Reduce('+',A4)
  integral2_hessian[(q/2+1):q,1:(q/2)] <- Reduce('+',A5)
  integral2_hessian[1:(q/2),(q/2+1):q] <- Reduce('+',A5)
  integral2_hessian[(q/2+1):q,(q/2+1):q] <- Reduce('+',A6)
  
  ## Prepare Output
  output <- -(integral1-integral2)
  attr(output,"gradient") <- -(integral1_derivative-integral2_derivative)
  attr(output,"hessian") <- -(integral1_hessian-integral2_hessian)
  
  #  output <- -(integral2)
  #  attr(output,"gradient") <- -(integral2_derivative)
  #  attr(output,"hessian") <- -(integral1_hessian-integral2_hessian)
  
  return(output)
}

Khalf <- function(x) {
  if(x>=0 | x<(-1)) {
    return(0)
  } else {
    return(2*(1+x))
  }
}

## Supporting functions called by the likelihood routine for locally linear estimation.
## These are usually not called by the user.
sf <- function(covs,lb,ub,mu0,mu1,h,mode) {
  return(Reduce('+',lapply(covs,iLL,lb=lb,ub=ub,mu0=mu0,mu1=mu1,h=h,mode=mode)))
}

iLL <- function(x,lb,ub,mu0,mu1,h,mode) {
  alpha <- sum(mu1*x)*h
  fac <- exp(sum(mu0*x))
  
  if(mode==1) {
    p <- 0
  }
  if(mode==2) {
    p <- 0
    fac <- x*fac
  }
  if(mode==3) {
    p <- 1
    fac <- h*x*fac
  }
  if(mode==4) {
    p <- 0
    fac <- as.matrix(x)%*%t(as.matrix(x))*fac
  }
  if(mode==5) {
    p <- 1
    fac <- h*as.matrix(x)%*%t(as.matrix(x))*fac
  }
  if(mode==6) {
    p <- 2
    fac <- h^2*as.matrix(x)%*%t(as.matrix(x))*fac
  }
  return(fac*(KintEx(ub,alpha,p)-KintEx(lb,alpha,p)))
}

KintEx_os <- function(x,alpha,p) {
  if(p==0) {
    if(alpha==0) {
      return(1/2*x^2+x+1/2)
    } else {
      return((x+1)/alpha*exp(alpha*x)-1/alpha^2*(exp(alpha*x)-exp(-alpha)))
    }
  }
  if(p==1) {
    if(alpha==0) {
      return(1/3*x^3+1/2*x^2-1/6)
    } else{
      return(x*(x+1)/alpha*exp(alpha*x)-(2*x+1)/alpha^2*exp(alpha*x)-1/alpha^2*exp(-alpha)+2/alpha^3*exp(alpha*x)-2/alpha^3*exp(-alpha))
    }
  }
  if(p==2) {
    if(alpha==0) {
      return(1/4*x^4+1/3*x^3+1/12)
    } else {
      return((x+1)*x^2/alpha*exp(alpha*x)-(3*x^2+2*x)/alpha^2*exp(alpha*x)+1/alpha^2*exp(-alpha)+(6*x+2)/alpha^3*exp(alpha*x)+4/alpha^3*exp(-alpha)-6/alpha^4*(exp(alpha*x)-exp(-alpha)))
    }
  }
}

KintEx <- function(x,alpha,p) {
  if(x<=-1) {
    return(0)
  } else if(x<=0) {
    return(2*KintEx_os(x,alpha,p))
  } else {
    return(2*KintEx_os(0,alpha,p))
  }
}

event_estimates_internal <- function(t,parameter,covariates,size) {
  t_1 <- t[1]
  t_2 <- t[2]
  
  pi_lower <- max(which(parameter$t_0<=t_1))
  pi_upper <- max(which(parameter$t_0<t_2))
  
  ci_lower <- max(which(covariates$time<=t_1))
  ci_upper <- max(which(covariates$time<t_2))
  
  if(pi_lower==-Inf | ci_lower==-Inf) {stop("Parameters or covariates do not cover the prediction range.\n"); return(1)}
  
  time_span <- unique(sort(c(parameter$t_0[pi_lower:pi_upper],covariates$time[ci_lower:ci_upper])))
  time_span <- c(t_1,time_span[which(time_span>t_1)],t_2)
  #  time_span <- c(time_span[which(time_span>=t_1)],t_2)
  
  
  information_matrix <- matrix(0,ncol=4,nrow=length(time_span)-1)
  run_parameter <- pi_lower
  run_covariate <- ci_lower
  for(k in 1:(length(time_span)-1)) {
    if(run_parameter<pi_upper) {
      if(parameter$t_0[run_parameter+1]<=time_span[k]) {
        run_parameter <- run_parameter+1
      }
    }
    if(run_covariate<ci_upper) {
      if(covariates$time[run_covariate+1]<=time_span[k]) {
        run_covariate <- run_covariate+1
      }
    }
    
    information_matrix[k,] <- c(time_span[k],time_span[k+1],run_parameter,run_covariate)
  }
  
  ee <- Reduce('+',apply(information_matrix,1,single_estimate,parameter=parameter,covariates=covariates,size=size))
  
  return(ee)
}

single_estimate <- function(info,parameter,covariates,size) {
  output <- sparseMatrix(covariates$edgelist[[info[4]]][,1],covariates$edgelist[[info[4]]][,2],x=0.5,dims=c(size,size))
  N <- dim(covariates$edgelist[[info[4]]])[1]
  
  p <- length(parameter$estimates[1,])/2
  for(k in 1:N) {
    ind <- covariates$edgelist[[info[4]]][k,]
    output[ind[1],ind[2]] <- exp(sum(parameter$estimates[info[3],1:p]*covariates$covars[[info[4]]][[k]]))*(info[2]-info[1])
  }
  
  return(output)
}

## Reteurns covariate index belonging to time point
covar_index <- function(t_0,covar_times) {
  return(length(covar_times[covar_times<=t_0]))
}

## Returns matrix of vector
sig_help <- function(x,theta) {
  xm <- matrix(x,ncol=1)
  return(xm%*%t(xm)*exp(sum(theta*x)))
}

## Computes some inner integral for test statistic
inner_integral_testA <- function(s,K,h,grid,Delta,w,snqpn,SigmaInv,cov_times) {
  all_weights <- 1/h*K((grid-s)/h)^2*w/snqpn
  
  q <- dim(SigmaInv[[1]])[1]
  out <- matrix(0,nrow=q,ncol=q)
  
  indis <- which(all_weights!=0)
  sigma_indis <- sapply(grid[indis],covar_index,covar_times=cov_times)
  sigma_indis_unique <- unique(sigma_indis)
  
  if(length(sigma_indis)>0) {
    for(k in 1:length(sigma_indis_unique)) {
      i <- indis[which(sigma_indis==k)]
      out <- out+sum(all_weights[i])*SigmaInv[[k]]%*%SigmaInv[[k]]*Delta
    }
  }
  
  return(out)
}