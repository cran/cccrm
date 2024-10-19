
ccc_sim_data_aux <- function(n = 30, nrep = 1,b = NULL, g = NULL, mu = 0, sa = 1,
                         sab = 0, sag = 0, bg = NULL, se = 1,  
                         extra.info=TRUE){
  
  options(dplyr.summarise.inform = FALSE)
  
  if ((extra.info %in% c(TRUE,FALSE))==FALSE) extra.info<-FALSE
  
# Checking and setting variances  
  
  # Times
  if(!is.null(g)){
    p <- length(g)
 
  }else{
    p <- 1
    sag <- 0
  }
  
  if(!is.null(b)){
    k <- length(b)
 
  }else{
    k <- 1
    sab <- 0
  }
  
  if ((length(bg) != k*p) & (!is.null(bg))) {
    stop("Warning: Method-time effect vector doesn't match its design length.")
    a <- ab <- ag <- e <- mu_jt <- y <- NULL
  }
  
  if(any(c(sa,sab,sag,se) < 0)){
    stop("Warning: Variance components should be positive")
  }
  
  id<- rep(1:n,k*p)
  met<- rep(1:k,n*p)
  
  if (p>1){
    times <- rep(1:p,n*k)
  } else {times<-NULL}

  
  dades <-  expand.grid(id = 1:n, met = 1:k, times = 1:p)
  
  # Overall mean & Random error
  
  dades <- dades |> dplyr::mutate(mu=mu)
  
  
  # Subject effect
  dades <- dades |> dplyr::group_by(id) |> dplyr::mutate(a= rnorm(1,0,sa))
  
  
  # Method effect
  if(!is.null(b)){
    dades.met <- data.frame(met=1:k,b=b)
    dades<-dades |> dplyr::left_join(dades.met,by="met")
  }else{
    dades <- dades |> dplyr::mutate(b = 0)
  }
  
  # Time effect
  if(!is.null(g)){
    dades.time <- data.frame(times=1:p,g=g)
    dades<-dades |> dplyr::left_join(dades.time,by="times")
  }else{
    dades <- dades |> dplyr::mutate(g = 0)
  }
  
  # Subject-method effect
  if(sab>0){
    dades<-dades |> dplyr::group_by(id,met) |> dplyr::mutate(ab=rnorm(1,0,sab))
  }else{
    dades <- dades |> dplyr::mutate(ab = 0)
  }
  
  # Subject-time effect
  if(sag>0){
    dades<-dades |> dplyr::group_by(id,times) |> dplyr::mutate(ag=rnorm(1,0,sag))
  }else{
    dades <- dades |> dplyr::mutate(ag = 0)
  }
  
  if(!is.null(bg)){
  # Method-time effect

    dades.met_time <- dades |> dplyr::group_by(met,times) |>
      dplyr::summarise(n=n()) |> dplyr::select(-n) |>
      dplyr::ungroup() |> dplyr::mutate(bg = bg)
    
    dades<-dades |> dplyr::left_join(dades.met_time,by=c("met","times"))
    
    # Cálculo var methods-time
    # dades.sbg <-dades.met_time |> dplyr::left_join(dades.met,by="met") |>
    #   dplyr::mutate(mu_jt=(b+bg)^2) |> dplyr::ungroup()  |> dplyr::summarise(sbg=sum(mu_jt)) 
    # 
    # sbg<-dades.sbg$sbg/(p*(k-1))
    dades.sbg <-dades.met_time |> dplyr::left_join(dades.met,by="met") |>
      dplyr::mutate(mu_jt=(b+bg)) |> group_by(times) |> 
      summarise_at(vars(mu_jt),diff) |> ungroup() |>
      summarise(sbg=sum(mu_jt^2))
    
    sbg<-dades.sbg$sbg/(p*k*(k-1))
    
  }else{
    dades <- dades |> dplyr::mutate(bg = 0)
    # Cálculo var methods
    if(!is.null(b)){ 
      dades.b <-dades.met |> 
        dplyr::mutate(mu_jt=b^2) |> dplyr::summarise(sb=sum(mu_jt)) |> dplyr::ungroup()
      sbg<-dades.b$sb/(k-1)}else{
        sbg<-0
      }
  }
  
  # Random error
  dades <- dades|> dplyr::group_by(id,met,times) |>
    dplyr::slice(rep(1:dplyr::n(),nrep)) |>
    dplyr::mutate(e=rnorm(nrep,0,se))
  
  # Outcome
  dades <-dades |> dplyr::mutate(y=mu+a+b+g+ab+ag+bg+e)
  
  
  if ((p>1) & (k>1)) {
    df <- dades |> dplyr::select(y,id,met,times)
  } 
  if ((p==1) & (k>1)){
    df <- dades |> dplyr::ungroup() |> dplyr::select(y,id,met)
  }
  if ((p>1) & (k==1)){
    df <- dades |> dplyr::ungroup() |> dplyr::select(y,id,times)
  }
  if ((p==1) & (k==1)){
    df <- dades |> dplyr::ungroup() |> dplyr::select(y,id)
  }
  
  # Compute ICC
  if (extra.info){
    pc <- (sa^2+sag^2)/(sa^2+sag^2+sab^2+sbg+se^2)
    
    cat("ICC/CCC simulated:",pc,"\n")
    varcomp=data.frame(var_a=sa^2,var_ag=sag^2,
                       var_ab=sab^2,
                       var_bg=sbg,
                       var_se=se^2)
    
    cat("Variance Components simulated:","\n")
    print(varcomp)
  }
  

  return(df)
}

#' Data simulation using fixed and random effects
#'
#' @export
#' @description
#' The fixed effects and standard deviations of random effects can be set to specific values or, alternatively, obtained from an object of class \code{lme}. 
#' @param n Integer. Number of subjects
#' @param nrep Integer. Number of replicates
#' @param nsim Integer. Number of data sets simulated.
#' @param model Object of class \code{lme}.
#' @param b Vector. Method fixed effects.
#' @param g Vector. Time fixed effects.
#' @param mu Integer. Overall mean.
#' @param sa Integer. Standard deviation of subject's random effect. 
#' @param sab Integer. Standard deviation of subject-method interaction's random effect.
#' @param sag Integer. Standard deviation of subject-time interaction's random effect.
#' @param bg Vector. Method-time interaction's fixed effects. The vector of effects have to be ordered by method and time.
#' @param se Integer. Standard deviation of random error effect.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in \code{\link[furrr]{furrr_options}}.
#' @param workers Integer. Number of cores to be used for parallellization. Default is 15. Capped to number of available cores minus 1.
#' @param extra.info Logical. Should the information about CCC/ICC and variance components simulated be shown?  Default is set to TRUE.
#' @param ... To pass further arguments.
#' @return A data frame with the simulated data. 
#' @details
#' Random effects are simulated as normal distributions with mean 0 and the correspondign standard deviations. The simulated data is obtained
#' as the addition of the simulated values and the fixed efffects. Parallel computation is used except if data is simulated from an object of class `lme`. In this case.
#' data is simulated using the \code{\link[nlmeU]{simulateY}} function from \code{nlmeU} package.
#' @seealso \code{\link{ccc_vc}}
#' @examples
#'\donttest{
#' # # Reliability data: 
#' # 50 subjects, one method, one time, 2 replicates
#' # Overall mean: -0.25; Subjects standard deviation: 1.5, Random error standard deviation: 1
#' set.seed(101)
#' df <- ccc_sim_data(n=50, b = NULL, g = NULL, mu = -0.25, sa = 1.5, se = 1, nrep=2)
#' 
#' # Method comparison data (non-longitudinal)
#' # 50 subjects, two methods, 2 replicates
#' # Overall mean: -0.25; Subjects standard deviation: 1.5, Random error standard deviation: 1
#' # Difference of means between methods 2 and 1: 1
#' # Three data sets simulated
#' 
#' set.seed(202)
#' df <- ccc_sim_data(n=50, nsim=3,b = c(0,1), mu = -0.25, sa = 1.5, se = 1, nrep=2)
#' 
#' # Method comparison data (longitudinal)
#' # 50 subjects, two methods, 3 times, 1 replicate, 
#' # Overall mean: -0.25; Subjects standard deviation: 1.5, Random error standard deviation: 1
#' # Difference of means between methods 2 and 1: 1
#' # Difference of means between times 3,2 and 1 respectively: 0.5 and 0.25.
#' # Subject-methods interaction standard deviation: 0.25
#' # Subject-times interaction standard deviation: 0.5
#' # Same difference of means at each time
#' 
#' set.seed(202)
#' df <- ccc_sim_data(n=50, b = c(0,1), g=c(0,0.25,0.5), mu = -0.25, sa = 1.5, 
#' sab=0.25,sag=0.5,se = 1, nrep=2)
#' 
#' # Simulate data using the estimates of a linear mixed model
#' set.seed(2024)
#' df3 <- ccc_sim_data(n=50, b = c(0,1), g=c(0,0.25,0.5), mu = -0.25, sa = 1.5, 
#'                     sab=0.25,sag=0.5,bg=c(0,0.5,0.75,0,1,1),se = 1, nrep=2)
#' mod3 <- lme_model(df3,"y","id","times","met",control.lme=nlme::lmeControl(opt = 'optim'))
#' ccc_sim_data(nsim=10,model=mod3)
#' }
#' 
ccc_sim_data <- function(n = 30, nrep = 1, nsim=1, model=NULL, b = NULL, g = NULL, mu = 0, sa = 1,
                         sab = 0, sag = 0, bg = NULL, se = 1,  
                         future_seed=TRUE,workers=15,extra.info=TRUE,...
                         ){
  
  options(dplyr.summarise.inform = FALSE)

if (!is.null(model)){
  if(inherits(model,"lme")){
    
    if (!is.null(model$data$time)){
      
      orig.data<-model$data |> select(ind,met,time)
    
      }else{
      
        if (!is.null(model$data$met)){
          orig.data<-model$data |> select(ind,met)
      
        }else{
        
          if (!is.null(model$data$ind)){
            orig.data<-model$data |> select(ind)
        
            }else{
          stop("Model has to be generated with cccrm package")
          }
        }
      }
  }
    
    df <- nlmeU::simulateY(model, nsim = nsim, verbose = F)
    #row.names(sim_data) <- 1:model$dims$N
    df<-data.frame(df)
    new.data<-cbind(df,orig.data)
    
    
    sim_data<-new.data |>
      tidyr::pivot_longer(
        cols = starts_with("X"),
        names_to = "sample",
        names_prefix = "X",
        values_to = "y"
      ) |> arrange(sample)
    
    
    
    
  
  
}else{
  
  if (sum(g)!=0){
    message("Time effects are transformed to sum 0")
    g<-g-mean(g)
    message("New g values:","\n",paste(as.character(g),collapse=",")
    )
  }
  if (sum(b)!=0){
    message("Method effects are transformed to sum 0")
    b<-b-mean(b)
    message("New b values:","\n",paste(as.character(b),collapse=",")
    )
  }
  if (sum(bg)!=0){
    message("Method-Time effects are transformed to sum 0")
    bg<-bg-mean(bg)
    message("New bg values:","\n",paste(as.character(bg),collapse=",")
    )
  }
  
  # Parallellization
  ncores <- parallelly::availableCores(omit = 1)
  if(workers >= ncores){
    workers <- ncores
  }
    with_progress({
      p <- progressr::progressor(steps = nsim)
      sim_data <- furrr::future_map_dfr(1:nsim,~{
        p()
        Sys.sleep(.2)
        if (.x>1) extra.info<-FALSE
        df<-ccc_sim_data_aux(n = n, b = b, g = g, mu = mu, sa = sa,
                             sab = sab, sag = sag, bg = bg, se = se,  
                             nrep = nrep,extra.info=extra.info)
        df<-df |> mutate(sample=.x)
        
        
        return(df)
      }, .options = furrr::furrr_options(seed = future_seed),
      p = p)
    })
}

return(sim_data)

}
