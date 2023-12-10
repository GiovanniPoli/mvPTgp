### Library & functions  -----
library(mvtnorm)
library(cowplot)
library(telegram.bot)
library(reshape2)
library(readr)
library(survival)
library(cowplot)
library(reshape2)
library(BayesLogit)
library(metafor)
library(mixmeta)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)

erfinv = function(x) qnorm((1+x)/2)/sqrt(2)
expit  = function(z) 1/(1+exp(-z))
logit  = function(p) if(all(p > 1e-16)){log(p)-log(1-p)}else{ log(p + 1e-16) - log( 1 - p - 1e-16) }

g0     = function(x, scale = 3.5){exp(log(2 * scale) - log(pi * { x * x + scale * scale }))}
G0     = function(q, scale = 3.5){ {2/pi} * atan(q/scale) }
G0inv  = function(p, scale = 3.5){ scale * tan({ pi * p }/2) }

GP_parameters.gen     = function(B, Q, correlated_levels, star_model_deep, correlation_star_matrix, c = 1,
                                 G0 = NULL, G0inv = NULL){
  Us    = G0(Q)  
  n     = nrow(Us)
  
  Pi_star = matrix( 0:(2^star_model_deep-1)/2^star_model_deep,
                    ncol = 2^star_model_deep, nrow = n, byrow = T)
  
  Pi_u    = vector("list", correlated_levels)
  
  Pi_u[[1]] = rbind(0, Us[,2], 1)
  Pi_u[[2]] = rbind(0, Us[,1], Us[,2], Us[,3], 1)
  Pi_u[3 :correlated_levels] = lapply(3:correlated_levels, function(j) apply(Us, 1,  function(u){
    Neps = (0:(2^(j-2)-1)) / 2^(j-2)
    c(Neps*u[1], u[1] + Neps*(u[2]-u[1]),
      u[2] + Neps*(u[3]-u[2]),
      u[3] + Neps*(1-u[3]),1)
  }))
  
  Index_list = vector("list",correlated_levels)
  
  for(deep in 1:correlated_levels) {
    list_to_append =  vector("list",2^(deep-1))
    for(comb in 1:2^(deep-1)){
      lows  = Pi_u[[deep]][comb*2 -1,]
      highs = Pi_u[[deep]][comb*2,]
      list_to_append[[comb]] = apply(Pi_star,2, function(col) col >= lows & col < highs  )
    }
    Index_list[[deep]] = list_to_append
  }
  
  Return_list = lapply(1:correlated_levels, function(deep) lapply(1:2^(deep-1), function(r.) NULL ))
  
  max.iter =  B
  width = options()$width  
  
  if(width>=114){
    barwidth = 100
  }else{
    barwidth = width-14
  }
  
  step  = round(seq(from = 1, to= max.iter, length.out = barwidth))
  perc  = round(step/max.iter*100)
  len   = as.numeric(perc>=10)
  n.cat = 1
  
  t0 = Sys.time()
  for(iter in 1:max.iter){
    Ye = matrix(1, nrow = n, ncol = 2^star_model_deep)
    for(deep in 1:star_model_deep){
      if(deep <= correlated_levels){ 
        Ye0 = expit( sqrt(2*trigamma(deep^2*c))*
                       rmvnorm(n = 2^(deep-1), mean = rep(0,n), sigma = correlation_star_matrix))
        Ye  = Ye * matrix( apply(Ye0,1, function(col) 
          t(sapply(col, function(p) rep(c(p,1-p),each = 2^(star_model_deep-deep))))),
          nrow = n, ncol = 2^star_model_deep)
      }else{
        Ye0 = matrix(rbeta(n = n * 2^(deep-1), shape1 = deep^2*c, 
                           shape2 = deep^2*c),
                     ncol = n, nrow = 2^(deep-1))
        Ye  = Ye * matrix( apply(Ye0,1, function(col) 
          t(sapply(col, function(p) rep(c(p,1-p),each = 2^(star_model_deep-deep))))),
          nrow = n, ncol = 2^star_model_deep)        
      }
    }
    
    ## Impute 1,...,B G_1,..,G_n
    divide = matrix(1,ncol = 1, nrow = n)
    for(deep in 1:correlated_levels){
      next_divide = NULL
      for(comb in 1:2^(deep-1)){
        probs       = rowSums(Index_list[[deep]][[comb]]*Ye)/divide[,comb]
        next_divide = cbind(next_divide, probs, 1-probs)
        Return_list[[deep]][[comb]] = rbind(Return_list[[deep]][[comb]], probs)
      }
      divide = next_divide * Reduce(cbind,apply( divide, 2, function(p) cbind(p,p), simplify = FALSE ))
    }
    
    
    if(step[n.cat]==iter){
      t1 = Sys.time()
      iter.time = difftime(t1,t0,units='secs')/iter
      time.to.end = as.numeric(iter.time)*(max.iter-iter)
      if(time.to.end>=3600*24){
        time.print = "+ 1d"
      }else if(time.to.end>3600){
        time.print =  paste0(round((time.to.end)/3600,1)," h")
      }else if(time.to.end>60){
        time.print =  paste0(round((time.to.end)/60,1)," m")
      }else{        
        time.print =  paste0(round((time.to.end),1)," s")
      }
      bar_print = paste0("|",paste0(rep('=', n.cat), collapse = ''), paste0(rep(' ', barwidth-n.cat), collapse = ''),
                         "|",paste0(rep(' ', 2-len[n.cat]), collapse = ''),perc[n.cat],"%"," (",time.print,")")
      bar_print = paste0(bar_print, paste0(rep(' ', abs(width - nchar(bar_print) +1) ), collapse = ''),collapse = '')
      
      cat('\r',bar_print)
      n.cat= n.cat + 1
    }
    
  }
  
  Ret = lapply(1:correlated_levels, function(deep) lapply(1:2^(deep-1), function(comb){
    data =    logit(Return_list[[deep]][[comb]])
    Ex_Y = colMeans(Return_list[[deep]][[comb]])
    
    Ex_Z   = colMeans(data)
    Var_Z  = cov(data)
    Cov_Z  = cor(data)
    
    cat("\r\tAlmost done! Inverting covariances matrix...\t",paste0("[", 2^(deep-1)-1+comb,"/",2^correlated_levels,"]"))
    list("Ex" = Ex_Z, "Var" = Var_Z, "Cov" = Cov_Z, "Y_Ex" = Ex_Y)
  }))
  Ret 
}
GP_parameters.gen_par = function(B, Q, correlated_levels, star_model_deep, correlation_star_matrix, c = 1,
                                 G0 = NULL, G0inv = NULL, n.cores = 20){
  Us    = G0(Q)  
  n     = nrow(Us)
  
  expit  = function(z) 1/(1+exp(-z))
  logit  = function(p) if(all(p>0)){log(p)-log(1-p)}else{ log(p + 1e-16) - log( 1 - p - 1e-16) }
  
  
  Pi_star = matrix( 0:(2^star_model_deep-1)/2^star_model_deep,
                    ncol = 2^star_model_deep, nrow = n, byrow = T)
  
  Pi_u    = vector("list", correlated_levels)
  
  Pi_u[[1]] = rbind(0, Us[,2], 1)
  Pi_u[[2]] = rbind(0, Us[,1], Us[,2], Us[,3], 1)
  Pi_u[3 :correlated_levels] = lapply(3:correlated_levels, function(j) apply(Us, 1,  function(u){
    Neps = (0:(2^(j-2)-1)) / 2^(j-2)
    c(Neps*u[1], u[1] + Neps*(u[2]-u[1]),
      u[2] + Neps*(u[3]-u[2]),
      u[3] + Neps*(1-u[3]),1)
  }))
  
  Index_list = vector("list",correlated_levels)
  
  for(deep in 1:correlated_levels) {
    list_to_append =  vector("list",2^(deep-1))
    for(comb in 1:2^(deep-1)){
      lows  = Pi_u[[deep]][comb*2 -1,]
      highs = Pi_u[[deep]][comb*2,]
      list_to_append[[comb]] = apply(Pi_star,2, function(col) col >= lows & col < highs  )
    }
    Index_list[[deep]] = list_to_append
  }
  
  n.cores = min(n.cores,detectCores(logical = TRUE)-2)
  cluster = makeCluster(n.cores)
  registerDoSNOW(cluster)
  
  cat("Detected cores:",detectCores(logical = TRUE),
      "\nUsing:",min(n.cores,detectCores(logical = TRUE)-2),
      "\nCheck:", getDoParWorkers(),
      "\nSimulations:", B)
  
  pb       = txtProgressBar(max = B, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts     = list(progress = progress)
  
  store = foreach(it = 1:B, .options.snow = opts,
                  .packages = c("mvtnorm")) %dopar%{
                    
                    Return_list = lapply(1:correlated_levels, function(deep) lapply(1:2^(deep-1), function(r.) NULL ))
                    
                    Ye = matrix(1, nrow = n, ncol = 2^star_model_deep)
                    for(deep in 1:star_model_deep){
                      if(deep <= correlated_levels){ 
                        Ye0 = expit( sqrt(2*trigamma(deep^2*c))*
                                       rmvnorm(n = 2^(deep-1), mean = rep(0,n), sigma = correlation_star_matrix,
                                               method = "svd"))
                        Ye  = Ye * matrix( apply(Ye0,1, function(col) 
                          t(sapply(col, function(p) rep(c(p,1-p),each = 2^(star_model_deep-deep))))),
                          nrow = n, ncol = 2^star_model_deep)
                      }else{
                        Ye0 = matrix(rbeta(n = n * 2^(deep-1), shape1 = deep^2*c, 
                                           shape2 = deep^2*c),
                                     ncol = n, nrow = 2^(deep-1))
                        Ye  = Ye * matrix( apply(Ye0,1, function(col) 
                          t(sapply(col, function(p) rep(c(p,1-p),each = 2^(star_model_deep-deep))))),
                          nrow = n, ncol = 2^star_model_deep)        
                      }
                    }
                    
                    ## Impute 1,...,B G_1,..,G_n
                    divide = matrix(1,ncol = 1, nrow = n)
                    for(deep in 1:correlated_levels){
                      next_divide = NULL
                      for(comb in 1:2^(deep-1)){
                        probs       = rowSums(Index_list[[deep]][[comb]]*Ye)/divide[,comb]
                        next_divide = cbind(next_divide, probs, 1-probs)
                        Return_list[[deep]][[comb]] = rbind(Return_list[[deep]][[comb]], probs)
                      }
                      divide = next_divide * Reduce(cbind,apply( divide, 2, function(p) cbind(p,p), simplify = FALSE ))
                    }
                    Return_list
                    
                  }
  
  stopCluster(cluster) 
  cat("\n")
  Ret = lapply(1:correlated_levels, function(deep) lapply(1:2^(deep-1), function(comb){
    
    Ye   = t(sapply(store, function(obj) obj[[deep]][[comb]]))
    data =    logit(Ye)
    Ex_Y = colMeans(Ye)
    
    Ex_Z   = colMeans(data)
    Var_Z  = cov(data)
    Cov_Z  = cor(data)
    
    cat("\r\tAlmost done! Inverting covariances matrix...\t",paste0("[", 2^(deep-1)+comb,"/",2^correlated_levels,"]"))
    list("Ex" = Ex_Z, "Var" = Var_Z, "Cov" = Cov_Z, "Y_Ex" = Ex_Y)
  }))
  Ret 
}

Pi.gen               = function(Q, truncation = 4, G0 = NULL, G0inv = NULL){
  
  Q = matrix(Q, ncol = 3)
  
  if(is.null(G0) | is.null(G0inv)){
    G0inv = function(q) punif(q, min = 0, max = 1)
    G0    = function(p) punif(p, min = 0, max = 1)
  }
  U = G0(Q)
  
  Pi = lapply(1:truncation, function(d){
    if(d == 1){
      layer = cbind(0, U[,2], 1)
    }else if(d == 2){
      layer = cbind(0, U[,1], U[,2], U[,3], 1)
    }else if(d>=3){
      B00 = U[,1]
      B01 = U[,2] - B00 
      B10 = U[,3] - B00 - B01
      B11 = rep(1, nrow(U)) - (B00 + B01 + B10 )
      
      layer = cbind( cbind(0,sapply( 1:2^(d-2), function(i) B00                   * i/2^(d-2))),
                     sapply( 1:2^(d-2), function(i) B00 +             B01 * i/2^(d-2)),
                     sapply( 1:2^(d-2), function(i) B00 + B01 +       B10 * i/2^(d-2)),
                     sapply( 1:2^(d-2), function(i) B00 + B01 + B10 + B11 * i/2^(d-2)))
    } 
    G0inv(layer)})
  Pi
}

### Simulation data      -----

data_s = data.frame( study  = paste0("study ",rep(1:31, each = 2)),
                     cohort = paste0("cohort ",1:62),
                     marker = rep(c("positive","negative"), 31),
                     tumor  = paste0("tumor ",c(rep(1,20),rep(2,20),rep(3,10),
                                              c(1,1,1,1,2,2,2,2,3,3,3,3))),
                     agent  = paste0("agent ", c(rep(1,10), rep(2,10),
                                               rep(1,10), rep(2,10),
                                               rep(1,6),  rep(2,4),
                                               c(1,1,2,2,1,1,2,2,1,1,2,2))))


median_positive_exp = 3
median_negative_exp = 2.5

median_positive_hn = 3.5
median_negative_hn = 3

median_positive_mix = 4
median_negative_mix = 3.5

keep = c(ls(),"keep","seed")

for(seed in 1:50){
  rm(list = ls()[!(ls() %in% keep)])
  set.seed(seed)

  data = data_s
  
  Names = paste0("sim_8_09_2023_seed_",seed,".rds")
  
  Qs  = NULL
  N_data = NULL
  
  true_medians = NULL
  study_eff    = NULL
  
  studies_sample_size = 20
  delta = rep(1,studies_sample_size)
  
  for(cohort in  1:25){
    
    reff = runif(1,.8,1.2)
    
    pc = data[data$study == paste0("study ",cohort) & data$marker == "positive",]
    nc = data[data$study == paste0("study ",cohort) & data$marker == "negative",]
    
    
    if(nc$agent == pc$agent){
      ag.eff = ifelse(nc$agent == "agent 2", 1, 0) 
    }else{
      cat("Error in data")
    }
    
    if(pc$tumor == "tumor 1"){
      
      lambda_positive = log(2)/ (reff *(median_positive_exp + ag.eff))
      lambda_negative = log(2)/ (reff *(median_negative_exp + ag.eff))
      
      t_pm = sort( rexp(studies_sample_size, rate = lambda_positive))
      t_nm = sort( rexp(studies_sample_size, rate = lambda_negative))
      
      true_medians = c(true_medians,
                       median_positive_exp + ag.eff,
                       median_negative_exp + ag.eff)
      study_eff = c(study_eff,reff,reff)
      
    }else if(pc$tumor == "tumor 2"){
      
      sigma_positive  = reff * (median_positive_hn + ag.eff)/erfinv(.5)/sqrt(2)
      sigma_negative  = reff * (median_negative_hn + ag.eff)/erfinv(.5)/sqrt(2)
      
      t_pm = sort(abs(rnorm(studies_sample_size, sd = sigma_positive)))
      t_nm = sort(abs(rnorm(studies_sample_size, sd = sigma_negative)))
      
      true_medians = c(true_medians,
                       median_positive_hn + ag.eff,
                       median_negative_hn + ag.eff)
      study_eff = c(study_eff,reff,reff)
      
    }else if(pc$tumor == "tumor 3"){
      mix_comp = rbinom(20,1,.5)
      
      lambda_mix_positive = log(2)/ (reff *(median_positive_mix + ag.eff))
      lambda_mix_negative = log(2)/ (reff *(median_negative_mix + ag.eff))
      
      sigma_mix_positive  = reff * (median_positive_mix + ag.eff)/erfinv(.5)/sqrt(2)
      sigma_mix_negative  = reff * (median_negative_mix + ag.eff)/erfinv(.5)/sqrt(2)
      
      t_pm = sort(    mix_comp  *      rexp(studies_sample_size, rate = lambda_mix_positive) +
                        (1-mix_comp) * abs(rnorm(studies_sample_size,   sd =  sigma_mix_positive)))
      
      t_nm = sort(    mix_comp  *      rexp(studies_sample_size, rate = lambda_mix_negative) +
                        (1-mix_comp) * abs(rnorm(studies_sample_size,   sd =  sigma_mix_negative)))
      
      true_medians = c(true_medians,
                       median_positive_mix + ag.eff,
                       median_negative_mix + ag.eff)
      
      study_eff = c(study_eff,reff,reff)
      
    }else{
      cat("Error in data")
    }
    
    km_fit   = survfit(Surv(t, d) ~ 1, data= data.frame(t = 1:studies_sample_size, d = delta))
    
    Qs = rbind( Qs,
                c(t_pm[which.max(km_fit$lower<.5)],
                  t_pm[which.max(km_fit$surv <.5)],
                  t_pm[which.max(km_fit$upper<.5)]),
                c( t_nm[which.max(km_fit$lower<.5)],
                   t_nm[which.max(km_fit$surv <.5)],
                   t_nm[which.max(km_fit$upper<.5)] ))
    
    N_data = rbind(N_data,   c(sum(delta[t_pm < t_pm[which.max(km_fit$lower<.5)]]),
                               sum(delta[t_pm >= t_pm[which.max(km_fit$lower<.5)] &
                                           t_pm <  t_pm[which.max(km_fit$surv<.5)]]),
                               sum(delta[t_pm >= t_pm[which.max(km_fit$surv<.5)] &
                                           t_pm <  t_pm[which.max(km_fit$upper<.5)]]),
                               sum(delta[t_pm >= t_pm[which.max(km_fit$upper<.5)]])),
                   c(sum(delta[t_nm <  t_nm[which.max(km_fit$lower<.5)]]),
                     sum(delta[t_nm >= t_nm[which.max(km_fit$lower<.5)] &
                                 t_nm <  t_nm[which.max(km_fit$surv<.5)]]),
                     sum(delta[t_nm >= t_nm[which.max(km_fit$surv<.5)] &
                                 t_nm <  t_nm[which.max(km_fit$upper<.5)]]),
                     sum(delta[t_nm >= t_nm[which.max(km_fit$upper<.5)]])))
    
  }
  
  
  rownames(data) = c(paste0("study.",rep(1:31,each = 2), rep(c(".+",".-"), 31)))
  
  K_star = matrix(0, ncol = nrow(data), nrow = nrow(data))
  
  for(r in 1:nrow(data)){
    for(c in (r):nrow(data)){
      K_star[r,c] =  K_star[c,r] =   ifelse(data$marker[r] == data$marker[c], 1, 0) +
        ifelse(data$tumor[r]  == data$tumor[c],  2, 0)  +
        ifelse(data$agent[r]  == data$agent[c],  2, 0)
      K_star[r,c] =  K_star[c,r] =   ifelse(data$marker[r] == data$marker[c], K_star[c,r],0 )
      K_star[r,c] =  K_star[c,r] =   ifelse(data$study[r]  == data$study[c], K_star[c,r] + 1, K_star[c,r])
    }
  }

  K_star = K_star / 6
  
  Qs = rbind(Qs, matrix( apply(Qs,2, function(x) median(x)), ncol = 3, nrow = 12,byrow = TRUE ))
  
  
  data[,c(6:8)] = Qs
  colnames(Qs) = colnames(data)[6:8] = c("lower","median","upper")
  colnames(N_data) = c("N.00","N.01","N.10","N.11")
  
  PAR_1 = GP_parameters.gen_par(B = 10000,
                                Q = Qs, 
                                correlated_levels = 6,
                                star_model_deep   = 12,
                                correlation_star_matrix = K_star,
                                c = 5,
                                G0     = G0,
                                G0inv  = G0inv,
                                n.cores = 25)
  
  Pi = Pi.gen(Q = G0(as.matrix(data[,6:8])), truncation = 14)

  G_list = list()
  for(study in 1:50){
    Gi = rep(1/2^8, 2^14)
    for(deep in 1:6){
      Zi = NULL
      for(comb in 1:2^(deep-1)){
        p = PAR_1[[deep]][[comb]]$Y_Ex[study]
        Zi = c(Zi, rep( c(p, 1-p), each = 2^(14-deep) ))
      }
      Gi = Gi * Zi
    }
    
    G_list[[study]] = Gi
  }
  
  G_Ex_list = list()
  for(study in 1:50){
    Gi. = rep(1/2^8, 2^14)
    for(deep in 3:6){
      Zi = NULL
      for(comb in 1:2^(deep-1)){
        p = PAR_1[[deep]][[comb]]$Y_Ex[study]
        Zi = c(Zi, rep( c(p, 1-p), each = 2^(14-deep) ))
      }
      Gi. = Gi. * Zi
    }
    G_Ex_list[[study]] = Gi.
  }
  
  PAR = list()
  for(deep in 1:6){
    l = list()
    for(c in 1:2^(deep-1)){
      it = PAR_1[[deep]][[c]]
      inv = solve(it$Var[1:50,1:50])
      solve(inv)
      l[[c]] =list("Ex"   = it$Ex[1:50], 
                   "Var"  = it$Var[1:50,1:50],
                   "Prec" = inv)
    }
    PAR[[deep]] = l
  }
  
  PRIORS = list("G_star"  = PAR_1,
                "GP"      = PAR,
                "D"       = 6,
                "Pi"      = Pi,
                "G0"      = G0,
                "G0inv"   = G0inv,
                "Ex_list" = G_Ex_list)
  
  ## Starting values & stuffs 
  n = 50
  
  Delta_start = vector(mode = "list", length = 20)
  Delta_start = lapply(rep(20,20), function(x) rep(1,20))
  Ye0_start   = lapply(1:2, function(deep) matrix(0.5, nrow = n, ncol = 2^(deep-1) ))
  
  # N_data = rbind(N_data, matrix(0, ncol  = 5, nrow = 4))
  
  N_data = N_data[1:n,]
  Neps_start  = vector(mode = "list", length = PRIORS$D)
  Neps_start[1:2]  = list(cbind(N_data[,1] + N_data[,2], N_data[,3] + N_data[,4]),
                          cbind(N_data[,1],  N_data[,2], N_data[,3],  N_data[,4]))
  
  for(deep in 3:PRIORS$D){
    Neps_start[[deep]] = Reduce(cbind, 
                                apply(Neps_start[[deep-1]], 2, function(x) 
                                  cbind( ceiling(x/2), floor(x/2)), simplify = FALSE))
  }
  
  
  STATE = list( "Ye0"     = Ye0_start,
                "medians" = data[,6])
  
  
  
  thinning = 25
  size     = 2500
  burn.in  = 1000
  max.iter = burn.in + size * thinning  
  
  width = options()$width  
  
  if(width>=114){
    barwidth = 100
  }else{
    barwidth = width-14
  }
  
  step  = round(seq(from = 1, to= max.iter, length.out = barwidth))
  perc  = round(step/max.iter*100)
  len   = as.numeric(perc>=10)
  n.cat = 1
  
  count = 0
  
  MCMC_list = list()
  
  cat("\n")
  
  t0    = Sys.time()
  for(iter in 1:max.iter){
    for(deep in 1:2){
      for(combination in 1:2^(deep-1)){
        
        Ne         =  rowSums(Neps_start[[deep]][,c(2*combination-1,2*combination)]) 
        empty_set  =  Ne == 0  
        
        Ze0_new = vector(length = n, "numeric")
        
        if(sum(empty_set)==0){
          
          try   =  Ne 
          sucs  =  c(Neps_start[[deep]][,c(2*combination-1)])
          
          k     = sucs - try / 2
          w     = BayesLogit::rpg(n, try, logit(STATE$Ye0[[deep]][,combination][!empty_set]) )
          
          
          Omega     = diag(w)
          Sigma_b   = solve(Omega      +   PRIORS$GP[[deep]][[combination]]$Prec)
          Mean_b    = c( Sigma_b%*%( k + c(PRIORS$GP[[deep]][[combination]]$Prec%*%PRIORS$GP[[deep]][[combination]]$Ex)) )
          
          Zi_inc = c(rmvnorm(n = 1, mean = Mean_b, sigma = Sigma_b))
          expit(Zi_inc)
          Ze0_new = Zi_inc
          
        }else if(sum(empty_set)==1){
          
          m2     = PRIORS$GP[[deep]][[combination]]$Ex[empty_set]
          V22    = PRIORS$GP[[deep]][[combination]]$Var[empty_set, empty_set]
          
          Zi_empty = rnorm(1, mean = m2, sd = sqrt(V22))
          
          try   =   Ne[!empty_set]
          sucs  =  c(Neps_start[[deep]][,c(2*combination-1)])[!empty_set]
          
          k     = sucs - try / 2
          w     = BayesLogit::rpg(n-1, try, logit( STATE$Ye0[[deep]][,combination][!empty_set] ))
          
          Omega     = diag(w)
          Sigma_b   = solve(Omega + PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,!empty_set])
          Mean_b    = c(Sigma_b%*% c( k +
                                        PRIORS$GP[[deep]][[combination]]$Ex[!empty_set] %*%  PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,!empty_set] -
                                        PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,empty_set] * (Zi_empty - m2)))       
          
          Zi_inc = c(rmvnorm(n = 1, mean = Mean_b, sigma = Sigma_b))
          
          Ze0_new[!empty_set] = Zi_inc
          Ze0_new[empty_set]  = Zi_empty
          
        }else if( sum(empty_set)== (n-1) ){
          
          m2     = PRIORS$GP[[deep]][[combination]]$Ex[empty_set]
          V22    = PRIORS$GP[[deep]][[combination]]$Var[empty_set, empty_set]
          
          Zi_empty = c(mvn_sampler_sigma(1, mean = m2, sigma = V22))
          
          try   = Ne[!empty_set]
          sucs  =  c(Neps_start[[deep]][,c(2*combination-1)])[!empty_set]
          
          k     = sucs - try / 2
          wi    = BayesLogit::rpg(sum(!empty_set), try, logit( STATE$Ye0[[deep]][,combination][!empty_set] ))
          
          Sigma_b   = 1 / (wi + PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,!empty_set])
          Mean_b    = Sigma_b *  ( k + c(PRIORS$GP[[deep]][[combination]]$Ex[!empty_set] *  
                                           PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,!empty_set]) -
                                     c(PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,empty_set]%*% (Zi_empty - m2)))     
          
          Zi_inc = rnorm(n = 1, mean = Mean_b, sd = sqrt(Sigma_b))
          
          Ze0_new[!empty_set] = Zi_inc
          Ze0_new[empty_set]  = Zi_empty
          
        }else if( sum(empty_set) == n ){
          
          Zi_empty = c(rmvnorm(n = 1, mean = PRIORS$GP[[deep]][[combination]]$Ex , sigma = PRIORS$GP[[deep]][[combination]]$Var))
          Ze0_new = Zi_empty
          
        }else{
          
          m2     = PRIORS$GP[[deep]][[combination]]$Ex[empty_set]
          V22    = PRIORS$GP[[deep]][[combination]]$Var[empty_set, empty_set]
          
          Zi_empty = c( rmvnorm(1, mean = m2, sigma = V22))
          
          try   = Ne[!empty_set]
          sucs  =  c(Neps_start[[deep]][,c(2*combination-1)])[!empty_set]
          
          k     = sucs - try / 2
          w     = BayesLogit::rpg(sum(!empty_set), try, logit( STATE$Ye0[[deep]][,combination][!empty_set] ))
          
          Omega     = diag(w)
          Sigma_b   = solve(Omega + PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,!empty_set])
          Mean_b    = c(Sigma_b%*% ( k +
                                       c(PRIORS$GP[[deep]][[combination]]$Ex[!empty_set] %*%  PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,!empty_set]) -
                                       c(PRIORS$GP[[deep]][[combination]]$Prec[!empty_set,empty_set]%*% (Zi_empty - m2))))     
          
          Zi_inc = c(rmvnorm(n = 1, mean = Mean_b, sigma = Sigma_b))
          
          Ze0_new[!empty_set] = Zi_inc
          Ze0_new[empty_set]  = Zi_empty
        }
        
        STATE$Ye0[[deep]][,combination] = expit(Ze0_new)
      }
    }
    
    # SAMPLER STORE 
    
    if( ((max.iter - iter) %% thinning) == 0 & iter > burn.in ){
      
      count = count + 1
      
      Gs =  1/2^(14-6)
      
      for(deep in 3:6){
        G_deep = NULL
        for(combination in 1:2^(deep-1)){
          Ye0 = expit(c(rmvnorm(1,  mean = PRIORS$GP[[deep]][[combination]]$Ex,
                                   sigma = PRIORS$GP[[deep]][[combination]]$Var)))
          Ye1 = 1 - Ye0
          G_deep = cbind(G_deep, matrix(rep(Ye0, 2^(14-deep)), ncol = 2^(14-deep), nrow = 50),
                                 matrix(rep(Ye1, 2^(14-deep)), ncol = 2^(14-deep), nrow = 50))
          
        }
        Gs = G_deep * Gs
      }
      
      for(study in 1:n){
        Yi0   = STATE$Ye0[[1]][study,1]
        Yi00  = STATE$Ye0[[2]][study,1]
        Yi10  = STATE$Ye0[[2]][study,2]
        
        Gi = Gs[study,] * rep(c(Yi0*Yi00, Yi0*(1-Yi00),  (1-Yi0)*Yi10,  (1-Yi0)*(1-Yi10)), each = 2^12)
        STATE$medians[study] = G0inv(PRIORS$Pi[[14]][study, which.max(cumsum(Gi)>.5)])
      }
      
      MCMC_list[[count]] = STATE
    }
    
    # PROGRESS BAR
    
    if(step[n.cat]==iter){
      t1 = Sys.time()
      iter.time = difftime(t1,t0,units='secs')/iter
      time.to.end = as.numeric(iter.time)*(max.iter-iter)
      if(time.to.end>=3600*24){
        time.print = "+ 1d"
      }else if(time.to.end>3600){
        time.print =  paste0(round((time.to.end)/3600,1)," h")
      }else if(time.to.end>60){
        time.print =  paste0(round((time.to.end)/60,1)," m")
      }else{        
        time.print =  paste0(round((time.to.end),1)," s")
      }
      bar_print = paste0("|",paste0(rep('=', n.cat), collapse = ''), paste0(rep(' ', barwidth-n.cat), collapse = ''),
                         "|",paste0(rep(' ', 2-len[n.cat]), collapse = ''),perc[n.cat],"%"," (",time.print,")")
      bar_print = paste0(bar_print, paste0(rep(' ', abs(width - nchar(bar_print))+1), collapse = ''),collapse = '')
      
      cat('\r',bar_print)
      n.cat= n.cat + 1
    }
  }
  tend  = Sys.time()  
  
  i1 = seq(1,49,2)
  i2 = seq(2,50,2)
  i3_obs = seq(from = 1, to = 49, by = 2)
  
  
  
  v.meta    = rep(1/20 + 1/20,25)
  y.meta    = log(data[i1,6] / data[i2,6])
  

  
  data.meta = data.frame(logHR  = y.meta, sd.HR = sqrt(v.meta), 
                         agent  = factor(data$agent[i3_obs]),
                         tumor  = factor(data$tumor[i3_obs]),
                         marker = factor(data$marker[i3_obs]))
  
 
  saveRDS(file = paste0("PTGP_simulation/",Names),
          list("seed"      = seed,
               "PRIORS"    = PRIORS,
               "data"      = data,
               "data.meta" = data.meta,
               "chain"     = MCMC_list,
               "time"      = tend - t0))
  
}

q(save = "no")
