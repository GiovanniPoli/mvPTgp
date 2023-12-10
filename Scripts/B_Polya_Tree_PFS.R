### Packages & Functions            -----
library(readr)
library(readxl)
library(stringr)
library(mvtnorm)
library(cowplot)
library(reshape2)
library(survival)
library(BayesLogit)
library(metafor)
library(meta)
library(mvmeta)
library(gdata)
library(varhandle)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)
library(MASS)
library(HDInterval)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggridges)
library(latex2exp)
library(ggforce)

My_gg.Heatmap.cor = function(Matrix, from = NULL, to = NULL, x_lab  = NULL, y_lab = NULL, gg_title=NULL){
  
  if(is.null(from)){from = min(c(Matrix))} 
  if(is.null(to)){to = max(c(Matrix))}
  
  label = colnames(Matrix)
  
  xlim = ncol = ncol(Matrix)
  ylim = nrow = nrow(Matrix)
  
  df.heatmap = data.frame("value" = c(Matrix))
  
  df.heatmap$x         = rep(1:ncol, each  = nrow)
  df.heatmap$y         = rep(1:nrow, times = ncol)
  
  plot = ggplot(df.heatmap, aes(x=x, y=y, fill=value)) +
    geom_tile(col="grey") +
    scale_fill_viridis_c("", limits = c(from,to)) +
    scale_x_continuous(expand = c(0,0), breaks = 1:ncol, labels = label) +
    scale_y_continuous(expand = c(0,0), breaks = 1:nrow, labels = label) +
    xlab(x_lab)+ylab(y_lab)+ggtitle(gg_title) + theme_minimal() +
    theme(text = element_text(family ="serif"),
          axis.text.y = element_text(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  plot
  
}


erfinv = function(x) qnorm((1+x)/2)/sqrt(2)
expit  = function(z) 1/(1+exp(-z))
logit  = function(p) if(all( p > 1e-16)){log(p)-log(1-p)}else{ log(p + 1e-16) - log( 1 - p - 1e-16) }

rPTmarg.med = function(correlated_levels, truncation, c = 1, G0 = NULL, G0inv = NULL){
  
  PT = list()
  
  if(is.null(G0) | is.null(G0inv)){
    G0inv = function(q) punif(q, min = 0, max = 1)
    G0    = function(p) punif(p, min = 0, max = 1)
  }
  
  Pi_u = sapply( 1:truncation, function(x) c(0, 1:(2^x)) /2^x)
  
  G = 1
  for(deep in 1:truncation){
    if(deep <=  correlated_levels ){
      Ye0 = expit(rnorm( 
        n = 2^(deep-1), mean   = 0,          sd = sqrt(2 * trigamma(deep^2*c))
      ))
    }else{
      Ye0 = rbeta( n = 2^(deep-1), shape1 = deep^2*c, shape2 = deep^2*c)
    }
    
    Ye = c(sapply(Ye0, function(p) c(p, 1-p)))
    G  = G * c(sapply(Ye, function(x) rep(x, each = 2^truncation / 2^deep  )))
    PT[[deep]] = list( "Be" = G0inv( Pi_u[[deep]] ),
                       "Ye" = Ye )
  }
  
  return( list("tree" = PT, "G" = G))
}


rGPmvPT_v1.vanilla   = function(n,    correlated_levels, truncation, correlation_matrix,  c = 1, G0inv = NULL){
  
  mvPT = vector("list", truncation)
  
  if(is.null(G0inv)){
    G0inv = function(q) punif(q, min = 0, max = 1)
  }
  
  Pi = lapply(1:truncation, function(deep){
    seq(from = 0, to = 1, length.out = 2^deep+1)
  })
  
  for(deep in 1:truncation){
    if(deep <=  correlated_levels){
      Ye0 = t(expit( sqrt(2*trigamma(deep^2* c)) * 
                       rmvnorm(n     = 2^(deep-1), 
                               mean  = rep(0,n),
                               sigma = correlation_matrix)))
    }else{
      Ye0 = matrix(rbeta(n      = n * 2^(deep-1), 
                         shape1 = deep^2*c,
                         shape2 = deep^2*c),
                   ncol   =  2^(deep-1),
                   nrow   = n)
    }
    mvPT[[deep]] = list("Be" = G0inv(Pi[[deep]]),
                        "Ye" = t(apply(Ye0,1,  function(r) sapply(r,
                                                                  function(p) c(p,1-p)))))
    
  }
  mvPT
} 
rGPmvPT_v2.quantiles = function(n, Q, correlated_levels, truncation, GP_parameters_list,  c = 1, G0 = NULL, G0inv = NULL){
  
  mvPT = list()
  
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
      
      layer = cbind( cbind(0,sapply( 1:2^(d-2), function(i)           B00 * i/2^(d-2))),
                     sapply( 1:2^(d-2), function(i) B00 +             B01 * i/2^(d-2)),
                     sapply( 1:2^(d-2), function(i) B00 + B01 +       B10 * i/2^(d-2)),
                     sapply( 1:2^(d-2), function(i) B00 + B01 + B10 + B11 * i/2^(d-2)))
    } 
    layer})
  
  Y0  = expit(c(rmvnorm( 1, 
                         mean  = GP_parameters_list[[1]][[1]]$Ex, 
                         sigma = GP_parameters_list[[1]][[1]]$Var
  )))
  Y00 = expit(c(rmvnorm( 1, 
                         mean  = GP_parameters_list[[2]][[1]]$Ex,
                         sigma = GP_parameters_list[[1]][[1]]$Var
  )))
  Y10 = expit(c(rmvnorm( 1,
                         mean  = GP_parameters_list[[2]][[2]]$Ex,
                         sigma = GP_parameters_list[[2]][[2]]$Var
  )))
  Ye0 = t(rbind(Y00,Y10))
  
  mvPT[[1]] = list("Be" = G0inv(Pi[[1]]), 
                   "Ye" = t(sapply(Y0,  function(p) c(p,1-p) )))
  mvPT[[2]] = list("Be" = G0inv(Pi[[2]]),
                   "Ye" = t(apply(Ye0,1, function(r) sapply(r, function(p) c(p,1-p)))))
  
  if(J>=3){
    for(deep in 3:J){
      if(deep <=  D){
        Ye0 = NULL
        for(Nn in 1: 2^(deep-1))
          Ye0 = cbind(Ye0, c(expit(rmvnorm( 1, 
                                            mean  = GP_parameters_list[[deep]][[Nn]]$Ex, 
                                            sigma = GP_parameters_list[[deep]][[Nn]]$Var))))
      }else{
        Ye0 = matrix(rbeta( n = n * 2^(deep-1), 
                            shape1 = deep^2*c,
                            shape2 = deep^2*c), ncol =  2^(deep-1), nrow = n)
      }
      mvPT[[deep]] = list("Be" = G0inv(Pi[[deep]]),
                          "Ye" = t(apply(Ye0,1, function(r) sapply(r, function(p) c(p,1-p)))))
    }
  }
  mvPT
}

GP_parameters.gen_par = function(B, Q, correlated_levels, star_model_deep, correlation_star_matrix, c = 1,
                                 G0 = NULL, G0inv = NULL, n.cores = 20){
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
                    
                    expit  = function(z) 1/(1+exp(-z))
                    logit  = function(p) if(all(p>0)){log(p)-log(1-p)}else{ log(p + 1e-16) - log( 1 - p - 1e-16) }
                    
                    
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
  
  Ret = lapply(1:correlated_levels, function(deep) lapply(1:2^(deep-1), function(comb){
    
    Ye   = t(sapply(store, function(obj) obj[[deep]][[comb]]))
    data =    logit(Ye)
    Ex_Y = colMeans(Ye)
    
    Ex_Z   = colMeans(data)
    Var_Z  = var(data)
    Cor_Z  = cor(data)
    
    cat("\r\tAlmost done! Inverting covariances matrix...\t",paste0("[", 2^(deep-1)+comb,"/",2^correlated_levels,"]"))
    list("Ex" = Ex_Z, "Var" = Var_Z, "Cor" = Cor_Z, "Y_Ex" = Ex_Y)
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

plot.mvPT_picewise  = function(mvtree, truncation, quantile_unique = TRUE, plot = TRUE){
  
  df_1 = NULL
  df_2 = NULL
  n = nrow(mvtree[[1]]$Ye)
  for(id in 1:n){
    probs = apply(sapply(1:truncation, function(x) rep(mvtree[[x]]$Ye[id,], each = 2^(truncation-x))),
                  1, function(y) prod(y) )
    
    if(quantile_unique){
      l = mvtree[[truncation]]$Be[1:(2^truncation-1)] 
      r = mvtree[[truncation]]$Be[2:(2^truncation)]    
      
      width = r - l
      
      df   = data.frame(ys = rep(probs[1:(2^truncation-1)])/width,
                        ye = rep(probs[1:(2^truncation-1)])/width,
                        xs = mvtree[[truncation]]$Be[1:(2^truncation-1)], 
                        xe = mvtree[[truncation]]$Be[2:(2^truncation)])
      
      df$col = id
      
      df.2 = data.frame(ye = c(0, probs[1:(2^truncation-2)])/c(1,width[1:(2^truncation-2)]),
                        ys = probs[1:(2^truncation-1)]/width,
                        xs = mvtree[[truncation]]$Be[1:(2^truncation-1)],
                        xe = mvtree[[truncation]]$Be[1:(2^truncation-1)])
      df.2$col = id
      
      df_1 = rbind(df_1, df)
      df_2 = rbind(df_2, df.2)
      
      
    }else{
      l = mvtree[[truncation]]$Be[id, 1:(2^truncation-1)] 
      r = mvtree[[truncation]]$Be[id, 2:(2^truncation)]
      
      width = r - l
      
      df   = data.frame(ys = rep(probs[1:(2^truncation-1)])/width,
                        ye = rep(probs[1:(2^truncation-1)])/width,
                        xs = mvtree[[truncation]]$Be[id, 1:(2^truncation-1)], 
                        xe = mvtree[[truncation]]$Be[id, 2:(2^truncation)  ])
      
      df$col = id
      
      df.2 = data.frame(ye = c(0, probs[1:(2^truncation-2)])/c(1,width[1:(2^truncation-2)]),
                        ys = probs[1:(2^truncation-1)]/width,
                        xe = mvtree[[truncation]]$Be[id, 1:(2^truncation-1)],
                        xs = mvtree[[truncation]]$Be[id, 1:(2^truncation-1)])
      df.2$col = id
      
      df_1 = rbind(df_1, df)
      df_2 = rbind(df_2, df.2)
    }
  }
  df_1$col = as.factor(df_1$col)
  df_2$col = as.factor(df_2$col)
  
  if(plot){
    plot = ggplot() +
      geom_segment(data = df_1,   aes(x = xs, y = ys, xend = xe, yend = ye, col= col)) +
      geom_segment(data = df_2,   aes(x = xs, y = ys, xend = xe, yend = ye, col= col)) +
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() + xlab("Time (month)") + ylab("density") +
      scale_color_viridis_d() +
      theme(text = element_text(family = "serif"),legend.position = "none")
    plot
    return(plot)
  }else{
    return(list(df_density = df_1, df_bar = df_2))
  }
}
plot.mvPT_survival  = function(mvtree, truncation, quantile_unique = TRUE, plot = TRUE){
  
  df_1 = NULL
  df_2 = NULL
  n = nrow(mvtree[[1]]$Ye)
  
  for(id in 1:n){
    probs = apply(sapply(1:truncation, function(x) rep(mvtree[[x]]$Ye[id,], each = 2^(truncation-x))),
                  1, function(y) prod(y) )
    
    Surv_fun = (1 - cumsum(probs))
    
    if(quantile_unique){
      
      df   = data.frame(ys = Surv_fun,
                        ye = Surv_fun,
                        xs = mvtree[[truncation]]$Be[1:(2^truncation)], 
                        xe = mvtree[[truncation]]$Be[2:(2^truncation+1)])
      
      df$col = id
      
      df.2 = data.frame(ye = c(1,Surv_fun),
                        ys = c(Surv_fun,0),
                        xs = c(mvtree[[truncation]]$Be[1:(2^truncation)],Inf),
                        xe = c(mvtree[[truncation]]$Be[1:(2^truncation)],Inf))
      df.2$col = id
      
      df_1 = rbind(df_1, df)
      df_2 = rbind(df_2, df.2)
      
      
    }else{
      l = mvtree[[truncation]]$Be[id, 1:(2^truncation)] 
      r = mvtree[[truncation]]$Be[id, 2:(2^truncation+1)]
      
      width = r - l
      
      df   = data.frame(ys = Surv_fun,
                        ye = Surv_fun,
                        xs = mvtree[[truncation]]$Be[id, 1:(2^truncation)], 
                        xe = mvtree[[truncation]]$Be[id, 2:(2^truncation+1)])
      
      df$col = id
      
      df.2 = data.frame(ye = c(1,Surv_fun),
                        ys = c(Surv_fun,0),
                        xe = c(mvtree[[truncation]]$Be[id, 1:(2^truncation)],Inf),
                        xs = c(mvtree[[truncation]]$Be[id, 1:(2^truncation)],Inf))
      df.2$col = id
      
      df_1 = rbind(df_1, df)
      df_2 = rbind(df_2, df.2)
    }
  }
  df_1$col = as.factor(df_1$col)
  df_2$col = as.factor(df_2$col)
  
  if(plot){
    plot = ggplot() +
      geom_segment(data = df_1,   aes(x = xs, y = ys, xend = xe, yend = ye, col= col)) +
      geom_segment(data = df_2,   aes(x = xs, y = ys, xend = xe, yend = ye, col= col)) +
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() + xlab("Time (month)") + ylab("density") +
      scale_color_viridis_d() +
      theme(text = element_text(family = "serif"),legend.position = "none")
    plot
    return(plot)
  }else{
    return(list(df_density = df_1, df_bar = df_2))
  }
}

### Import data                     -----
# Metdata

metadata       = read_table("Dati/metadata.txt")
metadata_excel = read_excel("Dati/Feb12022.xlsx")

metadata_excel$study_labels = apply(metadata_excel , 1, function(x){
  if(!is.na(x[11])){
    strs = unlist(stringr::str_split(x[11], ",",n = 4))
    strs = sapply(strs, function(x) str_remove_all(x, "[1234567890#]"), USE.NAMES = FALSE)
    strs = sapply(strs, function(x) str_squish(x), USE.NAMES = FALSE)
    strs = sapply(strs, FUN = function(x) str_split(x," "), USE.NAMES = FALSE)
    year = str_split(x[13], " ")[[1]][1]
    if(length(strs) == 1){
      string = paste0(strs[[1]][length(strs[[1]])]," (",year,")")
    }else if(length(strs) == 2){
      string = paste0(strs[[1]][length(strs[[1]])]," and ",
                      strs[[2]][length(strs[[2]])]," (",year,")")
    }else{
      string = paste0(strs[[1]][length(strs[[1]])]," et al. (",year,")")
    }
    string
  }else{
    return(NA)
  }})

# remove duplicate od Doi et al. 2019
metadata = metadata[metadata$k != "167",]

# OS
# osr_import =   read_table("C:/Users/jovi1/Dropbox/giovanniFirenze/data/osr_import.txt")
# osr_positive  = cbind(osr_import[,c(1,2,4:6)] ,"positive")
# osr_negative  = cbind(osr_import[,c(1,3,7:9)] ,"negative")
# colnames(osr_negative)    = colnames(osr_positive) = c("k","n","median","lower","upper","marker")
# osr_data  = rbind(osr_positive,  osr_negative)     

# PFS
pfsr_import   = read_table("Dati/pfsr_import.txt")
# remove duplicate of Doi et al. 2019
pfsr_import   = pfsr_import[pfsr_import$k != "167",]

pfsr_positive = cbind(pfsr_import[,c(1,2,4:6)],"positive")
pfsr_negative = cbind(pfsr_import[,c(1,3,7:9)],"negative")

colnames(pfsr_negative) = colnames(pfsr_positive)  = c("k","n","median","lower","upper","marker")

pfsr_data = rbind(pfsr_positive, pfsr_negative) 

pfsr_data = pfsr_data[rep(1:nrow(pfsr_positive),    each = 2) + 
                        rep(c(0,nrow(pfsr_positive)), nrow(pfsr_positive)),]

pfsr_data = pfsr_data[!is.na(pfsr_data$n),]     # NA remove for n
pfsr_data[3:5][pfsr_data[3:5] ==  99.9 | pfsr_data[3:5] == -99.9 | pfsr_data[3:5] ==  0 ] = NA
pfsr_data[3:5][pfsr_data[3:5] == -99.9] = NA
pfsr_data$statistics = NA
pfsr_data$statistics   = ifelse(is.na(pfsr_data$lower)  & is.na(pfsr_data$upper),  "none",  pfsr_data$statistics )
pfsr_data$statistics   = ifelse(!is.na(pfsr_data$lower) & is.na(pfsr_data$upper),  "lower", pfsr_data$statistics )
pfsr_data$statistics   = ifelse(!is.na(pfsr_data$lower) & !is.na(pfsr_data$upper), "both",  pfsr_data$statistics )
pfsr_data$censinfo     = "none"
pfsr_data$event        = NA
pfsr_data$ci           = .95
pfsr_data$ci.type      = "log"
#plain_study = c("12","27","95","137","142","149","166","167","153")
plain_study = c("12","27","95","137","142","149","166","153")

pfsr_data$ci.type[pfsr_data$k %in% plain_study] = "plain"
pfsr_data$cens.pattern = NA
rownames(pfsr_data)  = paste0("cohort.",1:86)

## 0) STUDY CORRECTIN
## Study 12 positive marker impossible upper limit with 11 observations  // Sutdy 12 just have a "plain" style
# pfsr_data$upper      [pfsr_data$k == "12" & pfsr_data$marker == "positive"] = NA
# pfsr_data$statistics [pfsr_data$k == "12" & pfsr_data$marker == "positive"] = "lower"
## Study 27, info ok, ci = 80%
pfsr_data$ci[pfsr_data$k == "27"] = .8
## Study 38 Upper limit not reached (our mistake)  // coded xls 12+ as 12
pfsr_data$upper     [pfsr_data$k == "38"] = NA
pfsr_data$statistics[pfsr_data$k == "38"] = "lower"
## Study 45, n = 27 in negative  //  must check it
# pfsr_data$n         [pfsr_data$k == "45" & pfsr_data$marker == "negative"] = 19
## Study 165, n = 27 in negative // must check it
pfsr_data$n         [pfsr_data$k == "165" & pfsr_data$marker == "negative"] = 27


## 1) ADD KNOW NUMBER OF EVENTS 
# k = 7, from fig .2
pfsr_data$censinfo[pfsr_data$k == "7"] = "events"
pfsr_data$event   [pfsr_data$k == "7"] = c(147-32,109-12) 
# k = 37, from appendix
pfsr_data$censinfo[pfsr_data$k == "37"] = "events"
pfsr_data$event   [pfsr_data$k == "37"] = c(28,37, 57,21)
# k = 137, from appendix, if corect events = n-1
pfsr_data$censinfo[pfsr_data$k == "137"] = "events"
pfsr_data$event   [pfsr_data$k == "137"] = c(45,11)
# k = 149, from appendix
pfsr_data$censinfo[pfsr_data$k == "149"] = "events"
pfsr_data$event   [pfsr_data$k == "149"] = c(11,11)
# k = 153, from appendix, if corect events = n-1
pfsr_data$censinfo[pfsr_data$k == "153"] = "events"
pfsr_data$event   [pfsr_data$k == "153"] = c(13,37)
# k = 170, from appendix plot
pfsr_data$censinfo[pfsr_data$k == "170"] = c("events","events","none","none")
pfsr_data$event   [pfsr_data$k == "170"] = c(18,34,NA,NA)

## 2) ADD KNOW CENSURE PATTERN 
# k = 2 from fig .3
pfsr_data$censinfo[pfsr_data$k == "2"] = "complete"
pfsr_data$cens.pattern[pfsr_data$k == "2" & pfsr_data$marker == "positive" ] = list(rep(1,27))
pfsr_data$cens.pattern[pfsr_data$k == "2" & pfsr_data$marker == "negative" ] = list(rep(1,24))
# k = 12 from word online appendix  KM,  
pfsr_data$censinfo[pfsr_data$k == 12] = "complete"
pfsr_data$cens.pattern[pfsr_data$k == "12" & pfsr_data$marker == "positive" ] = list(c(rep(1,10),0))
pfsr_data$cens.pattern[pfsr_data$k == "12" & pfsr_data$marker == "negative" ] = list(rep(1,27))
# k = 100 from Km, no censure,  
pfsr_data$censinfo[pfsr_data$k == "100"] = "complete"
pfsr_data$cens.pattern[pfsr_data$k == "100" & pfsr_data$marker == "positive" ] = list(rep(1,12))
pfsr_data$cens.pattern[pfsr_data$k == "100" & pfsr_data$marker == "negative" ] = list(rep(1,12))

## Study with errors

pfsr_data = pfsr_data[c(1:86)[!c(c(1:86) %in% c(13,14))],]
pfsr_data["cohort.31",] = data.frame("k" = 73,"n" = 20, "median" = 9,"lower" = NA, "upper" = NA, "marker" = "positive",
                                     "statistics" = "none", "censinfo" = "none","event" =  NA, "ci" = 0.95, 
                                     "ci.type"    = "log",  "cens.pattern" = NA)
pfsr_data["cohort.32",] = data.frame("k" = 73,"n" = 19, "median" = 1,"lower" = NA, "upper" = NA, "marker" = "negative",
                                     "statistics" = "none", "censinfo" = "none","event" =  NA, "ci" = 0.95, 
                                     "ci.type"    = "log",  "cens.pattern" = NA)

## LABELS
for(i in 1:length(metadata_excel$trialNumbering)){
  if(is.na(metadata_excel$trialNumbering[i])){
    metadata_excel$trialNumbering[i] = old
  }else{
    old = metadata_excel$trialNumbering[i]
  }
}
index  = !is.na(metadata_excel$study_labels)
labels = metadata_excel$study_labels[index]
names(labels) =  as.character(metadata_excel$trialNumbering[index])
pfsr_data$label = labels[as.character(pfsr_data$k)]

labels = metadata_excel$study_labels
labels = labels[as.character(metadata$k)]

index  = !is.na(metadata_excel$study_labels)
labels = metadata_excel$study_labels[index]
names(labels) =  as.character(metadata_excel$trialNumbering[index])
pfsr_data$label  = labels[as.character(pfsr_data$k)]
pfsr_data$label2 = NA

start = 1

for(i in 1:(nrow(pfsr_data)/2) ){
  ix1 = i*2-1
  ix2 = i*2
  
  pfsr_data$label2[ix1] = paste0("[+|",start,"] ",pfsr_data$label[ix1])
  pfsr_data$label2[ix2] = paste0("[-|",start+1,"] ",pfsr_data$label[ix1])
  
  
  if(pfsr_data$label[ix2] == pfsr_data$label[ix2+1] & ix2 < nrow(pfsr_data)){
    start = start + 2
  }else{
    start = 1
  }
}


pfsr_data$label[pfsr_data$k %in% metadata$k[metadata$agent == "pembrolizumab" ]]

remove(i,index, ix1, ix2, labels, old, start, pfsr_positive, pfsr_negative, pfsr_import,
       metadata_excel,plain_study)

### Star Model (Prior of ref model) -----


g0     = function(x, scale = 3.5){ exp(log(2 * scale) - log(pi * { x * x + scale * scale }))}
G0     = function(q, scale = 3.5){ {2/pi} * atan(q/scale) }
G0inv  = function(p, scale = 3.5){ scale * tan({ pi * p }/2) }

weights = c(2,1,2,2,
            2,5,
            1,5,
            2,1,
            2,4,
            1,2,11,
            1)

tumors = c("NSCLC","breast","melanoma","other",
           "melanoma","other",
           "NSCLC","other",
           "NSCLC","other",
           "melanoma","other",
           "NSCLC","breast","other",
           "other")

agents = c("atezolizumab","atezolizumab","atezolizumab","atezolizumab",
           "avelumab","avelumab",
           "durvalumab","durvalumab",
           "ipilimumab/nivolumab","ipilimumab/nivolumab",
           "nivolumab","nivolumab",
           "pembrolizumab", "pembrolizumab", "pembrolizumab",
           "pembrolizumab-or-nivolumab")

metadata_future = data.frame(cbind(rep(1:16, each = 2) + 300,
                                   rep(tumors, each = 2),
                                   rep(agents, each = 2),
                                   rep(c("positive","negative"),16),
                                   "monotherapy","≥2","1"))
n_futures_study  = 32
n_observed_study = nrow(pfsr_data)


K.star = matrix(NA, ncol = n_observed_study + n_futures_study,
                nrow = n_observed_study + n_futures_study)




colnames(metadata_future) = c("k","tumor","agent","marker","mono","firstline","phase")

colnames(K.star) = rownames(K.star) = c(pfsr_data$label2, paste0(metadata_future$k,".",paste0(metadata_future$marker)))

for(row in 1:(n_observed_study+n_futures_study-1)){
  
  for(col in (row+1):(n_futures_study+n_observed_study)){
    
    if(row <= n_observed_study & col <= n_observed_study){
      
      K.star[row,col] = K.star[col,row] = 0 + 
        ifelse(pfsr_data$marker[row] == pfsr_data$marker[col],                                                           1/2,0) + # same marker
        ifelse(metadata$firstline[metadata$k == pfsr_data$k[row]] == metadata$firstline[metadata$k == pfsr_data$k[col]], 1/2, 0) + # same lines
        ifelse(metadata$phase[metadata$k == pfsr_data$k[row]]     == metadata$phase[metadata$k == pfsr_data$k[col]],     1/2, 0) + # same phase
        ifelse(metadata$mono[metadata$k == pfsr_data$k[row]]     ==  metadata$mono[metadata$k == pfsr_data$k[col]],      1/2, 0) + # same theraphy
        ifelse(metadata$agent[metadata$k == pfsr_data$k[row]]     == metadata$agent[metadata$k == pfsr_data$k[col]],       2, 0) + # same agent
        ifelse(metadata$tumor[metadata$k == pfsr_data$k[row]]     == metadata$tumor[metadata$k == pfsr_data$k[col]] &
                 metadata$tumor[metadata$k == pfsr_data$k[row]]     != "other", 2, 0)                                              # same tumor
      
      K.star[row,col] = K.star[col,row] = ifelse(pfsr_data$marker[row] != pfsr_data$marker[col], 0 ,K.star[row,col])  # different marker correction
      
      K.star[row,col] = K.star[col,row] = ifelse(pfsr_data$k[row] == pfsr_data$k[col],  K.star[row,col] + 1,
                                                 K.star[row,col])  # same study
      
    }
    
    if(row <= n_observed_study & col >  n_observed_study){
      row. = col - n_observed_study
      K.star[row,col] = K.star[col,row] = 0 + 
        ifelse(metadata_future$marker[row.]    == pfsr_data$marker[row],                                  1/2, 0) + # same marker
        ifelse(metadata_future$firstline[row.] == metadata$firstline[metadata$k == pfsr_data$k[row]],     1/2, 0) + # same lines
        ifelse(metadata_future$phase[row.]     ==     metadata$phase[metadata$k == pfsr_data$k[row]],     1/2, 0) + # same phase
        ifelse( metadata_future$mono[row.]     ==      metadata$mono[metadata$k == pfsr_data$k[row]],     1/2, 0) + # same theraphy
        ifelse(metadata_future$agent[row.]     ==     metadata$agent[metadata$k == pfsr_data$k[row]],     2, 0) + # same agent
        ifelse(metadata_future$tumor[row.]     ==     metadata$tumor[metadata$k == pfsr_data$k[row]] &
                 metadata$tumor[metadata$k == pfsr_data$k[row]] != "other"   ,                              2, 0)                            # same tumor
      
      
      K.star[row,col] = K.star[col,row] = ifelse(metadata_future$marker[row.] != pfsr_data$marker[row], 0 ,K.star[row,col])  # different marker
      K.star[row,col] = K.star[col,row] = ifelse(metadata_future$k[row.] == pfsr_data$k[row], K.star[row,col] + 1, 
                                                 K.star[row,col])               # same study
      
    }
    
    if(row > n_observed_study & col > n_observed_study)  {
      row. = col - n_observed_study
      col. = row - n_observed_study
      
      K.star[row,col] = K.star[col,row] = 0 + 
        ifelse(metadata_future$marker[row.]    == metadata_future$marker[col.],    1/2, 0) + # same marker
        ifelse(metadata_future$firstline[row.] == metadata_future$firstline[col.], 1/2, 0) + # same lines
        ifelse(metadata_future$phase[row.]     == metadata_future$phase[col.],     1/2, 0) + # same phase
        ifelse(metadata_future$mono[row.]      == metadata_future$mono[col.],      1/2, 0) + # same theraphy
        ifelse(metadata_future$agent[row.]     == metadata_future$agent[col.],     2, 0) + # same agent
        ifelse(metadata_future$tumor[row.]     == metadata_future$tumor[col.] &
                 metadata_future$tumor[row.]     != "other",     2, 0)   # same tumor  
      
      
      K.star[row,col] = K.star[col,row] = ifelse(metadata_future$marker[row.] != metadata_future$marker[col.],
                                                 0,
                                                 K.star[row,col])  
      
      K.star[row,col] = K.star[col,row] =  ifelse(metadata_future$k[row.] == metadata_future$k[col.],
                                                  K.star[row,col] + 1, 
                                                  K.star[row,col])  # same study
      
    }
  }
}
max(K.star, na.rm = TRUE)

K.star = K.star / 8 
max(K.star, na.rm = TRUE)

diag(K.star) = 1

K.plot = K.star
colnames(K.plot) = rownames(K.plot) = NULL
R.heatmap = My_gg.Heatmap.cor(K.plot[ 
  c(seq(from = 1, to = 115, by =2), seq(from = 2, to = 116, by =2)),
  c(seq(from = 1, to = 116, by =2), seq(from = 2, to = 116, by =2)) ])
R.heatmap

Qs = pfsr_data[,c(4,3,5)]
Qs$lower[is.na(Qs$lower)] = round(G0inv( (    G0(pfsr_data[,3])) /2)[is.na(Qs$lower)],2)
Qs$upper[is.na(Qs$upper)] = round(G0inv( G0(pfsr_data[,3]) + (1 - G0(pfsr_data[,3])) /2)[is.na(Qs$upper)],2)
Qs = rbind(as.matrix(Qs),
           matrix( c(median(Qs$lower),
                     median(Qs$median),
                     median(Qs$upper)), byrow= TRUE, ncol = 3, nrow = n_futures_study)  )

PAR = readRDS("Star_PFSr_02112023.rds")


### PRIORS                          -----

GP_prior     = list() # Need: Yex for Prior, Zex, Var Prec for update, Var

for(deep in 1:6){
  l  = list()
  l2 = list()
  for(comb in 1:2^(deep-1)){
    it  = PAR[[deep]][[comb]]
    
    inv = chol2inv(chol(it$Var[ 1:n_observed_study, 1:n_observed_study]))
    
    l[[comb]] = list("Z_Ex"   = it$Ex[1:n_observed_study],
                     "Z_Prec"   = inv,
                     "Z_Var"    = it$Var[1:n_observed_study, 
                                         1:n_observed_study],
                     "Y_Ex"    = it$Y_Ex[1:n_observed_study])
  }
  GP_prior[[deep]] = l
}

MCMC_data = data.frame(Qs[1:n_observed_study,],
                       "si" =pfsr_data$statistics[1:n_observed_study])

InvCDF.exp    = function(n, a, b = Inf, rate = .01) -log(1-runif(n, 
                                                                 min = pexp(a, rate = rate), 
                                                                 max = pexp(b, rate = rate))) / rate

InvCDF.G0     = function(n, a, b = Inf) G0inv(runif(n, min = G0(a), max =G0(b)))

KM_start =list()

PRIORS = list("GP"               = GP_prior,
              "D"                = 6,
              "Qs"               = Qs,
              "observed"         = n_observed_study,
              "predictive"       = n_futures_study,
              "Pi"               = Pi.gen(Q = Qs[1:n_observed_study,], truncation =  14, G0 = G0,G0inv = G0inv),
              "G0"               = G0, 
              "G0inv"            = G0inv,
              "rH0"              = function(n) rexp(n = n, rate = .01),
              "rH0tr"            = function(n,a,b) InvCDF.exp(n = n, a = a, b = b, rate = .05),
              "H0"               = function(q) pexp(q = q, rate = .05))

### STARTING VALUES                 -----

for(cohort in 1:n_observed_study){
  
  if(pfsr_data$censinfo[cohort] == "none"){
    delta = rep(1, pfsr_data$n[cohort])
  }else if(pfsr_data$censinfo[cohort] == "events"){
    delta = c(rep(1,pfsr_data$event[[cohort]]),
              rep(0,pfsr_data$n[[cohort]] - pfsr_data$event[[cohort]]))
  }else if(pfsr_data$censinfo[cohort] == "complete"){
    delta = pfsr_data$cens.pattern[[cohort]]
  }
  
  if(cohort == 13) delta = c(rep(1,6),0,0,0)
  if(cohort == 19) delta = c(rep(1,15),rep(0,4))
  if(cohort == 28) delta = c(rep(1,13),0,0,0)
  
  df    = data.frame(t = 1:pfsr_data$n[cohort], d = delta)
  km    = survfit(Surv(t, d) ~ 1, data= df, conf.int  = pfsr_data$ci[cohort],
                  conf.type = pfsr_data$ci.type[cohort])
  q = c(which.max(km$lower <= .5),
        which.max(km$surv  <= .5),  
        which.max(km$upper <= .5))
  
  # if(pfsr_data$n[cohort]%%2 == 0) q[1] = df$t[which.min(df$t<q[1])- 1] 
  
  
  
  if(pfsr_data$statistics[cohort] == "both"){
    
    t1 = InvCDF.G0(length(df$t[df$t < q[1]]              ), a = 0, b = Qs[cohort,1])
    t2 = InvCDF.G0(length(df$t[df$t < q[2] & df$t > q[1]]), a = Qs[cohort,1], b = Qs[cohort,2])
    t3 = InvCDF.G0(length(df$t[df$t > q[2] & df$t < q[3]]), a = Qs[cohort,2], b = Qs[cohort,3]) 
    t4 = InvCDF.G0(length(df$t[df$t > q[3] & df$d != 0]  ), a = Qs[cohort,3])
    
    times    = c(sort(t1),Qs[cohort,1],sort(t2),Qs[cohort,2],sort(t3),Qs[cohort,3],sort(t4))
    censure  = InvCDF.exp(length(times), times)
    
    
    if(pfsr_data$n[cohort] - length(times) > 0){
      observed      = length(times)
      censored      = pfsr_data$n[cohort] - length(times)
      last_obs_time = times[length(times)]
      censure       = c(censure, InvCDF.exp(censored, a = last_obs_time))
      times         = c(times,   InvCDF.G0(n = censored,
                                           a = censure[(observed+1):(observed+censored)]))
    }
    
  }else if(pfsr_data$statistics[cohort] == "lower"){
    
    t1 = InvCDF.G0(length(df$t[df$t < q[1]]              ), a = 0, b = Qs[cohort,1])
    t2 = InvCDF.G0(length(df$t[df$t < q[2] & df$t > q[1]]), a = Qs[cohort,1], b = Qs[cohort,2])
    t3 = InvCDF.G0(length(df$t[df$t > q[2] & df$d != 0 ]),  a = Qs[cohort,2]) 
    
    times    = c(sort(t1),Qs[cohort,1],sort(t2),Qs[cohort,2],sort(t3))
    censure  = InvCDF.exp(length(times), times, b = Inf)
    
    if(pfsr_data$n[cohort] - length(times) > 0){
      observed      = length(times)
      censored      = pfsr_data$n[cohort] - length(times)
      last_obs_time = times[length(times)]
      censure       = c(censure, InvCDF.exp(censored, a = last_obs_time))
      times         = c(times,   InvCDF.G0(n = censored,
                                           a = censure[(observed+1):(observed+censored)]))
    }
    
  }else if(pfsr_data$statistics[cohort] == "none"){
    
    t1 = InvCDF.G0(length(df$t[df$t < q[2]]              ), a = 0, b = Qs[cohort,2])
    t2 = InvCDF.G0(length(df$t[df$t > q[2] & df$d != 0 ]),  a = Qs[cohort,2]) 
    
    times    = c(sort(t1),Qs[cohort,2],sort(t2))
    censure  = InvCDF.exp(length(times), times, b = Inf)
    
    if(pfsr_data$n[cohort] - length(times) > 0){
      observed      = length(times)
      censored      = pfsr_data$n[cohort] - length(times)
      last_obs_time = times[length(times)]
      censure       = c(censure, InvCDF.exp(censored, a = last_obs_time))
      times         = c(times,   InvCDF.G0(n = censored,
                                           a = censure[(observed+1):(observed+censored)]))
    }
  }
  cens_count = c(0, 0, 0, sum(delta==0))
  time_count = c(sum(times  < Qs[cohort,1]),
                 sum(times >= Qs[cohort,1] & times < Qs[cohort,2]),
                 sum(times >= Qs[cohort,2] & times < Qs[cohort,3]),
                 sum(times >= Qs[cohort,3] & delta !=0 ) )
  tobs    = pmin(times,censure)
  order   = order(tobs)
  tobs    = tobs[order]
  times   = times[order]
  censure = censure[order] 
  
  KM_start[[cohort]] = list("times"   = times,
                            "censure" = censure,
                            "t_obs"   = tobs,
                            "delta"   = delta)
}

starting_values = list("Y0"  = matrix(.5, ncol = 1, nrow =  n_observed_study),
                       "Ye0" = matrix(.5, ncol = 2, nrow =  n_observed_study),
                       "KM"  = KM_start)


samples.size_step = function(KM,Q){
  t(sapply(1:length(KM), function(ind){
    c(sum( KM[[ind]]$times  <  Q[ind,1]),
      sum( KM[[ind]]$times >=  Q[ind,1] & KM[[ind]]$times < Q[ind,2]),
      sum( KM[[ind]]$times >=  Q[ind,2] & KM[[ind]]$times < Q[ind,3]),
      sum( KM[[ind]]$times >=  Q[ind,3]))}))
}

### Gibbs                           -----

Meta_sampler.v1 = function(data, 
                           PRIORS,
                           Qs = Qs,
                           starting_values,
                           burn.in  = 10,
                           thinning = 50,
                           sample   = 1000,
                           seed = 1){
  
  MCMC_list = list()
  ## STARTING VALUES & VARIABLES
  set.seed(seed)
  
  count = 0
  
  observed_sample = PRIORS$observed
  
  ar1 = vector("numeric", observed_sample)
  ar2 = vector("numeric", observed_sample)
  
  STATE   = starting_values  
  STATE$N = samples.size_step(STATE$KM, Qs)
  STATE$medians = vector(mode = "numeric", length = observed_sample)
  
  max.iter = burn.in + sample * thinning
  
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
  
  t0    = Sys.time()
  
  for(iter in 1:max.iter){
    
    ## STEP I - SAMPLER Y0,Y00, Y01 | impossible empty set (always lower limit)
    ## 1.1) Y0 
    
    try   =  rowSums(STATE$N)
    sucs  =  STATE$N[,1] + STATE$N[,2]
    
    k     = sucs - try / 2
    w     = BayesLogit::rpg.devroye(observed_sample , try, logit(STATE$Y0))
    
    
    Omega     = diag(w)
    Sigma_b   = solve(Omega      +   PRIORS$GP[[1]][[1]]$Z_Prec)
    Mean_b    = c( Sigma_b%*%( k + c(PRIORS$GP[[1]][[1]]$Z_Prec%*%PRIORS$GP[[1]][[1]]$Z_Ex)) )
    
    STATE$Y0 = expit(c(rmvnorm(n = 1, mean = Mean_b, sigma = Sigma_b)))
    
    # 1.2) Y00
    
    try   =  STATE$N[,1] + STATE$N[,2]
    sucs  =  STATE$N[,1] 
    
    k     = sucs - try / 2
    w     = BayesLogit::rpg.devroye(observed_sample , try, logit( STATE$Ye0[,1] ))
    
    Omega     = diag(w)
    Sigma_b   = solve(Omega      +   PRIORS$GP[[2]][[1]]$Z_Prec)
    Mean_b    = c( Sigma_b%*%( k + c(PRIORS$GP[[2]][[1]]$Z_Prec%*%PRIORS$GP[[2]][[1]]$Z_Ex)) )
    
    STATE$Ye0[,1] = expit(c(rmvnorm(n = 1, mean = Mean_b, sigma = Sigma_b)))
    # 1.3) Y01 
    
    try   =  STATE$N[,3] + STATE$N[,4]
    sucs  =  STATE$N[,3] 
    
    k     = sucs - try / 2
    w     = BayesLogit::rpg.devroye(observed_sample , try, logit(STATE$Ye0[,2] ))
    
    
    Omega     = diag(w)
    Sigma_b   = solve(Omega      +   PRIORS$GP[[2]][[2]]$Z_Prec)
    Mean_b    = c( Sigma_b%*%( k + c(PRIORS$GP[[2]][[2]]$Z_Prec %*% PRIORS$GP[[2]][[2]]$Z_Ex)) )
    
    STATE$Ye0[,2] = expit(c(rmvnorm(n = 1, mean = Mean_b, sigma = Sigma_b)))
    
    ## STEP II - Kaplan Maier 
    # 2.0) Data Augmentation for G
    
    
    Yie = t(sapply(c(STATE$Y0), function(p) rep(c(p,1-p), each = 2^13))) 
    
    Yiee = cbind( t(sapply(STATE$Ye0[1:n_observed_study,1], function(p) rep(c(p,1-p), each = 2^12))) ,
                  t(sapply(STATE$Ye0[1:n_observed_study,2], function(p) rep(c(p,1-p), each = 2^12))))
    
    Gs = Yie * Yiee * 2^(-8) # Ex after lv 6
    
    ## Prior Sim for levels 3, ..., 6
    for(deep in 3:6){
      G_deep = NULL
      for(combination in 1:2^(deep-1)){
        Ye0 = expit(c(rmvnorm(1, mean = PRIORS$GP[[deep]][[combination]]$Z_Ex[1:n_observed_study],
                              sigma = PRIORS$GP[[deep]][[combination]]$Z_Var[1:n_observed_study, 1:n_observed_study])))
        Ye1 = 1 - Ye0
        G_deep = cbind(G_deep, matrix(rep(Ye0,2^(14-deep)), ncol = 2^(14-deep), nrow = n_observed_study),
                       matrix(rep(Ye1,2^(14-deep)), ncol = 2^(14-deep), nrow = n_observed_study))
        
      }
      Gs = G_deep * Gs
    }
    
    ## STEP III KM
    for(cohort in 1:observed_sample){
      
      if( data$censinfo[cohort] == "none" ){
        delta_step = sample(c("move","flip"),1)
      }else if( data$censinfo[cohort] == "events"){
        delta_step = "move"
      }else if(data$censinfo[cohort] == "complete"){
        delta_step = "none"
      }
      
      if(delta_step=="flip"){
        
        new_delta = STATE$KM[[cohort]]$delta
        new_times = STATE$KM[[cohort]]$t_obs
        
        who.to.flip = sample.int(n = data$n[cohort], size = 1)
        di          = 1 - new_delta[who.to.flip]
        
        lower = max(0,   STATE$KM[[cohort]]$t_obs[who.to.flip - 1])
        upper = min(Inf, STATE$KM[[cohort]]$t_obs[who.to.flip + 1], na.rm = TRUE)
        
        G_pos_1 = min( which.max( PRIORS$Pi[[14]][ cohort, ] > lower),2^14)
        G_pos_2 = min( which.max( PRIORS$Pi[[14]][ cohort, ] >= upper),2^14)
        G_pos_3 = min( which.max( PRIORS$Pi[[14]][ cohort, ] >  STATE$KM[[cohort]]$t_obs[who.to.flip]),2^14)  
        
        
        if(di == 1){
          
          t     = ifelse(G_pos_1 != G_pos_2,
                         sample(x = PRIORS$Pi[[14]][ cohort, G_pos_1:G_pos_2], prob = Gs[cohort, G_pos_1:G_pos_2 ], size = 1),
                         PRIORS$Pi[[14]][ cohort, G_pos_1])
          # ifelse for becouse if x is a scalar sample gives => 
          # incorrect number of probabilities
          
          cen   = PRIORS$rH0tr(1, t, Inf)
          
          p_t   = sum(Gs[cohort, G_pos_1:G_pos_2 ])
          p_c.t = 1 - PRIORS$H0(t)
          
          rvp_c   = PRIORS$H0(upper) - PRIORS$H0(lower)
          rvp_t.c = sum(Gs[cohort, G_pos_3 : 2^14 ])
          
          alpha = min(1, p_t * p_c.t / rvp_c / rvp_t.c)
          
        }else{
          cen  = PRIORS$rH0tr(1, lower, upper)
          t = sample(x = PRIORS$Pi[[14]][ cohort, G_pos_1:2^14], prob = Gs[cohort, G_pos_1 : 2^14 ], size = 1)
          
          p_c = PRIORS$H0(upper) - PRIORS$H0(lower)
          p_t.c = sum(Gs[cohort, G_pos_1 : 2^14 ])
          
          rvp_t = sum(Gs[cohort, G_pos_1 : G_pos_2 ])
          rvp_c.t = 1 - PRIORS$H0( STATE$KM[[cohort]]$t_obs[who.to.flip]  ) - 1e-16 # to solve 0
          
          alpha = min(1, p_c * p_t.c / rvp_t / rvp_c.t)
          
        }
        
        new_delta[who.to.flip] = di
        new_times[who.to.flip] = min(t,c)
        
        
        km = survfit(Surv(new_times,  new_delta) ~ 1, 
                     conf.type =  data$ci.type[cohort], 
                     conf.int  =  data$ci[cohort])
        
        
        q = new_times[ c(which.max(km$lower <= .5),
                         which.max(km$surv  <= .5),  
                         which.max(km$upper <= .5))]
        
        
        if(data$statistics[cohort]=="both"){
          
          Sx = ifelse(any(km$upper > .5, na.rm = TRUE) & any(km$upper < .5, na.rm = TRUE), 1, 0) *  # obs up
            ifelse(any(km$lower > .5, na.rm = TRUE) & any(km$lower < .5, na.rm = TRUE), 1, 0) *  # obs lw
            (Qs[cohort,1]  == q[1] & Qs[cohort,2]  == q[2] & Qs[cohort,3]  == q[3])
          
        }else if(data$statistics[cohort]== "lower"){
          
          Sx = ifelse(all(km$upper > .5, na.rm = TRUE),1,0)                                      *  # no  up
            ifelse(any(km$lower > .5, na.rm = TRUE) & any(km$lower < .5, na.rm = TRUE), 1, 0)    *  # obs lw
            (Qs[cohort,1]  == q[1] & Qs[cohort,2]  == q[2])
          
          
        }else if(data$statistics[cohort] == "none"){        
          Sx = as.numeric(Qs[cohort,2]  == q[2 ])
          
        }
        
        
        if(runif(1) < alpha * Sx){
          STATE$KM[[cohort]]$censure[who.to.flip]  = cen
          STATE$KM[[cohort]]$times[who.to.flip]    = t
          STATE$KM[[cohort]]$delta[who.to.flip]    = di
          ar1[cohort] = ar1[cohort] + 1
        }
        
      }else if(delta_step=="move"){
        
        new_delta = STATE$KM[[cohort]]$delta
        new_times = STATE$KM[[cohort]]$t_obs
        
        s.and.r = sample.int(data$n[cohort], size = 2)
        
        s_tc = c( sample(x = PRIORS$Pi[[14]][cohort,1:2^14], prob = Gs[cohort,], size = 1),
                  PRIORS$rH0(1))
        r_tc = c( sample(x = PRIORS$Pi[[14]][cohort,1:2^14], prob = Gs[cohort,], size = 1),
                  PRIORS$rH0(1))
        
        new_times[s.and.r] = c(min(s_tc), min(r_tc))
        
        new_delta[s.and.r] = as.numeric(c(s_tc[1]<s_tc[2],
                                          r_tc[1]<r_tc[2]))
        
        km = survfit(Surv(new_times,  new_delta) ~ 1, 
                     conf.type =  data$ci.type[cohort], 
                     conf.int  =  data$ci[cohort])
        
        
        q = sort(new_times)[ c(which.max(km$lower <= .5),
                               which.max(km$surv  <= .5),  
                               which.max(km$upper <= .5))]
        
        if(data$statistics[cohort]=="both"){
          
          Sx = ifelse(any(km$upper > .5, na.rm = TRUE) & any(km$upper < .5, na.rm = TRUE), 1, 0) *  # obs up
            ifelse(any(km$lower > .5, na.rm = TRUE) & any(km$lower < .5, na.rm = TRUE), 1, 0) *  # obs lw
            (Qs[cohort,1]  == q[1] & Qs[cohort,2]  == q[2] & Qs[cohort,3]  == q[3])
          
        }else if(data$statistics[cohort]== "lower"){
          
          Sx = ifelse(all(km$upper > .5, na.rm = TRUE),1,0)                                      *  # no  up
            ifelse(any(km$lower > .5, na.rm = TRUE) & any(km$lower < .5, na.rm = TRUE), 1, 0) *  # obs lw
            (Qs[cohort,1]  == q[1] & Qs[cohort,2]  == q[2])
          
          
        }else if(data$statistics[cohort] == "none"){        
          Sx = as.numeric(Qs[cohort,2]  == q[2])
          
        }
        
        
        if(Sx == 1){
          ord = order(new_times)
          
          STATE$KM[[cohort]]$t_obs = new_times[ord]
          
          STATE$KM[[cohort]]$times[s.and.r] = c(s_tc[1],r_tc[1])
          STATE$KM[[cohort]]$times          =  STATE$KM[[cohort]]$times[ord]
          
          STATE$KM[[cohort]]$censure[s.and.r] = c(s_tc[2],r_tc[2])
          STATE$KM[[cohort]]$censure          =  STATE$KM[[cohort]]$censure[ord]
          
          STATE$KM[[cohort]]$delta = new_delta[ord]
          
          ar2[cohort] = ar2[cohort] + 1
        }          
      }
      
    }
    
    
    
    STATE$N = samples.size_step(STATE$KM, Qs)
    
    
    # SAMPLER STORE 
    
    if( ((max.iter - iter) %% thinning) == 0 & iter > burn.in ){
      count = count + 1
      
      STATE$medians = sapply(1:n_observed_study, function(study) PRIORS$Pi[[14]][study, 
                                                                                 which.max(cumsum(Gs[study,])>.5)] )
      
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
  
  list("chain" = MCMC_list,
       "time"  = t0 - tend,
       "ar"    = cbind(ar1,ar2)/max.iter)
}

Model.data = readRDS("Model_02112023.rds")

### Plot 1 (observed G_i)           ----

observed_medians = sapply(Model.data$chain, function(x) x$medians)
cbind(t(apply(observed_medians,1, function(x) quantile(x, prob = c(.025,.5,.975)))) 
      #,pfsr_data$lower,  pfsr_data$median, pfsr_data$upper, pfsr_data$n 
)

cbind(t(apply(observed_medians,1, function(x){
  int = HDInterval::hdi(x)
  c(int[1],quantile(x, prob = .5), int[2])})) 
  #,pfsr_data$lower,  pfsr_data$median, pfsr_data$upper, pfsr_data$n 
)


CredI = cbind(t(apply(observed_medians,1,  function(x){
  int = HDInterval::hdi(x)
  c(int[1],quantile(x, prob = .5), int[2])}
)))


pfsr_plot = pfsr_data[,c(1,3,4,5,6,7,13)]

tumor_plot        = metadata$tumor
names(tumor_plot) = metadata$k

agent_plot        = metadata$agent
names(agent_plot) = metadata$k

pfsr_plot$tumor = tumor_plot[pfsr_plot$k]
pfsr_plot$agent = agent_plot[pfsr_plot$k]

pfsr_plot$lower[is.na(pfsr_plot$lower)] = 0
pfsr_plot$upper[is.na(pfsr_plot$upper)] = Inf

pfsr_plot = rbind(pfsr_plot,pfsr_plot) 
pfsr_plot[1:84,2:4] = CredI[,c(2,1,3)]
pfsr_plot$estimate       = "Observed"
pfsr_plot$estimate[1:84] = "Posterior"

index = rep(-1, 84*2)
even = seq(from = 1, to = 167, by = 2)
odd  = seq(from = 2, to = 168, by = 2)

index[even]  = 1:84
index[odd]   = 85:168
pfsr_plot = pfsr_plot[index,]
pfsr_plot$label = as.factor(pfsr_plot$label)
pfsr_plot$marker = factor(pfsr_plot$marker, 
                          levels = c("positive", "negative"),
                          labels = c("positive", "negative"))

pfsr_plot$tumor[pfsr_plot$tumor == "breast"] = "Breast cancer"
pfsr_plot$tumor[pfsr_plot$tumor == "NSCLC"] = "Non-small-cell lung cancer"
pfsr_plot$tumor[pfsr_plot$tumor == "melanoma"] = "Melanoma"
pfsr_plot$tumor[pfsr_plot$tumor == "other"] = "Other"
pfsr_plot$tumor[pfsr_plot$tumor == "Other" & pfsr_plot$label == "Tamura et al. (2019)"] = "Other (C.I. 80%)"

pfsr_plot$tumor = factor(pfsr_plot$tumor, 
                         levels = c("Non-small-cell lung cancer","Breast cancer","Melanoma", "Other", "Other (C.I. 80%)"),
                         labels = c("Non-small-cell lung cancer","Breast cancer","Melanoma", "Other", "Other (C.I. 80%)"))

pfsr_plot$estimate = factor(pfsr_plot$estimate,
                            levels = c("Posterior", "Observed"),
                            labels = c("Posterior", "Observed") )

pfsr_plot$agent = factor(pfsr_plot$agent,
                         levels = c("durvalumab","avelumab","nivolumab","pembrolizumab","atezolizumab","ipilimumab/nivolumab","pembrolizumab-or-nivolumab"),
                         labels = c("durvalumab","avelumab","nivolumab","pembrolizumab","atezolizumab","hybrid treat. 1","hybrid treat. 2"))


plot.1_0 <- pfsr_plot[pfsr_plot$tumor == "Non-small-cell lung cancer",] %>% 
  ggplot( aes(x = median, y = label,
              xmin = lower, xmax = upper, shape = estimate,
              linetype = estimate, color = marker )) +
  geom_pointrange(size = .25, position = position_dodge2(width = 0.8, preserve = "single")) +
  facet_wrap(~tumor, scales = "free_y", ncol = 1) +
  coord_cartesian(xlim = c(0,25)) +
  scale_shape_manual(values=c(3, 4))+
  xlab("time (month)") + ylab(NULL) + theme_bw() + 
  theme(legend.position = "none", text = element_text(family = "serif"))

plot.1_1 <- pfsr_plot[pfsr_plot$tumor != "Other",] %>% 
  ggplot( aes(x = median, y = label,
              xmin = lower, xmax = upper, shape = estimate,
              linetype = estimate, color = marker )) +
  geom_pointrange(size = .25, position = position_dodge2(width = 0.8, preserve = "single")) +
  ggforce::facet_col(vars(tumor),
                     space="free",
                     scale="free_y") +
  coord_cartesian(xlim = c(0,25)) +
  scale_shape_manual(values=c(3, 4))+
  xlab("months") + ylab(NULL) + theme_bw() + 
  theme(legend.position = "none", text = element_text(family = "serif"))




plot.1_2 <- pfsr_plot[pfsr_plot$tumor == "Other" & pfsr_plot$label != "Tamura et al. (2019)",] %>%
  ggplot( aes(x = median, y = label,
              xmin = lower, xmax = upper, shape = estimate,
              linetype = estimate, color = marker )) +
  geom_pointrange(size = .25, position = position_dodge2(width = 1 , preserve = "single")) +
  ggforce::facet_col(vars(tumor),
                     space="free",
                     scale="free_y") +
  coord_cartesian(xlim = c(0,25)) +
  scale_shape_manual(values=c(3, 4))+
  xlab("months") + ylab(NULL) + theme_bw() + 
  theme(legend.position = "none", text = element_text(family = "serif"))

plot.1 = ( plot.1_1 | plot.1_2) + plot_layout(widths = c(1, 1))

ggsave(filename = "plot_data_post_medians.pdf",
       plot = plot.1,
       width  = 21*1.2,
       height = 24*1.2,
       units = "cm",
       dpi = 800)


plot.1_0grey <- pfsr_plot[pfsr_plot$tumor == "Non-small-cell lung cancer",] %>% 
  ggplot( aes(x = median, y = label,
              xmin = lower, xmax = upper, shape = estimate,
              linetype = estimate, color = marker )) +
  geom_pointrange(size = .25, position = position_dodge2(width = 0.8, preserve = "single")) +
  facet_wrap(~tumor, scales = "free_y", ncol = 1) +
  coord_cartesian(xlim = c(0,25)) +
  scale_shape_manual(values=c(3, 4))+
  scale_color_manual(values=c("black","darkgrey")) +
  xlab("time (month)") + ylab(NULL) + theme_bw() + 
  theme(legend.position = "none", text = element_text(family = "serif"))

plot.1_1grey <- pfsr_plot[pfsr_plot$tumor != "Other",] %>% 
  ggplot( aes(x = median, y = label,
              xmin = lower, xmax = upper, shape = estimate,
              linetype = estimate, color = marker )) +
  geom_pointrange(size = .25, position = position_dodge2(width = 0.8, preserve = "single")) +
  ggforce::facet_col(vars(tumor),
                     space="free",
                     scale="free_y") +
  coord_cartesian(xlim = c(0,25)) +
  scale_shape_manual(values=c(3, 4))+
  scale_color_manual(values=c("black","darkgrey")) +
  xlab("months") + ylab(NULL) + theme_bw() + 
  theme(legend.position = "none", text = element_text(family = "serif"))




plot.1_2grey <- pfsr_plot[pfsr_plot$tumor == "Other" & pfsr_plot$label != "Tamura et al. (2019)",] %>%
  ggplot( aes(x = median, y = label,
              xmin = lower, xmax = upper, shape = estimate,
              linetype = estimate, color = marker )) +
  geom_pointrange(size = .25, position = position_dodge2(width = 1 , preserve = "single")) +
  scale_color_manual(values=c("black","darkgrey")) +
  ggforce::facet_col(vars(tumor),
                     space="free",
                     scale="free_y") +
  coord_cartesian(xlim = c(0,25)) +
  scale_shape_manual(values=c(3, 4))+
  xlab("months") + ylab(NULL) + theme_bw() + 
  theme(legend.position = "none", text = element_text(family = "serif"))

plot.1grey = ( plot.1_1grey | plot.1_2grey) + plot_layout(widths = c(1, 1))

ggsave(filename = "plot_data_post_medians_grey.pdf",
       plot = plot.1grey,
       width  = 21*1.2,
       height = 24*1.2,
       units  = "cm",
       dpi    = 800)


### Plot 2 & 3 post processing      ----


new_index = (n_observed_study+1) : (n_observed_study+n_futures_study)
old_index = 1:(n_observed_study)

Z1.0      = PAR[[1]][[1]]$Ex [ new_index ]
S11.0     = PAR[[1]][[1]]$Var[ new_index , new_index ]
S12.0     = PAR[[1]][[1]]$Var[ new_index , old_index ]
S21.0     = PAR[[1]][[1]]$Var[ old_index , new_index ]
S22.0     = PAR[[1]][[1]]$Var[ old_index , old_index ]

Z2.0      = PAR[[1]][[1]]$Ex[old_index]
S22inv.0  = chol2inv( chol( PAR[[1]][[1]]$Var[old_index, old_index]))

# Z00
Z1.00      = PAR[[2]][[1]]$Ex [ new_index ]
S11.00     = PAR[[2]][[1]]$Var[ new_index , new_index ]
S12.00     = PAR[[2]][[1]]$Var[ new_index , old_index ]
S21.00     = PAR[[2]][[1]]$Var[ old_index , new_index ]
S22.00     = PAR[[2]][[1]]$Var[ old_index , old_index ]

Z2.00      = PAR[[2]][[1]]$Ex[old_index]
S22inv.00  = chol2inv( chol( PAR[[2]][[1]]$Var[old_index, old_index]))

# Z10
Z1.10      = PAR[[2]][[2]]$Ex [ new_index ]
S11.10     = PAR[[2]][[2]]$Var[ new_index , new_index ]
S12.10     = PAR[[2]][[2]]$Var[ new_index , old_index ]
S21.10     = PAR[[2]][[2]]$Var[ old_index , new_index ]
S22.10     = PAR[[2]][[2]]$Var[ old_index , old_index ]

Z2.10      = PAR[[2]][[2]]$Ex[old_index]
S22inv.10  = chol2inv( chol( PAR[[2]][[2]]$Var[old_index, old_index]))

Pi_new = Pi.gen( Q = matrix( apply(Qs[ 1: n_observed_study, ],2, median),
                             nrow = 32, ncol = 3, byrow = TRUE) , truncation =  14, G0 = G0,G0inv = G0inv)

Medians = NULL
Mixture.Median = NULL

Ex_Medians = NULL
Ex_Mixture.Median = NULL
weights = weights / sum(weights)


W = matrix( rep(weights, each = 2^14),
            ncol = 2^14, nrow = 16, byrow = TRUE)

## Posterior predictive

set.seed(1414)

obj = readRDS(file = "Post_real_data.rds")

Ex_m_n = obj$Ex_m_n
Ex_m_p = obj$Ex_m_p
Ex_Medians = obj$Ex_Medians
Ex_negative = obj$Ex_negative
Ex_positive = obj$Ex_positive
m_n = obj$m_n
m_p = obj$m_p 
Medians = obj$Medians
Mixture.Median    = obj$Mixture.Median   
Ex_Mixture.Median = obj$Ex_Mixture.Median
### Plot 2 (futrue G_i and P_i)     ----
df.plot.2 = NULL

for(i in 1:32){
  df.append = data.frame( x = Medians[,i],
                          agent   = metadata_future[i,"agent"],
                          tumor   = metadata_future[i,"tumor"],
                          marker  = metadata_future[i,"marker"]) 
  df.plot.2 = rbind( df.plot.2, df.append)
}

df.plot.2$agent = factor(df.plot.2$agent,
                         levels = c("durvalumab","avelumab","nivolumab","pembrolizumab","atezolizumab","ipilimumab/nivolumab","pembrolizumab-or-nivolumab"),
                         labels = c("durvalumab","avelumab","nivolumab","pembrolizumab","atezolizumab","hybrid treatment (1)","hybrid treatment (2)"))


df.plot.3 = NULL
for(i in 1:16){
  
  df.append = data.frame( x = Medians[, i*2-1] - Medians[, i*2],
                          agent   = metadata_future[i*2,"agent"],
                          tumor   = metadata_future[i*2,"tumor"])
  
  
  df.plot.3 = rbind( df.plot.3, df.append)
}

for(ag in unique(metadata_future$agent)){
  for(tt in unique(metadata_future$tumor)){
    index = df.plot.3$tumor == tt & df.plot.3$agent == ag
    df.plot.3$lower[index] =  hdi(df.plot.3$x[index], credMass = 0.95)[1]
    df.plot.3$upper[index] =  hdi(df.plot.3$x[index], credMass = 0.95)[2]
  }
}


future_study_plot = ggplot(df.plot.2, aes(x = x, y = tumor , fill = factor(marker, levels = c("positive","negative")))) + 
  geom_density_ridges(scale = 1.2, alpha = .5, bandwidth = .3) + 
  facet_grid( tumor ~ agent, scales = "free_y") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_fill_discrete("") +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10)) + 
  theme(legend.position = "none", text = element_text(family = "serif"))

future_study_plot_diff = ggplot(df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 2, bandwidth = 0.3) +
  facet_grid( tumor ~ agent, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "red") +
  scale_fill_manual(values = c("transparent", "green", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,4)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

future_study_plot_diff_t = ggplot(df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 2, bandwidth = 0.3) +
  facet_grid( agent ~ tumor, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "red") +
  scale_fill_manual(values = c("transparent", "green", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,6), expand = c(0,0)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

ggsave2("future_study_plot.pdf",        future_study_plot        , units = c("cm"), width = 21, height = 15)
ggsave2("future_study_plot_diff.pdf",   future_study_plot_diff   , units = c("cm"), width = 21, height = 15)
ggsave2("future_study_plot_diff_t.pdf", future_study_plot_diff_t , units = c("cm"), width = 21, height = 15)

data.plot.overall = data.frame(Mixture.Median[,1:2])
colnames(data.plot.overall) = c("positive","negative")
data.plot.overall = melt(data.plot.overall)
data.plot.overall$title = "Target population (Mixture)"
plot.mix = ggplot(data = data.plot.overall, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.5, bw = .3) +
  facet_wrap(~title) +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits=c(2,6), expand = c(0,0)) +
  scale_fill_discrete(NULL)+
  theme(legend.position = "top", text = element_text(family = "serif"))
plot.mix     

ggsave2("Over_ALL.pdf", plot.mix , units = c("cm"), width = 10, height = 10)

### Ex plots

Ex_df.plot.2 = NULL

for(i in 1:32){
  df.append = data.frame( x = Ex_Medians[,i],
                          agent   = metadata_future[i,"agent"],
                          tumor   = metadata_future[i,"tumor"],
                          marker  = metadata_future[i,"marker"]) 
  Ex_df.plot.2 = rbind( Ex_df.plot.2, df.append)
}

Ex_df.plot.2$agent = factor(Ex_df.plot.2$agent,
                            levels = c("durvalumab","avelumab","nivolumab","pembrolizumab","atezolizumab","ipilimumab/nivolumab","pembrolizumab-or-nivolumab"),
                            labels = c("durvalumab","avelumab","nivolumab","pembrolizumab","atezolizumab","hybrid treatment (1)","hybrid treatment (2)"))



Ex_future_study_plot = ggplot(Ex_df.plot.2, aes(x = x, y = tumor , fill = factor(marker, levels = c("positive","negative")))) + 
  geom_density_ridges(scale = 1.2, alpha = .5, bandwidth = 0.3) + 
  facet_grid( tumor ~ agent, scales = "free_y") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_fill_discrete("") + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  theme(legend.position = "none", text = element_text(family = "serif"))


Ex_df.plot.3 = NULL
for(i in 1:16){
  
  df.append = data.frame( x = Ex_Medians[, i*2-1] - Ex_Medians[, i*2],
                          agent   = metadata_future[i*2,"agent"],
                          tumor   = metadata_future[i*2,"tumor"])
  
  
  Ex_df.plot.3 = rbind( Ex_df.plot.3, df.append)
}

for(ag in unique(metadata_future$agent)){
  for(tt in unique(metadata_future$tumor)){
    index = Ex_df.plot.3$tumor == tt & Ex_df.plot.3$agent == ag
    Ex_df.plot.3$lower[index] =  hdi(Ex_df.plot.3$x[index])[1]
    Ex_df.plot.3$upper[index] =  hdi(Ex_df.plot.3$x[index])[2]
  }
}


# 2500*32*2

Ex_future_study_plot = ggplot(Ex_df.plot.2, aes(x = x, y = tumor , fill = factor(marker, levels = c("positive","negative")))) + 
  geom_density_ridges(scale = 1.2, alpha = .5, bandwidth = 0.3) + 
  facet_grid( tumor ~ agent, scales = "free_y") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_fill_discrete("") +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10)) + 
  theme(legend.position = "none", text = element_text(family = "serif"))



Ex_future_study_plot_diff = ggplot(Ex_df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  #geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 2, bandwidth = 0.3) +
  facet_grid( tumor ~ agent, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "red") +
  scale_fill_manual(values = c("transparent", "green", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,4)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

Ex_future_study_plot_diff_t = ggplot(Ex_df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 2, bandwidth = 0.3) +
  facet_grid( agent ~ tumor, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "red") +
  scale_fill_manual(values = c("transparent", "green", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,4)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

ggsave2("Ex_future_study_plot.pdf",        Ex_future_study_plot        , units = c("cm"), width = 21, height = 15,
        dpi = 800)
ggsave2("Ex_future_study_plot_diff.pdf",   Ex_future_study_plot_diff   , units = c("cm"), width = 21, height = 15,
        dpi = 800)
ggsave2("Ex_future_study_plot_diff_t.pdf", Ex_future_study_plot_diff_t , units = c("cm"), width = 21, height = 15,
        dpi = 800)

Ex_data.plot.overall = data.frame(Ex_Mixture.Median[,1:2])
colnames(Ex_data.plot.overall) = c("positive","negative")
Ex_data.plot.overall = melt(Ex_data.plot.overall)
Ex_data.plot.overall$title = "Target population (Mixture)"
Ex_plot.mix = ggplot(data = Ex_data.plot.overall, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.5, bw = .3) +
  facet_wrap(~title) +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits=c(2,6), expand = c(0,0)) +
  scale_fill_discrete(NULL)+
  theme(legend.position = "top", text = element_text(family = "serif"))
Ex_plot.mix     

ggsave2("Ex_Over_ALL.pdf", Ex_plot.mix , units = c("cm"), width = 10, height = 10)


#### GREY SCALE PLOTS ----


future_study_plot_grey = ggplot(df.plot.2, aes(x = x, y = tumor , fill = factor(marker, levels = c("positive","negative")))) + 
  geom_density_ridges(scale = 1.2, alpha = .75, bandwidth = .3) + 
  facet_grid( tumor ~ agent, scales = "free_y") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=c("black","darkgrey")) +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10)) + 
  theme(legend.position = "none", text = element_text(family = "serif"))

future_study_plot_diff_grey = ggplot(df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 0, bandwidth = 0.3) +
  facet_grid( tumor ~ agent, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "black") +
  scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,4), expand = c(0,0)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

future_study_plot_diff_t_grey = ggplot(df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 0, bandwidth = 0.3) +
  facet_grid( agent ~ tumor, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "black") +
  scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,6), expand = c(0,0)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

ggsave2("future_study_plot_grey.pdf", future_study_plot_grey,
        units = c("cm"), width = 21, height = 15,
        dpi = 800)
ggsave2("future_study_plot_diff_grey.pdf", future_study_plot_diff_grey,
        units = c("cm"), width = 21, height = 15,
        dpi = 800)

ggsave2("future_study_plot_diff_t_grey.pdf", future_study_plot_diff_t_grey,
        units = c("cm"), width = 21, height = 15,
        dpi = 800)

plot.mix_grey = ggplot(data = data.plot.overall, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.75, bw = .3) +
  facet_wrap(~title) +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits=c(2,6), expand = c(0,0)) +
  scale_fill_manual(NULL, values=c("black","darkgrey")) +
  theme(legend.position = "top", text = element_text(family = "serif"))


ggsave2("Over_ALL_grey.pdf", plot.mix_grey , units = c("cm"), width = 10, height = 10,
        dpi = 800)

### Ex plots grey

Ex_future_study_plot_grey = ggplot(Ex_df.plot.2, aes(x = x, y = tumor , fill = factor(marker, levels = c("positive","negative")))) + 
  geom_density_ridges(scale = 1.2, alpha = .5, bandwidth = 0.3) + 
  facet_grid( tumor ~ agent, scales = "free_y") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_fill_manual("", values=c("black","darkgrey")) +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10)) + 
  theme(legend.position = "none", text = element_text(family = "serif"))


Ex_future_study_plot_diff_grey = ggplot(Ex_df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 0, bandwidth = 0.3) +
  facet_grid( tumor ~ agent, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "black") +
  scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,4)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

Ex_future_study_plot_diff_t_grey = ggplot(Ex_df.plot.3, aes(x = x, y = NA, fill = stat(quantile))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE, 
                               quantiles = c(0.05, 1),
                               vline_linetype = 0, bandwidth = 0.3) +
  facet_grid( agent ~ tumor, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, col = "black") +
  scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none") +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits = c(-1,4)) + 
  theme(legend.position = "top", text = element_text(family = "serif"))

ggsave2( "Ex_future_study_plot_grey.pdf",
         Ex_future_study_plot_grey, units = c("cm"), width = 21, height = 15,
         dpi = 800)
ggsave2( "Ex_future_study_plot_diff_grey.pdf",   Ex_future_study_plot_diff_grey,
         units = c("cm"), width = 21, height = 15,
         dpi = 800)
ggsave2( "Ex_future_study_plot_diff_t_grey.pdf", Ex_future_study_plot_diff_t_grey,
         units = c("cm"), width = 21, height = 15,
         dpi = 800)

Ex_plot.mix_grey = ggplot(data = Ex_data.plot.overall, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.5, bw = .3) +
  facet_wrap(~title) +
  theme_bw() + xlab(NULL) + ylab(NULL) + 
  scale_y_discrete(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(limits=c(2,6), expand = c(0,0)) +
  scale_fill_manual(NULL, values=c("black","darkgrey")) +
  theme(legend.position = "top", text = element_text(family = "serif"))
Ex_plot.mix_grey     

ggsave2("Ex_Over_ALL_grey.pdf", Ex_plot.mix_grey , units = c("cm"), width = 10, height = 10,
        dpi = 800)
