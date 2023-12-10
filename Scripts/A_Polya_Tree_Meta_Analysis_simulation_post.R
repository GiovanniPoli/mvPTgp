library(metafor)
library(survival)
library(mvtnorm)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)

set.seed(1)

file_list = list.files("~/PTGP_simulation/")[351:400]
dir       = "~/PTGP_simulation/"

BF_from_pvalues = function(p){
  ifelse(p < 1 /exp(1) , - 1 / (exp(1) * p * log(p)), 1 )}
gg_color_hue    = function(n){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

t0 = Sys.time()

p.values  = matrix(NA, ncol = 7, nrow = length(file_list))
post.prob = matrix(NA, ncol = 7, nrow = length(file_list))
post.prob_Ex = matrix(NA, ncol = 7, nrow = length(file_list))

freq.Est      = matrix(NA, ncol = 7, nrow = length(file_list))
Bayes.ExlogHR = matrix(NA, ncol = 7, nrow = length(file_list))
Bayes.logHR   = matrix(NA, ncol = 7, nrow = length(file_list))

erfinv = function(x) qnorm((1+x)/2)/sqrt(2)
expit  = function(z) 1/(1+exp(-z))
logit  = function(p) if(all(p>0)){log(p)-log(1-p)}else{ log(p + 1e-16) - log( 1 - p - 1e-16) }

g0     = function(x, scale = 3.5){exp(log(2 * scale) - log(pi * { x * x + scale * scale }))}
G0     = function(q, scale = 3.5){ {2/pi} * atan(q/scale) }
G0inv  = function(p, scale = 3.5){ scale * tan({ pi * p }/2) }

count = 0

keep = c(list = ls(),"keep")

### Correction of mistakes and summary 

for(file in file_list){
  
  count = count + 1
  
  SIM   = readRDS(paste0(dir,file))

  metareg = rma(yi     = SIM$data.meta$logHR,
                vi     = SIM$data.meta$sd.HR^2,
                mods   = ~0 + tumor * agent, data = SIM$data.meta,
                method = "DL")
  meta.avg = rma(yi     = SIM$data.meta$logHR,
                 vi     = SIM$data.meta$sd.HR^2,
                 data   = SIM$data.meta,
                 method = "DL")

  SIM$metareg  = metareg
  SIM$meta.avg = meta.avg
  
  z    = c(metareg$beta["tumortumor 1",],
           metareg$beta["tumortumor 1",] + metareg$beta["agentagent 2",],
           metareg$beta["tumortumor 2",],
           metareg$beta["tumortumor 2",] + metareg$beta["agentagent 2",] + metareg$beta["tumortumor 2:agentagent 2",],
           metareg$beta["tumortumor 3",],
           metareg$beta["tumortumor 3",] + metareg$beta["agentagent 2",] + metareg$beta["tumortumor 3:agentagent 2",])
  
  S = metareg$vb
  
  C = matrix(c(1,0,0,0,0,0,
               1,0,0,1,0,0,
               0,1,0,0,0,0,
               0,1,0,1,1,0,
               0,0,1,0,0,0,
               0,0,1,1,0,1), ncol = 6, nrow = 6, byrow = TRUE)
  
  sd = sqrt(diag(C%*% S %*% t(C)))
  
  p.value = pnorm( z, mean = 0, sd = sd, lower.tail = FALSE) 
  
  ### AVG -> Effect
  
  # Z0
  Z1.0      = SIM$PRIORS$G_star[[1]][[1]]$Ex[51:62]
  S11.0     = SIM$PRIORS$G_star[[1]][[1]]$Var[51:62, 51:62]
  S12.0     = SIM$PRIORS$G_star[[1]][[1]]$Var[51:62,  1:50]
  S21.0     = SIM$PRIORS$G_star[[1]][[1]]$Var[ 1:50, 51:62]
  
  Z2.0      = SIM$PRIORS$GP[[1]][[1]]$Ex
  S22inv.0  = SIM$PRIORS$GP[[1]][[1]]$Prec
  
  # Z00
  Z1.00      = SIM$PRIORS$G_star[[2]][[1]]$Ex[51:62]
  S11.00     = SIM$PRIORS$G_star[[2]][[1]]$Var[51:62, 51:62]
  S12.00     = SIM$PRIORS$G_star[[2]][[1]]$Var[51:62,  1:50]
  S21.00     = SIM$PRIORS$G_star[[2]][[1]]$Var[ 1:50, 51:62]
  Z2.00      = SIM$PRIORS$GP[[2]][[1]]$Ex
  S22inv.00  = SIM$PRIORS$GP[[2]][[1]]$Prec
  
  # Z10
  Z1.10      = SIM$PRIORS$G_star[[2]][[2]]$Ex[51:62]
  S11.10     = SIM$PRIORS$G_star[[2]][[2]]$Var[51:62, 51:62]
  S12.10     = SIM$PRIORS$G_star[[2]][[2]]$Var[51:62,  1:50]
  S21.10     = SIM$PRIORS$G_star[[2]][[2]]$Var[ 1:50, 51:62]
  Z2.10      = SIM$PRIORS$GP[[2]][[2]]$Ex
  S22inv.10  = SIM$PRIORS$GP[[2]][[2]]$Prec
  
  
  for(obs in 1:length(SIM$chain)){
    
    cat("\r[",obs,"/ 2500 ] - Sim:", count)
    
    Zobs.0  = logit(SIM$chain[[obs]]$Ye0[[1]][,1])
    Zobs.00 = logit(SIM$chain[[obs]]$Ye0[[2]][,1])
    Zobs.10 = logit(SIM$chain[[obs]]$Ye0[[2]][,2])
    
    Ex_Z0.new      = Z1.0  + S12.0  %*% S22inv.0  %*% (Zobs.0 - Z2.0) 
    Ex_Z00.new     = Z1.00 + S12.00 %*% S22inv.00 %*% (Zobs.00 - Z2.00) 
    Ex_Z10.new     = Z1.10 + S12.10 %*% S22inv.10 %*% (Zobs.10 - Z2.10) 
    
    Var_Z0.new     = SIM$PRIORS$G_star[[1]][[1]]$Var[51:62, 51:62] - S12.0  %*% S22inv.0  %*% t(S12.0)
    Var_Z00.new    = SIM$PRIORS$G_star[[2]][[1]]$Var[51:62, 51:62] - S12.00 %*% S22inv.00 %*% t(S12.00)
    Var_Z10.new    = SIM$PRIORS$G_star[[2]][[2]]$Var[51:62, 51:62] - S12.10 %*% S22inv.10 %*% t(S12.10) 
    
    Ex_Yi0   = colMeans(expit(rmvnorm(10000, mean = c(Ex_Z0.new),  sigma = Var_Z0.new)))
    Ex_Yi00  = colMeans(expit(rmvnorm(10000, mean = c(Ex_Z00.new), sigma = Var_Z00.new)))
    Ex_Yi10  = colMeans(expit(rmvnorm(10000, mean = c(Ex_Z10.new), sigma = Var_Z10.new)))
    
    Ex_Yie = t(sapply(Ex_Yi0, function(p) rep(c(p,1-p), each = 2^(14-1)))) 
    
    Ex_Yiee = cbind( t(sapply(Ex_Yi00, function(p) rep(c(p,1-p), each = 2^(14-2)))) ,
                     t(sapply(Ex_Yi10, function(p) rep(c(p,1-p), each = 2^(14-2)))))
    
    
    Ex_Gs = 1/2^(14-6) * Ex_Yie * Ex_Yiee
    
    Yi0   = c(expit(rmvnorm( 1, mean = c(Ex_Z0.new),  sigma = Var_Z0.new)))
    Yi00  = c(expit(rmvnorm(1,  mean = c(Ex_Z00.new), sigma = Var_Z00.new)))
    Yi10  = c(expit(rmvnorm(1,  mean = c(Ex_Z10.new), sigma = Var_Z10.new)))
    
    Yie = t(sapply( Yi0, function(p) rep(c(p,1-p), each = 2^(14-1)))) 
    
    Yiee = cbind( t(sapply( Yi00, function(p) rep(c(p,1-p), each = 2^(14-2)))) ,
                   t(sapply( Yi10, function(p) rep(c(p,1-p), each = 2^(14-2)))))
    
    
    Gs = 1/2^(14-6) * Yie * Yiee

    for(deep in 3:6){
      Ex_G_deep = NULL
      G_deep    = NULL
      for(combination in 1:2^(deep-1)){
        Ex_Ye0 = SIM$PRIORS$G_star[[deep]][[combination]]$Y_Ex[51:62]
        Ex_Ye1 = 1 - Ex_Ye0
        
        Ye0 = rmvnorm(1,  mean  = SIM$PRIORS$G_star[[deep]][[combination]]$Ex[51:62],
                          sigma = SIM$PRIORS$G_star[[deep]][[combination]]$Var[51:62, 51:62])
        Ye1 = 1 - Ye0
        
        Ex_G_deep = cbind(Ex_G_deep, matrix(rep(Ex_Ye0,2^(14-deep)), ncol = 2^(14-deep), nrow = 12),
                       matrix(rep(Ex_Ye1,2^(14-deep)), ncol = 2^(14-deep), nrow = 12))
        G_deep    = cbind(G_deep,    matrix(rep(Ye0,2^(14-deep)),    ncol = 2^(14-deep), nrow = 12),
                                     matrix(rep(Ye1,2^(14-deep)),    ncol = 2^(14-deep), nrow = 12))
      }
      Ex_Gs = Ex_G_deep * Ex_Gs
      Gs =    G_deep *    Gs
    }
    
    W = matrix ( rep(c(5/25, 5/25,
                       5/25, 5/25,
                       3/25, 2/25), each = 2^14 ) , ncol = 2^14, nrow = 6, byrow = TRUE)    
    
    
    
    Ex_positive = Ex_Gs[c(1,3,5,7, 9,11),] * W
    Ex_negative = Ex_Gs[c(2,4,6,8,10,12),] * W
    
    positive = Gs[c(1,3,5,7, 9,11),] * W
    negative = Gs[c(2,4,6,8,10,12),] * W
    
    # Checked: SIM$PRIORS$Pi[[14]][52,2] == SIM$PRIORS$Pi[[14]][51,2]

    Ex_m_p = G0inv(SIM$PRIORS$Pi[[14]][51, which.max(cumsum(colSums(Ex_positive)) > .5)]) # 51-62 Pi is Eq.
    Ex_m_l = G0inv(SIM$PRIORS$Pi[[14]][51, which.max(cumsum(colSums(Ex_negative)) > .5)]) #
    
    m_p = G0inv( SIM$PRIORS$Pi[[14]][51, which.max(cumsum(colSums( positive)) > .5)])   # 
    m_l = G0inv( SIM$PRIORS$Pi[[14]][51, which.max(cumsum(colSums( negative)) > .5)])   #
    
    SIM$chain[[obs]]$Ex_new_study_medians = c(apply(Ex_Gs,1,function(G) G0inv(SIM$PRIORS$Pi[[14]][51,which.max(cumsum(G)>.5)])))
    SIM$chain[[obs]]$Ex_Overall_medians   = c("pos." = Ex_m_p, "neg." = Ex_m_l)
    
    SIM$chain[[obs]]$new_study_medians = c(apply(Gs  ,1,function(G) G0inv(SIM$PRIORS$Pi[[14]][51,which.max(cumsum(G)>.5)])))
    SIM$chain[[obs]]$Overall_medians   = c(   "pos." = m_p,    "neg." = m_l)
    
  }

  
  ####
  # Ex_new_study_medians = median of Ex of G_n+1 G_I
  # Ex_Overall_medians   = median of Ex \sum_r=1^I \pi * G_n+r 
  # 
  # new_study_medians = median of a sample from                   G_n+i \mid data 
  # Overall_medians   = median of a sample from  \sum_r=1^I \pi * G_n+r \mid data 
  
  

  SIM$Ex_post_probs    = apply( sapply( SIM$chain, function(x) x$Ex_new_study_medians[ c(1, 3, 5, 7, 9 ,11)] <
                                                               x$Ex_new_study_medians[ c(2, 4, 6, 8, 10,12)]), 1, mean)
 
  SIM$Ex_post_prob_avg     = mean(  sapply( SIM$chain, function(x) x$Ex_Overall_medians[1] < 
                                                                   x$Ex_Overall_medians[2]))
  
  SIM$post_probs    = apply( sapply( SIM$chain, function(x) x$new_study_medians[ c(1, 3, 5, 7, 9 ,11)] <
                                                            x$new_study_medians[ c(2, 4, 6, 8, 10,12)]), 1, mean)
  
  SIM$post_prob_avg     = mean(  sapply( SIM$chain, function(x) x$Overall_medians[1] < 
                                                                x$Overall_medians[2]))
  
  SIM$Ex_logHRs      = apply( sapply( SIM$chain, function(x) log(x$Ex_new_study_medians[ c(1, 3, 5, 7, 9 , 11)] /
                                                                 x$Ex_new_study_medians[ c(2, 4, 6, 8, 10, 12)])), 1, median)
  
  SIM$Ex_logHR_avg   = median(  sapply( SIM$chain, function(x) log(x$Ex_Overall_medians[1] / x$Ex_Overall_medians[2]) ) )
  
  SIM$logHRs         = apply( sapply( SIM$chain, function(x) log(x$new_study_medians[ c(1, 3, 5, 7, 9 ,11)] / x$new_study_medians[ c(2, 4, 6, 8, 10,12)]) ), 1, median)
  
  SIM$logHR_avg      = median(  sapply( SIM$chain, function(x) log(x$Overall_medians[1] / x$Overall_medians[2])))
   

  post.prob_Ex[count,] = c(SIM$Ex_post_probs, SIM$Ex_post_prob_avg)
  post.prob[count,] = c(SIM$post_probs,       SIM$post_prob_avg)
  
  
  p.values[count, ] = c(p.value,"AVG" = pnorm(q = c(meta.avg$beta), sd = meta.avg$se, lower.tail = FALSE))
  
  freq.Est[count,]      = c(SIM$metareg$beta, SIM$meta.avg$beta)
  Bayes.ExlogHR[count,] = c(SIM$Ex_logHRs, SIM$Ex_logHR_avg)
  Bayes.logHR[count,]   = c(SIM$logHRs,    SIM$logHR_avg)
  
  
  saveRDS(object = SIM,file = paste0(dir,"post_proc_",file))
  
  SIM = NULL
  remove(list = ls()[! ls()%in% keep])
}



saveRDS(object = list("Gi"            = post.prob,
                      "Ex_Gi"         = post.prob_Ex,
                      "logHR"         = Bayes.ExlogHR,
                      "Ex_logHR"      = Bayes.logHR,
                      "p.values"      = p.values,
                      "freq.estimate" = freq.Est),
        file = "RESULTS_GPPT_SIM_v2.rds")



post.prob      =  obj$Gi    
post.prob_Ex   =  obj$Ex_Gi
Bayes.ExlogHR  =  obj$logHR
Bayes.logHR    =  obj$Ex_logHR
p.values       =  obj$p.values
freq.Est       =  obj$freq.estimate
 


names = c("(TT1, A0)",
          "(TT1, A1)",
          "(TT2, A0)",
          "(TT2, A1)",
          "(TT3, A0)",
          "(TT3, A1)",
          "Overall")

colnames(post.prob)    = names
colnames(post.prob_Ex) = names
colnames(p.values)     = names


df_pvalues  = melt(as.data.frame(p.values))
df_postprob_Ex = melt(as.data.frame(post.prob_Ex))
df_postprob = melt(as.data.frame(post.prob))


df_BF_Ex       = rbind( melt(as.data.frame( (1-post.prob_Ex) / post.prob_Ex)),
                     melt(as.data.frame( BF_from_pvalues(p.values))))
df_BF       = rbind( melt(as.data.frame( (1-post.prob) / post.prob)),
                      melt(as.data.frame( BF_from_pvalues(p.values))))

df_BF_Ex$variable = df_BF$variable = paste0(df_BF$variable, rep(c(" (Bayesian)"," (Converted)"), each = 350))
df_BF_Ex$type     = df_BF$type     = rep(c("Bayesian","Converted"), each = 350)

plot_11 = ggplot(data = df_pvalues, aes( y = variable, x = value, fill = variable)) +
          geom_boxplot( col = "#808080") + xlab(NULL) + ylab(NULL) + ggtitle("p-values") +
          theme_bw() +
          scale_x_continuous(limits = c(0,1))+
          scale_fill_viridis_d() +
          theme(text = element_text(family = "serif"), legend.position = "none")
plot_11


plot_11_v = ggplot(data = df_pvalues, aes( y = variable, x = value, fill = variable)) +
  geom_violin( col = "#808080") + xlab(NULL) + ylab(NULL) + ggtitle("p-values") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_viridis_d() +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_11_v

plot_12 = ggplot(data = df_postprob, aes( y = variable, x = value, fill = variable)) +
          geom_boxplot( col = "#000000") + xlab(NULL) + ylab(NULL) + ggtitle("posterior probability") +
          theme_bw() +
          scale_x_continuous(limits = c(0,1))+
          scale_fill_viridis_d() +
          theme(text = element_text(family = "serif"), legend.position = "none")
plot_12

plot_12_v = ggplot(data = df_postprob, aes( y = variable, x = value, fill = variable)) +
  geom_violin( col = "#000000") + xlab(NULL) + ylab(NULL) + ggtitle("posterior probability") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_viridis_d() +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_12_v

plot_12_Ex = ggplot(data = df_postprob_Ex, aes( y = variable, x = value, fill = variable)) +
  geom_boxplot( col = "#000000") + xlab(NULL) + ylab(NULL) + ggtitle("posterior probability") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_viridis_d() +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_12_Ex


plot_12_Ex_v = ggplot(data = df_postprob_Ex, aes( y = variable, x = value, fill = variable)) +
  geom_violin( col = "#000000") + xlab(NULL) + ylab(NULL) + ggtitle("posterior probability") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_viridis_d() +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_12_Ex_v

col = scales::viridis_pal()(7)
col = col[rep(1:7,each = 2)]

plot_2 = ggplot(data = df_BF, aes( y = variable, x = value, fill = variable, col = type)) +
         geom_boxplot() + xlab(NULL) + ylab(NULL) + ggtitle("BF (mvPT) vs BF bounds (meta-regression)") +
         theme_bw() +
         scale_x_log10() +
         scale_colour_manual( values = c("#000000", "#808080")) +
         scale_fill_manual(values = col)+
         theme(text = element_text(family = "serif"), legend.position = "none")
plot_2


plot_2_Ex = ggplot(data = df_BF_Ex, aes( y = variable, x = value, fill = variable, col = type)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) + ggtitle("BF (mvPT) vs BF bounds (meta-regression)") +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual( values = c("#000000", "#808080")) +
  scale_fill_manual(values = col)+
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_2_Ex

plot    = grid.arrange( arrangeGrob(plot_12, plot_11,ncol = 2), plot_2, nrow = 2)
plot_Ex = grid.arrange( arrangeGrob(plot_12_Ex, plot_11,ncol = 2), plot_2_Ex, nrow = 2)
plot_pairs = grid.arrange(plot_12_Ex + ggtitle(TeX("posterior probability - (expected values)")),
                          plot_12   + ggtitle(TeX("posterior probability - (future cohorts)")), 
                          ncol = 1)


ggsave("Sim_results_post_probs.pdf",
       plot = plot,  dpi = 800,
       width = 20, height = 20, units = "cm")

ggsave("Sim_results_post_probs_Ex.pdf",
       plot = plot_Ex, dpi = 800,
       width = 20, height = 20, units = "cm")

ggsave("Sim_Bar_vs_sampled.pdf",
       plot = plot_pairs, dpi = 800,
       width = 7, height = 15, units = "cm")


col_bw = hcl.colors(7, "Gray")
col_bw2 = col_bw[rep(1:7,each = 2)]


plot_11_bw = ggplot(data = df_pvalues, aes( y = variable, x = value, fill = variable)) +
  geom_boxplot( lty = 5) + xlab(NULL) + ylab(NULL) + ggtitle("p-values") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_manual(values = col_bw) +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_11_bw


plot_12_bw = ggplot(data = df_postprob, aes( y = variable, x = value, fill = variable)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) + ggtitle("posterior probability") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_manual(values = col_bw) +
  theme(text = element_text(family = "serif"), legend.position = "none")

plot_12_bw

plot_12_Ex_bw = ggplot(data = df_postprob_Ex, aes( y = variable, x = value, fill = variable)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) + ggtitle("posterior probability") +
  theme_bw() +
  scale_x_continuous(limits = c(0,1))+
  scale_fill_manual(values = col_bw) +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_12_Ex_bw



plot_2_bw = ggplot(data = df_BF, aes( y = variable, x = value, fill = variable, linetype = type)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) + ggtitle("BF (mvPT) vs BF bounds (meta-regression)") +
  theme_bw() +
  scale_x_log10() +
  scale_linetype_manual( values = c(1,5)) +
  scale_fill_manual(values = col_bw2) +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_2_bw


plot_2_Ex_bw = ggplot(data = df_BF_Ex, aes( y = variable, x = value, fill = variable, linetype = type)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) + ggtitle("BF (mvPT) vs BF bounds (meta-regression)") +
  theme_bw() +
  scale_x_log10() +
  scale_linetype_manual( values = c(1,5)) +
  scale_fill_manual(values = col_bw2) +
  theme(text = element_text(family = "serif"), legend.position = "none")
plot_2_Ex_bw

plot_bw    = grid.arrange( arrangeGrob(plot_12_bw, plot_11_bw,ncol = 2), plot_2_bw, nrow = 2)
plot_Ex_bw = grid.arrange( arrangeGrob(plot_12_Ex_bw, plot_11_bw,ncol = 2), plot_2_Ex_bw, nrow = 2)
plot_pairs_bw = grid.arrange(plot_12_Ex_bw, plot_12_bw, ncol = 1)


ggsave("Sim_results_post_probs_bw.pdf",
       plot = plot_bw  , dpi = 800,
       width = 20, height = 20, units = "cm")

ggsave("Sim_results_post_probs_Ex_bw.pdf",
       plot = plot_Ex_bw, dpi = 800,
       width = 20, height = 20, units = "cm")

ggsave("Sim_Bar_vs_sampled_bw.pdf",
       plot = plot_pairs_bw, dpi = 800,
       width = 7, height = 15, units = "cm")