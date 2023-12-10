library(metafor)
library(survival)
library(mvtnorm)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)
library(latex2exp)

set.seed(1)

file_list = list.files("PTGP_simulation/")
dir       = "PTGP_simulation/"

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


obj = readRDS("RESULTS_GPPT_SIM_v2.rds")


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