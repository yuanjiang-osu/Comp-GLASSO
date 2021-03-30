ROC_plot <- function()
{
  library(cowplot)
  n.rep <- 100
  num_cpu <- 10
  length_rholist <- 70

  # ROC curve of TP/FP from all rhos, in one simulation

  TP.1_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  TP.2_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  TP.3_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  TP.4_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  FP.1_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  FP.2_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  FP.3_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))
  FP.4_array <- array(0, c(num_cpu, n.rep / num_cpu, length_rholist))

for (i in 1:num_cpu) {
    for (j in 1:(n.rep / num_cpu)) {
      for (k in 1:length_rholist) {
        TP.1 <- results[[i]][[j]]$TP.1_s[k]
        TP.2 <- results[[i]][[j]]$TP.2_s[k]
        TP.3 <- results[[i]][[j]]$TP.3_s[k]
        TP.4 <- results[[i]][[j]]$TP.4_s[k]
        FP.1 <- results[[i]][[j]]$FP.1_s[k]
        FP.2 <- results[[i]][[j]]$FP.2_s[k]
        FP.3 <- results[[i]][[j]]$FP.3_s[k]
        FP.4 <- results[[i]][[j]]$FP.4_s[k]
        
        TP.1_array[i, j, k] <- TP.1
        TP.2_array[i, j, k] <- TP.2
        TP.3_array[i, j, k] <- TP.3
        TP.4_array[i, j, k] <- TP.4
        FP.1_array[i, j, k] <- FP.1
        FP.2_array[i, j, k] <- FP.2
        FP.3_array[i, j, k] <- FP.3
        FP.4_array[i, j, k] <- FP.4
      }
    }
}

avg_TP.1 <- apply(TP.1_array, 3, mean, na.rm = TRUE) 
avg_TP.2 <- apply(TP.2_array, 3, mean, na.rm = TRUE) 
avg_TP.3 <- apply(TP.3_array, 3, mean, na.rm = TRUE) 
avg_TP.4 <- apply(TP.4_array, 3, mean, na.rm = TRUE) 
avg_FP.1 <- apply(FP.1_array, 3, mean, na.rm = TRUE) 
avg_FP.2 <- apply(FP.2_array, 3, mean, na.rm = TRUE)
avg_FP.3 <- apply(FP.3_array, 3, mean, na.rm = TRUE) 
avg_FP.4 <- apply(FP.4_array, 3, mean, na.rm = TRUE)

roc <- data.frame(TPR = c(avg_TP.1, avg_TP.2, avg_TP.3), FPR = c(avg_FP.1, avg_FP.2, avg_FP.3), 
                  method = c(rep("Compo-glasso", length_rholist), rep("glasso", length_rholist), rep("MB", length_rholist)))
roc_curve <- ggplot(data = roc, aes(x = FPR, y = TPR, group = method)) +
  geom_path(aes(colour = method, linetype = method)) +
  scale_color_manual(name = "Method", breaks = c("Compo-glasso", "glasso", "MB"), values =c("blue", "red", "black")) +
  scale_linetype_manual(name = "Method", breaks= c("Compo-glasso", "glasso", "MB"), values=c("solid", "longdash", "dotted")) + 
  xlim(0, 1) + 
  ylim(0, 1.05)
# saveRDS(roc_curve, "plot1")
return(roc_curve)
gc()
}

library(tidyverse)
# Type 1
results = read_rds("T1_n_100_K_200_l_l_MB.rds")
T1_K_200_l_info_l_var <- ROC_plot() + 
  labs(caption = "h, h") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("T1_n_100_K_200_l_m_MB.rds")
T1_K_200_l_info_m_var <- ROC_plot() +
  labs(caption = "h, l") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("T1_n_100_K_200_m_l_MB.rds")
T1_K_200_m_info_l_var <- ROC_plot() + 
  labs(caption = "l, h") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("Results/T1_n_100_K_200_m_m_MB.rds")
T1_K_200_m_info_m_var <- ROC_plot() +
  labs(caption = "l, l") + 
  theme(legend.position="none", legend.title=element_blank(), legend.text=element_text(size = 14),
        plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

panel_T1 <- plot_grid(T1_K_200_l_info_l_var, T1_K_200_l_info_m_var, T1_K_200_m_info_l_var, T1_K_200_m_info_m_var, nrow = 1)
title_1 <- ggdraw() + 
  draw_label("Chain", size = 20, fontface = 'bold')
T1 <- plot_grid(title_1, panel_T1, ncol = 1, rel_heights = c(0.2, 1.5)) 


# Type 3
results = read_rds("T3_n_100_K_200_l_l_MB.rds")
T3_K_200_l_info_l_var <- ROC_plot() + 
  labs(caption = "h, h") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("T3_n_100_K_200_l_m_MB.rds")
T3_K_200_l_info_m_var <- ROC_plot() +
  labs(caption = "h, l") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("T3_n_100_K_200_m_l_MB.rds")
T3_K_200_m_info_l_var <- ROC_plot() + 
  labs(caption = "l, h") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("Results/T3_n_100_K_200_m_m_MB.rds")
T3_K_200_m_info_m_var <- ROC_plot() +
  labs(caption = "l, l") + 
  theme(legend.position=c(0.2, 0.6), legend.title=element_blank(), legend.text=element_text(size = 14),
        plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold')) + 
  theme(legend.position= "none", legend.title=element_blank(), legend.text=element_text(size = 14),
        plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

panel_T3 <- plot_grid(T3_K_200_l_info_l_var, T3_K_200_l_info_m_var, T3_K_200_m_info_l_var, T3_K_200_m_info_m_var, nrow = 1)
title_3 <- ggdraw() + 
  draw_label("Random", size = 20, fontface = 'bold')
T3 <- plot_grid(title_3, panel_T3, ncol = 1, rel_heights = c(0.2, 1.5)) 

# Type 4
results = read_rds("Results/T4_n_100_K_200_l_l_MB.rds")
T4_K_200_l_info_l_var <- ROC_plot() + 
  labs(caption = "h, h") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("Results/T4_n_100_K_200_l_m_MB.rds")
T4_K_200_l_info_m_var <- ROC_plot() +
  labs(caption = "h, l") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("Results/T4_n_100_K_200_m_l_MB.rds")
T4_K_200_m_info_l_var <- ROC_plot() + 
  labs(caption = "l, h") + 
  theme(legend.position="none", plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

results = read_rds("Results/T4_n_100_K_200_m_m_MB.rds")
T4_K_200_m_info_m_var <- ROC_plot() +
  labs(caption = "l, l") + 
  theme(legend.position="bottom", legend.text=element_text(size = 14),
        plot.caption = element_text(hjust = 0.5, size = 18, face = 'bold'))

legend <- get_legend(T4_K_200_m_info_m_var + theme(legend.position = "bottom",
                     legend.text = element_text(size = 18), legend.title = element_text(size = 18)))
T4_K_200_m_info_m_var <- T4_K_200_m_info_m_var + theme(legend.position="None")

panel_T4 <- plot_grid(T4_K_200_l_info_l_var, T4_K_200_l_info_m_var, T4_K_200_m_info_l_var, T4_K_200_m_info_m_var, nrow = 1)
title_4 <- ggdraw() + 
  draw_label("Hub", size = 20, fontface = 'bold')
T4 <- plot_grid(title_4, panel_T4, ncol = 1, rel_heights = c(0.2, 1.5)) 
T4

T4_down <- plot_grid(T4, legend, ncol = 1, rel_heights = c(2, 0.3))

plot_grid(T1, T3, T4_down, nrow = 3, rel_heights = c(2, 2, 2.3))

rm("results")

save.image("ROC_MB.RData")


################################################ Run from here
library(tidyverse)
library(cowplot)
load("ROC_MB.RData")
# plot_grid(T1, T3, T4_down, nrow = 3, rel_heights = c(2, 2, 2.3))

# Remove legend:
plot_grid(T1, T3, T4, nrow = 3, rel_heights = c(2, 2, 2))


