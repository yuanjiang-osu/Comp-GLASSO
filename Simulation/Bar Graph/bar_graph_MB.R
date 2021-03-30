library(tidyverse)
library(reshape2)
library(cowplot)

results <- readRDS("results_StARS_sim.rds")
# Type 1
type1_K200_R <- results %>% filter(type == 1, K == 200) %>% group_by(method, var_read) %>% summarise(TP_mean = mean(TP), TP_sd = sd(TP), 
                                                                                                     l_bar = max(0, TP_mean-TP_sd), h_bar = TP_mean+TP_sd)
T1_R <- ggplot(type1_K200_R, aes(x = var_read, y = TP_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "Recall") +
  ylim(0, 1) +
  theme_light() + 
  theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 14))

type1_K200_P <- results %>% filter(type == 1, K == 200) %>% group_by(method, var_read) %>% summarise(PR_mean = mean(PR), PR_sd = sd(PR),
                                                                                                     l_bar = max(0, PR_mean-PR_sd), h_bar = PR_mean+PR_sd)
T1_P <- ggplot(type1_K200_P, aes(x = var_read, y = PR_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "Precision") +
  ylim(0, 1) +
  theme_light() + 
  theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 14))

type1_K200_F <- results %>% filter(type == 1, K == 200) %>% group_by(method, var_read) %>% summarise(F1_mean = mean(F1), F1_sd = sd(F1),
                                                                                                     l_bar = max(0, F1_mean-F1_sd), h_bar = F1_mean+F1_sd)
T1_F <- ggplot(type1_K200_F, aes(x = var_read, y = F1_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "F1") +
  ylim(0, 1) +
  theme_light() + 
  # theme(legend.position = c(0.75, 0.8),plot.title = element_text(hjust = 0.5, size = 16),
  #     axis.text=element_text(size = 14),
  #     legend.key = element_rect(colour = NA, fill = NA), legend.text=element_text(size=12)) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16),
        axis.text=element_text(size = 14),
        legend.key = element_rect(colour = NA, fill = NA), legend.text=element_text(size=12))

panel_T1 <- plot_grid(T1_R, T1_P, T1_F, nrow = 1)
title_1 <- ggdraw() + 
  draw_label("Chain", size = 20, fontface = 'bold')
x_label <- ggdraw() + 
  draw_label("Multinomial Variance, Sequencing Depth", size = 14)
# T1 <- plot_grid(title_1, panel_T1, x_label, ncol = 1, rel_heights = c(0.2, 1.5, 0.17)) 
T1 <- plot_grid(title_1, panel_T1, ncol = 1, rel_heights = c(0.2, 1.5)) 

# Type 3
type3_K200_R <- results %>% filter(type == 3, K == 200) %>% group_by(method, var_read) %>% summarise(TP_mean = mean(TP), TP_sd = sd(TP), 
                                                                                                     l_bar = max(0, TP_mean-TP_sd), h_bar = TP_mean+TP_sd)
T3_R <- ggplot(type3_K200_R, aes(x = var_read, y = TP_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "Recall") +
  ylim(0, 1) +
  theme_light() + 
  theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 14))

type3_K200_P <- results %>% filter(type == 3, K == 200) %>% group_by(method, var_read) %>% summarise(PR_mean = mean(PR), PR_sd = sd(PR),
                                                                                                     l_bar = max(0, PR_mean-PR_sd), h_bar = PR_mean+PR_sd)
T3_P <- ggplot(type3_K200_P, aes(x = var_read, y = PR_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "Precision") +
  ylim(0, 1) +
  theme_light() + 
  theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 14))

type3_K200_F <- results %>% filter(type == 3, K == 200) %>% group_by(method, var_read) %>% summarise(F1_mean = mean(F1), F1_sd = sd(F1),
                                                                                                     l_bar = max(0, F1_mean-F1_sd), h_bar = F1_mean+F1_sd)
T3_F <- ggplot(type3_K200_F, aes(x = var_read, y = F1_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "F1") +
  ylim(0, 1) +
  theme_light() + 
  # theme(legend.position = c(0.75, 0.8),plot.title = element_text(hjust = 0.5, size = 16),
  #       axis.text=element_text(size = 14),
  #       legend.key = element_rect(colour = NA, fill = NA), legend.text=element_text(size=12)) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16),
        axis.text=element_text(size = 14),
        legend.key = element_rect(colour = NA, fill = NA), legend.text=element_text(size=12))

panel_T3 <- plot_grid(T3_R, T3_P, T3_F, nrow = 1)
title_3 <- ggdraw() + 
  draw_label("Random", size = 20, fontface = 'bold')
x_label <- ggdraw() + 
  draw_label("Multinomial Variance, Sequencing Depth", size = 14)
# T3 <- plot_grid(title_3, panel_T3, x_label, ncol = 1, rel_heights = c(0.2, 1.5, 0.17)) 
T3 <- plot_grid(title_3, panel_T3, ncol = 1, rel_heights = c(0.2, 1.5)) 

# Type 4
type4_K200_R <- results %>% filter(type == 4, K == 200) %>% group_by(method, var_read) %>% summarise(TP_mean = mean(TP), TP_sd = sd(TP), 
                                                                                                     l_bar = max(0, TP_mean-TP_sd), h_bar = TP_mean+TP_sd)
T4_R <- ggplot(type4_K200_R, aes(x = var_read, y = TP_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "Recall") +
  ylim(0, 1) +
  theme_light() + 
  theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 14))

type4_K200_P <- results %>% filter(type == 4, K == 200) %>% group_by(method, var_read) %>% summarise(PR_mean = mean(PR), PR_sd = sd(PR),
                                                                                                     l_bar = max(0, PR_mean-PR_sd), h_bar = PR_mean+PR_sd)
T4_P <- ggplot(type4_K200_P, aes(x = var_read, y = PR_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", fill = "", title = "Precision") +
  ylim(0, 1) +
  theme_light() + 
  theme(legend.position = "None", plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 14))

type4_K200_F <- results %>% filter(type == 4, K == 200) %>% group_by(method, var_read) %>% summarise(F1_mean = mean(F1), F1_sd = sd(F1),
                                                                                                     l_bar = max(0, F1_mean-F1_sd), h_bar = F1_mean+F1_sd)
T4_F <- ggplot(type4_K200_F, aes(x = var_read, y = F1_mean, fill = method, width = 0.8)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar(aes(ymin= l_bar, ymax= h_bar),                  # Width of the error bars
                position=position_dodge()) + 
  labs(x = "", y = "", title = "F1") +
  scale_fill_discrete(name = "Method", labels = c('Compo-glasso', "glasso", "MB")) + 
  ylim(0, 1) +
  theme_light() + 
  # theme(legend.position = c(0.75, 0.8),plot.title = element_text(hjust = 0.5, size = 16),
  #       axis.text=element_text(size = 14),
  #       legend.key = element_rect(colour = NA, fill = NA), legend.text=element_text(size=12)) + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 16),
        axis.text=element_text(size = 14))

legend <- get_legend(T4_F + theme(legend.position = c(0.23, 0.7), legend.direction = "horizontal",
                                  legend.text=element_text(size = 18), legend.title = element_text(size = 18)))
T4_F <- T4_F + theme(legend.position="none")

panel_T4 <- plot_grid(T4_R, T4_P, T4_F, nrow = 1)
title_4 <- ggdraw() + 
  draw_label("Hub", size = 20, fontface = 'bold')
x_label <- ggdraw() + 
  draw_label("Multinomial Variance, Sequencing Depth", size = 14)
# T4 <- plot_grid(title_4, panel_T4, x_label, ncol = 1, rel_heights = c(0.2, 1.5, 0.17)) 
T4 <- plot_grid(title_4, panel_T4, ncol = 1, rel_heights = c(0.2, 1.5)) 

T4_down <- plot_grid(T4, legend, ncol = 1, rel_heights = c(2, 0.2))


# Remove legend:
all_3 = plot_grid(T1, T3, T4, nrow = 3, rel_heights = c(2, 2, 2))
all_3


