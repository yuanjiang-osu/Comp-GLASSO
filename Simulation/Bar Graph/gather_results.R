library(tidyverse)
library(reshape2)
# library(cowplot)

# Count the number of settings/method combination
num_n <- 1
num_K <- 2
num_method <- 3
num_var <- 2
num_seq_depth <- 2
num_type <- 4
num_comb <- num_n * num_K * num_method * num_var * num_seq_depth * num_type
n_rep <- 50

# Initialization
TP <- rep(0, num_comb * n_rep)
FP <- rep(0, num_comb * n_rep)
PR <- rep(0, num_comb * n_rep)
F1 <- rep(0, num_comb * n_rep)
method <- rep(0, num_comb * n_rep)
multi_var <- rep(0, num_comb * n_rep)
seq_depth <- rep(0, num_comb * n_rep)
K <- rep(0, num_comb * n_rep)
n <- rep(0, num_comb * n_rep)
type <- rep(0, num_comb * n_rep)

counter <- 0
# T1, n = 100, K = 200, m, m
results <- read_rds("T1_n_100_K_200_m_m_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

# T1, n = 100, K = 200, m, l
results <- read_rds("T1_n_100_K_200_m_l_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

# T1, n = 100, K = 200, l, m
results <- read_rds("T1_n_100_K_200_l_m_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

# T1, n = 100, K = 200, l, l
results <- read_rds("T1_n_100_K_200_l_l_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 1
}
counter <- counter + 1

# T3, n = 100, K = 200, m, m
results <- read_rds("T3_n_100_K_200_m_m_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

# T3, n = 100, K = 200, m, l
results <- read_rds("T3_n_100_K_200_m_l_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

# T3, n = 100, K = 200, l, m
results <- read_rds("T3_n_100_K_200_l_m_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

# T3, n = 100, K = 200, l, l
results <- read_rds("T3_n_100_K_200_l_l_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 3
}
counter <- counter + 1


# T4, n = 100, K = 200, m, m
results <- read_rds("T4_n_100_K_200_m_m_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

# T4, n = 100, K = 200, m, l
results <- read_rds("T4_n_100_K_200_m_l_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "l"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

# T4, n = 100, K = 200, l, m
results <- read_rds("T4_n_100_K_200_l_m_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "l"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

# T4, n = 100, K = 200, l, l
results <- read_rds("T4_n_100_K_200_l_l_StARS_MB.rds")
for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.1
  FP[i + counter * n_rep] <- results[[i]]$FP.1
  PR[i + counter * n_rep] <- results[[i]]$PR.1
  F1[i + counter * n_rep] <- results[[i]]$F1.1
  method[i + counter * n_rep] <- "Compo-glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.2
  FP[i + counter * n_rep] <- results[[i]]$FP.2
  PR[i + counter * n_rep] <- results[[i]]$PR.2
  F1[i + counter * n_rep] <- results[[i]]$F1.2
  method[i + counter * n_rep] <- "glasso"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

for (i in 1:n_rep) {
  TP[i + counter * n_rep] <- results[[i]]$TP.3
  FP[i + counter * n_rep] <- results[[i]]$FP.3
  PR[i + counter * n_rep] <- results[[i]]$PR.3
  F1[i + counter * n_rep] <- results[[i]]$F1.3
  method[i + counter * n_rep] <- "mb"
  multi_var[i + counter * n_rep] <- "h"
  seq_depth[i + counter * n_rep] <- "h"
  K[i + counter * n_rep] <- 200
  n[i + counter * n_rep] <- 100
  type[i + counter * n_rep] <- 4
}
counter <- counter + 1

result_collection <- data.frame(n = n[1:1800], K = K[1:1800], type = type[1:1800], multi_var = multi_var[1:1800], seq_depth = seq_depth[1:1800],
                                method = method[1:1800], TP = TP[1:1800], FP = FP[1:1800], PR = PR[1:1800], F1 = F1[1:1800])
result_collection
results <- result_collection %>% unite(var_read, c(multi_var, seq_depth), sep = ", ", remove = FALSE)
saveRDS(results, file = "results_StARS_sim.rds")
