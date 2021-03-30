# Aim 1: simulation 
rm(list = ls())

simulation_comp <- function(seed)  
{
  gc()
  n.rep <- 100  ##of replicates
  num_cpu <- 10
  library(MASS)
  library(glasso)
  library(huge)
  source("CompoGlasso.R")
  
  n <- 100  # number of samples
  K <- 200 # number of OTUs (reference OTU not counted)
  length_rholist <- 70 # number of penalty parameter
  
  # Type of inverse-covariance matrix
  type <- 1 # chain
  # type <- 3 # random
  # type <- 4 # cluster
  
  # Size of sequencing depth
  seq_depth <- 'M' 
  # seq_depth <- 'L'
  
  # Size of multivariate normal variance
  # z_var <- 'M'
  z_var <- 'L'
  
  offset <- K + 1   ##20170406:To see if this influence a lot
  option <- 2 # Doesn't account for the cases when reference OTU = 0?

  
  results <- list() ##Setting list to save results
  
  ## data generation
  ###Setting for loop:
  for (n_pa in 1:(n.rep / num_cpu)) {
    temp <- generate_cov(K, type)
    Sigma.True <- temp$Sigma.True
    Omega.True <- temp$Omega.True
    
    avg <- 0
    mu <- rep(avg, K)
    if (z_var == 'M') err <- 1 else if (z_var == 'L') err <- 5
  
    z <- matrix(mvrnorm(n, mu, err * Sigma.True), nrow = n, ncol = K)

  p <- matrix(0, nrow = n, ncol = K + 1)
  p[, -(K + 1)] <- exp(z)/((apply(exp(z), 1, sum) + 1) %*% matrix(1, ncol = K))
  p[, (K + 1)] <- 1/(apply(exp(z), 1, sum) + 1)
  cat("Sum of Multinomial variance:", sum(apply(p, MARGIN = 1, FUN = var)), "\n")
  
  # Setting sequencing depth
  if (seq_depth == 'M') M <- runif(n, 20 * K, 40 * K) else if (seq_depth == 'L') M <- runif(n, 100 * K, 200 * K)
  cat("M's simulated \n")
  x <- matrix(0, n, K + 1)
  for(i in 1 : n)
  {
    x[i, ] <- rmultinom(1, size = M[i], prob = p[i, ])
  }
  cat("X's initiated \n")
  
  z.hat <- z_hat_offset(x, offset, option)
  Sigma.2 <- cov(z.hat)
  cat("Sigma 2 comupted \n")
  path.2 <- huge(Sigma.2, method = "glasso")
  cat("Glasso path created \n")
  
  rho.list.2 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 1000), length = length_rholist))
  
  # SPEIC-EASI gl
  cat("rho.list.2:", "\n")
  print(rho.list.2)

  n2.rho <- length(rho.list.2)
  TP.2_s <- rep(NA, n2.rho)
  FP.2_s <- rep(NA, n2.rho)
  ER.2_s <- rep(NA, n2.rho)
  RC.2_s <- rep(NA, n2.rho)
  
  path.2 <- huge(Sigma.2, method = "glasso", lambda = rho.list.2)
  Omegas.2 <- path.2$icov
  for(i.rho.2 in 1 : n2.rho)
  {
    cat("i.rho.2 = :", i.rho.2, "\n")
    Omega.2 <- as.matrix(Omegas.2[[i.rho.2]])
    TP.2_s[i.rho.2] <- sum((Omega.2 - diag(diag(Omega.2))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
    FP.2_s[i.rho.2] <- sum(Omega.2 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
    ER.2_s[i.rho.2] <- sum((Omega.2 - Omega.True)^2)  #Estimation errors
    RC.2_s[i.rho.2] <- sum((Omega.2 - diag(diag(Omega.2))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.2 != 0) - K) #Precision
  }
  
  # SPEIC-EASI mb
  rho.list.3 <- exp(seq(log(max(Sigma.2) * 2), log(max(Sigma.2) / 4000), length = length_rholist))
  
  cat("rho.list.3:", "\n")
  print(rho.list.3)
  
  n3.rho <- length(rho.list.3)
  TP.3_s <- rep(NA, n3.rho)
  FP.3_s <- rep(NA, n3.rho)
  ER.3_s <- rep(NA, n3.rho)
  RC.3_s <- rep(NA, n3.rho)
  
  path.3 <- huge(Sigma.2, method = "mb", lambda = rho.list.3)
  Omegas.adj.3 <- path.3$path
  print(Omegas.3)
  for(i.rho.3 in 1 : n3.rho)
  {
    cat("i.rho.3 = :", i.rho.3, "\n")
    Omega.3 <- as.matrix(Omegas.adj.3[[i.rho.3]]) + diag(K)
    Omegas.3[, , i.rho.3] = Omega.3
    TP.3_s[i.rho.3] <- sum((Omega.3 - diag(diag(Omega.3))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
    FP.3_s[i.rho.3] <- sum(Omega.3 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
    ER.3_s[i.rho.3] <- sum((Omega.3 - Omega.True)^2)  #Estimation errors
    RC.3_s[i.rho.3] <- sum((Omega.3 - diag(diag(Omega.3))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.3 != 0) - K) #Precision
  }
  
  # Compo-glasso
  rho.list.1 = rho.list.2
  n1.rho <- length(rho.list.1)
  TP.1_s <- rep(NA, n1.rho)
  FP.1_s <- rep(NA, n1.rho)
  ER.1_s <- rep(NA, n1.rho)
  RC.1_s <- rep(NA, n1.rho)
  
  Omegas.1 <- array(0, dim = c(K, K, n1.rho))
  flag.1 <- rep(0, n1.rho)
  Sigmas_rho <- list()
  Omegas_rho <- list()
  for(i.rho.1 in 1 : n1.rho)    
  {
    cat("i.rho.1 =", i.rho.1, "\n")
    rho <- rho.list.1[i.rho.1]
    Sigmas_iter <- list()
    Omegas_iter <- list()
    
    Sigma.0 <- cov(z.hat)
    if (i.rho.1 == 1) {
      Omega.0 <- as.matrix(Omegas.2[[length_rholist]])
      }
    else if (i.rho.1 > 1) {
      Omega.0 <- as.matrix(Omegas.2[[1]])
      }
    Omega.1 <- matrix(0, K, K)
    iter <- 0
    z.0 <- z.hat
    z.1 <- matrix(0, n, K)
    mu.0 <- apply(z.0, 2, mean)
    z.start <- z.hat
    z.end <- matrix(0, n, K)
    while (mean((Omega.0 - Omega.1) ^ 2) > 0.0000001 || mean((z.start - z.end) ^ 2) > 0.1)
    {
      cat("iter = ", iter + 1, "mean((Omega.0 - Omega.1) ^ 2) = ", mean((Omega.0 - Omega.1) ^ 2), "\n")
      if (iter != 0) {
        Omega.0 <- Omega.1
        }
      
      if (mean((z.start - z.end) ^ 2) > 0.1) {
        if (iter != 0) {
        z.start <- z.end
        }
        iter <- iter + 1
        # cat("iter = ", iter, "\n")
        z.0 <- z.start
        mu.0 <- apply(z.0, 2, mean)
        z.1 <- matrix(0, n, K)
        z.iter <- 0
        while (mean((z.0 - z.1) ^ 2) > 0.0000001){
          if (z.iter != 0) {
            z.0 <- z.1
          }
          z.iter <- z.iter + 1
          for (j in 1:n) {
            dipi <- M[j] * exp(z.0[j,]) / (t(rep(1, K)) %*% exp(z.0[j,]) + 1) - x[j, 1:K] +  as.vector(Omega.0 %*% (z.0[j,] - mu.0))
            tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - 
              M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)) ^ 2 + Omega.0
            z.1[j,] <- z.0[j,] - solve(tripi) %*% dipi
          }
          cat("z iteration = ", z.iter, "mean square difference in z's = ", mean((z.0 - z.1) ^ 2), "\n")
        }
      z.end <- z.1
      Sigma.1 <- cov(z.end)
    
      mod <- huge(x = Sigma.1, lambda = rho, method = "glasso")
    
      Omega.1 = as.matrix(mod$icov[[1]])
  
      Sigmas_iter[[iter]] <- Sigma.1
      Omegas_iter[[iter]] <- Omega.1
      }
    }
    
    Sigmas_rho[[i.rho.1]] <- Sigmas_iter
    Omegas_rho[[i.rho.1]] <- Omegas_iter
  
    Omegas.1[, , i.rho.1] <- Omega.1
    TP.1_s[i.rho.1] <- sum((Omega.1 - diag(diag(Omega.1))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.True != 0) - K)
    FP.1_s[i.rho.1] <- sum(Omega.1 != 0 & Omega.True == 0)/sum(Omega.True == 0) #False positive rates
    ER.1_s[i.rho.1] <- sum((Omega.1 - Omega.True)^2)  #Estimation errors
    RC.1_s[i.rho.1] <- sum((Omega.1 - diag(diag(Omega.1))) & (Omega.True- diag(diag(Omega.True))))/(sum(Omega.1 != 0) - K)
  }

  results[[n_pa]] <- list(Omega.True, Omegas.1, Omegas.2, Omegas.3,
                          TP.1_s, TP.2_s, TP.3_s,  
                          FP.1_s, FP.2_s, FP.3_s, 
                          ER.1_s, ER.2_s, ER.3_s, 
                          RC.1_s, RC.2_s, RC.3_s,
                          rho.list.1, rho.list.2, rho.list.3)
  names(results[[n_pa]]) <- c("Omega.True", "Omegas.1", "Omegas.2", "Omegas.3",
                              "TP.1_s", "TP.2_s", "TP.3_s", 
                              "FP.1_s", "FP.2_s", "FP.3_s", 
                              "ER.1_s", "ER.2_s", "ER.3_s", 
                              "RC.1_s", "RC.2_s", "RC.3_s", 
                              "rho.list.1", "rho.list.2", "rho.list.3")
  }
  return(results)
}

### Basic parellel computing:
start_time <- proc.time()
library(snow)
seed_list <- sample(1:n.rep, num_cpu, replace = FALSE)
cluster <- makeCluster(num_cpu, type = "SOCK", outfile = "")
res_list <- clusterApply(cluster, seed_list, simulation_comp)
stopCluster(cluster)

warnings()
saveRDS(res_list, "T1_n_100_K_200_m_l_MB.rds") 
end_time <- proc.time()
print(end_time - start_time)
