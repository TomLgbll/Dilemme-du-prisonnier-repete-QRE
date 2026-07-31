#-------------------------------------------------------------------------------
#                        Monte Carlo 
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Étude Monte Carlo avec hétérogénéité des préférences
# (Les paramètres R, S, T, P sont déjà définis)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Fonction U avec sécurité intégrée
# ------------------------------------------------------------------------------

# Fonction d'utilité des agents

U <- function(x, r){
  # Sécurité : x doit être strictement positif
  x <- max(x, 1e-6)
  
  # Sécurité : r ne doit pas être trop proche de 1
  if(abs(r - 1) < 1e-6){
    return(log(x))
  } else {
    return(x^(1 - r) / (1 - r))
  }
}

# ------------------------------------------------------------------------------
# Simulation avec apprentissage (on stocke les croyances)
# ------------------------------------------------------------------------------
# Vrais paramètres 
mu1_true <- 0.7
mu2_true <- 0.5
r1_true  <- 0.3
r2_true  <- 0.6

simulate_data_hetero <- function(mu1, mu2, r1, r2, n = 2000){
  
  p <- 0.5
  q <- 0.5
  out <- data.frame() 
  
  for(t in 1:n){
    
    EU1_C <- q * U(R, r1) + (1 - q) * U(S, r1)
    EU1_D <- q * U(T, r1) + (1 - q) * U(P, r1)
    
    EU2_C <- p * U(R, r2) + (1 - p) * U(S, r2)
    EU2_D <- p * U(T, r2) + (1 - p) * U(P, r2)
    
    p1 <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
    p2 <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
    
    p1 <- max(min(p1, 1 - 1e-8), 1e-8)
    p2 <- max(min(p2, 1 - 1e-8), 1e-8)
    
    a1 <- rbinom(1, 1, p1)
    a2 <- rbinom(1, 1, p2)
    
    # Stockage des croyances AVANT mise à jour
    out <- rbind(out, data.frame(
      a1 = a1, 
      a2 = a2,
      p_before = p,
      q_before = q
    ))
    
    # Mise à jour des croyances
    p <- 0.9 * p + 0.1 * a1
    q <- 0.9 * q + 0.1 * a2
  }
  
  return(out)   
}

# ------------------------------------------------------------------------------
# Estimation STATIQUE 
# ------------------------------------------------------------------------------
estimate_model_static_hetero <- function(dat_hetero){   
  
  loglik <- function(par){
    
    mu1 <- par[1]; mu2 <- par[2]
    r1  <- par[3]; r2  <- par[4]
    
    # Sécurité : paramètres hors bornes ou invalides
    if(any(par <= 0) || any(is.na(par)) || any(is.infinite(par)) ||
       r1 < 0.01 || r1 > 0.99 || r2 < 0.01 || r2 > 0.99 ||
       mu1 < 0.01 || mu1 > 3 || mu2 < 0.01 || mu2 > 3){
      return(1e10)
    }
    
    p <- 0.5
    q <- 0.5
    
    for(i in 1:100){
      
      EU1_C <- q * U(R, r1) + (1 - q) * U(S, r1)
      EU1_D <- q * U(T, r1) + (1 - q) * U(P, r1)
      
      EU2_C <- p * U(R, r2) + (1 - p) * U(S, r2)
      EU2_D <- p * U(T, r2) + (1 - p) * U(P, r2)
      
      # Sécurité : valeurs invalides
      if(any(is.na(c(EU1_C, EU1_D, EU2_C, EU2_D))) ||
         any(is.infinite(c(EU1_C, EU1_D, EU2_C, EU2_D)))){
        return(1e10)
      }
      
      p_best <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
      q_best <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
      
      p_new <- 0.5 * p_best + 0.5 * p
      q_new <- 0.5 * q_best + 0.5 * q
      
      if(abs(p_new - p) < 1e-8 && abs(q_new - q) < 1e-8){
        p <- p_new
        q <- q_new
        break
      }
      
      p <- p_new
      q <- q_new
    }
    
    p <- max(min(p, 1 - 1e-8), 1e-8)
    q <- max(min(q, 1 - 1e-8), 1e-8)
    
    # Log-vraisemblance
    ll <- sum(dat_hetero$a1 * log(p) + (1 - dat_hetero$a1) * log(1 - p)) +
      sum(dat_hetero$a2 * log(q) + (1 - dat_hetero$a2) * log(1 - q))
    
    if(is.na(ll) || is.infinite(ll) || is.nan(ll)) return(1e10)
    
    return(-ll)
  }
  
  # Bornes très serrées pour éviter les problèmes
  opt <- optim(c(0.6, 0.4, 0.3, 0.3), loglik,
               method = "L-BFGS-B",
               lower = c(0.1, 0.1, 0.05, 0.05),
               upper = c(2, 2, 0.8, 0.8),
               control = list(maxit = 200, trace = 0))
  
  return(opt$par)
}


# ------------------------------------------------------------------------------
# Monte Carlo
# ------------------------------------------------------------------------------
Nsim <- 50
n <- 2000

results <- matrix(0, nrow = Nsim, ncol = 4)

for(s in 1:Nsim){
  cat("Simulation", s, "/", Nsim, "\r")
  dat_mc <- simulate_data_hetero(mu1_true, mu2_true, r1_true, r2_true, n = n)   
  est <- estimate_model_static_hetero(dat_mc)                                    
  results[s, ] <- est
}

colnames(results) <- c("mu1_hat", "mu2_hat", "r1_hat", "r2_hat")

cat("\n==============================\n")
cat("      RESULTATS MONTE CARLO\n")
cat("==============================\n\n")

cat("Moyennes des estimateurs :\n")
print(colMeans(results))

cat("\nEcarts-types des estimateurs :\n")
print(apply(results, 2, sd))

#-------------------------------------------------------------------------------
#                  Monte Carlo : effet de la taille d'échantillon n
#-------------------------------------------------------------------------------
n_values <- c(20, 50, 100, 200, 500)
results_n <- data.frame()

for (n in n_values) {
  
  mu1_hat <- numeric(Nsim)
  mu2_hat <- numeric(Nsim)
  r1_hat  <- numeric(Nsim)
  r2_hat  <- numeric(Nsim)
  
  for (s in 1:Nsim) {
    
    dat_mc <- simulate_data_hetero(mu1_true, mu2_true, r1_true, r2_true, n = n)
    est <- estimate_model_static_hetero(dat_mc)
    
    mu1_hat[s] <- est[1]
    mu2_hat[s] <- est[2]
    r1_hat[s]  <- est[3]
    r2_hat[s]  <- est[4]
  }
  
  results_n <- rbind(results_n, data.frame(
    n = n,
    mu1_mean = mean(mu1_hat),
    mu2_mean = mean(mu2_hat),
    r1_mean = mean(r1_hat),
    r2_mean = mean(r2_hat),
    mu1_sd = sd(mu1_hat),
    mu2_sd = sd(mu2_hat),
    r1_sd = sd(r1_hat),
    r2_sd = sd(r2_hat)
  ))
}
cat("\n==============================\n")
cat(" EFFET DE LA TAILLE n\n")
cat("==============================\n\n")

print(results_n)

# Graphique des écarts-types en fonction de n
plot(results_n$n, results_n$mu1_sd, 
     type = "b",
     pch = 19, 
     col = "blue",
     xlab = "Taille d'échantillon n",
     ylab = "Écart-type des estimateurs",
     ylim = c(0, max(results_n$mu1_sd, results_n$mu2_sd, 
                     results_n$r1_sd, results_n$r2_sd) * 1.2),
     main = "Effet de n sur la précision des estimateurs")

lines(results_n$n, results_n$mu2_sd, 
      type = "b",
      pch = 17, 
      col = "cyan")

lines(results_n$n, results_n$r1_sd, 
      type = "b",
      pch = 15, 
      col = "red")

lines(results_n$n, results_n$r2_sd, 
      type = "b",
      pch = 18, 
      col = "darkred")

legend("topright",
       legend = c(expression(hat(mu)[1]), expression(hat(mu)[2]),
                  expression(hat(r)[1]), expression(hat(r)[2])),
       col = c("blue", "cyan", "red", "darkred"),
       pch = c(19, 17, 15, 18),
       lty = 1)


#-------------------------------------------------------------------------------
#                  Monte Carlo : effet du paramètre de bruit mu
#-------------------------------------------------------------------------------
mu1_values <- c(0.2, 0.5, 0.7, 1, 1.5, 2)
results_mu <- data.frame()

for (mu1 in mu1_values) {
  
  # On fixe les autres paramètres
  mu2 <- 0.5
  r1  <- 0.3
  r2  <- 0.6
  
  mu1_hat <- numeric(Nsim)
  mu2_hat <- numeric(Nsim)
  r1_hat  <- numeric(Nsim)
  r2_hat  <- numeric(Nsim)
  
  for (s in 1:Nsim) {
    
    dat_mc <- simulate_data_hetero(mu1, mu2, r1, r2, n = 2000)
    est <- estimate_model_static_hetero(dat_mc)
    
    mu1_hat[s] <- est[1]
    mu2_hat[s] <- est[2]
    r1_hat[s]  <- est[3]
    r2_hat[s]  <- est[4]
  }
  
  results_mu <- rbind(results_mu, data.frame(
    mu1_true = mu1,
    mu1_mean = mean(mu1_hat),
    mu2_mean = mean(mu2_hat),
    r1_mean = mean(r1_hat),
    r2_mean = mean(r2_hat),
    mu1_sd = sd(mu1_hat),
    mu2_sd = sd(mu2_hat),
    r1_sd = sd(r1_hat),
    r2_sd = sd(r2_hat)
  ))
}
cat("\n==============================\n")
cat(" EFFET DU PARAMETRE DE BRUIT μ\n")
cat("==============================\n\n")
print(results_mu)


# =========================
# Graphique : effet de mu1 sur la précision des 4 estimateurs
# =========================
plot(results_mu$mu1_true, results_mu$mu1_sd,
     type = "b",
     pch = 19,
     col = "blue",
     xlab = expression(mu[1]~"(vrai)"),
     ylab = "Écart-type des estimateurs",
     ylim = c(0, max(results_mu$mu1_sd, results_mu$mu2_sd,
                     results_mu$r1_sd, results_mu$r2_sd) * 1.2),
     main = "Effet du bruit sur la précision des estimateurs")

lines(results_mu$mu1_true, results_mu$mu2_sd,
      type = "b",
      pch = 17,
      col = "cyan")

lines(results_mu$mu1_true, results_mu$r1_sd,
      type = "b",
      pch = 15,
      col = "red")

lines(results_mu$mu1_true, results_mu$r2_sd,
      type = "b",
      pch = 18,
      col = "darkred")

legend("topleft",
       legend = c(expression(hat(mu)[1]), expression(hat(mu)[2]),
                  expression(hat(r)[1]), expression(hat(r)[2])),
       col = c("blue", "cyan", "red", "darkred"),
       pch = c(19, 17, 15, 18),
       lty = 1)