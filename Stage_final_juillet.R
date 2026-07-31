#-------------------------------------------------------------------------------
#            Simulation exploratoire des données
#-------------------------------------------------------------------------------
set.seed(42)

nbr_tour = 100        # nombre de périodes
pairs = 10

# paramètres du jeu (notation standard du dilemme du prisonnier)
R = 3          # récompense mutuelle
S = 0          # sucker payoff
T = 4         # tentation de dévier 
P = 1          # punition mutuelle

data_list = list()
k = 1

for(p in 1:pairs){
  
  # hétérogénéité des comportements
  p_j1 = runif(1, 0.45, 0.70)
  p_j2 = runif(1, 0.60, 0.85)
  
  for(t in 1:nbr_tour){   # boucle sur les périodes
    
    j1_C = rbinom(1, 1, p_j1)
    j2_C = rbinom(1, 1, p_j2)
    
    # actions finales
    action_joueur1 = ifelse(j1_C == 1, "C", "D")
    action_joueur2 = ifelse(j2_C == 1, "C", "D")
    
    # gains joueur 1
    gain1 = ifelse(action_joueur1=="C" & action_joueur2=="C", R,
                   ifelse(action_joueur1=="C" & action_joueur2=="D", S,
                          ifelse(action_joueur1=="D" & action_joueur2=="C", T, P)))
    
    # gains joueur 2
    gain2 = ifelse(action_joueur1=="C" & action_joueur2=="C", R,
                   ifelse(action_joueur1=="C" & action_joueur2=="D", T,
                          ifelse(action_joueur1=="D" & action_joueur2=="C", S, P)))
    
    data_list[[k]] = data.frame(
      pair = p,
      tour = t,
      action_joueur1 = action_joueur1,
      action_joueur2 = action_joueur2,
      gain1 = gain1,
      gain2 = gain2
    )
    
    k = k + 1
  }
}

data = do.call(rbind, data_list)

View(data)
#-------------------------------------------------------------------------------
#               Fréquence de coopérations
#-------------------------------------------------------------------------------

prop_C_j1 = mean(data$action_joueur1 == "C")
prop_C_j2 = mean(data$action_joueur2 == "C")

cat("Proportion coopération joueur 1 :", prop_C_j1, "\n")
cat("Proportion coopération joueur 2 :", prop_C_j2, "\n")

aggregate(
  cbind(
    mean_j1 = as.numeric(data$action_joueur1 == "C"),
    mean_j2 = as.numeric(data$action_joueur2 == "C")
  ),
  by = list(pair = data$pair),
  FUN = mean
)

#-------------------------------------------------------------------------------
#                          QRE + RESOLUTION
#-------------------------------------------------------------------------------

# paramètres du jeu
R = 3
S = 0
T = 4
P = 1

# fonction à minimiser
objective = function(x, mu) {
  p = x[1]
  q = x[2]
  
  logit_p = exp((q*R + (1-q)*S)/mu) / 
    (exp((q*R + (1-q)*S)/mu) + exp((q*T + (1-q)*P)/mu))
  
  logit_q = exp((p*R + (1-p)*S)/mu) / 
    (exp((p*R + (1-p)*S)/mu) + exp((p*T + (1-p)*P)/mu))
  
  F1 = (p - logit_p)^2
  F2 = (q - logit_q)^2
  
  return(F1 + F2)
}

# grille de mu
mu_values = c(0.05, 0.1, 0.3, 0.5, 0.7, 1, 2)

# stockage des résultats
results = data.frame(mu = numeric(),
                     p = numeric(),
                     q = numeric())

for (m in mu_values) {
  # optimisation
  opt = optim(par = c(0.5, 0.5), 
              fn = objective, 
              mu = m,
              method = "L-BFGS-B",
              lower = c(0, 0),
              upper = c(1, 1))
  
  results = rbind(results,
                  data.frame(mu = m,
                             p = opt$par[1],
                             q = opt$par[2]))
}

print(results)

#------------------------------------------------------------
# GRAPH : effet de mu sur le QRE
#------------------------------------------------------------
# Graphique
plot(results$mu, results$p,
     type = "b",
     pch = 19,
     col = "blue",
     ylim = c(0, 1),
     xlab = expression(mu),
     ylab = "Probabilité de coopération",
     main = "Effet de μ sur l'équilibre QRE")

lines(results$mu, results$q,
      type = "b",
      pch = 17,
      col = "red")

legend("bottomright",
       legend = c("p* (joueur 1)", "q* (joueur 2)"),
       col = c("blue", "red"),
       pch = c(19, 17),
       lty = 1)

#-------------------------------------------------------------------------------
#               QRE AVEC AVERSION AU RISQUE (CRRA)
#-------------------------------------------------------------------------------
set.seed(42)

nbr_tour = 100        # nombre de périodes
pairs = 10

# paramètres du jeu (notation standard du dilemme du prisonnier)
R = 3          # récompense mutuelle
S = 0          # sucker payoff
T = 4         # tentation de dévier 
P = 1          # punition mutuelle

data_list = list()
k = 1



U = function(x, r){
  x = max(x, 1e-8) #car log(x) avec x=0 impossible or S=0 
  if(r == 1){
    return(log(x))
  } else {
    return(x^(1 - r) / (1 - r))
  }
}

for(p in 1:pairs){
  
  mu_j1 = runif(1, 0.2, 1)
  mu_j2 = runif(1, 0.2, 1)
  
  r_j1 = runif(1, 0, 0.8)
  r_j2 = runif(1, 0, 0.8)
  
  j1_C_prob = 0.5
  j2_C_prob = 0.5
  
  for(t in 1:nbr_tour){
    
    EU1_C = j2_C_prob * U(R, r_j1) + (1 - j2_C_prob) * U(S, r_j1)
    EU1_D = j2_C_prob * U(T, r_j1) + (1 - j2_C_prob) * U(P, r_j1)
    
    EU2_C = j1_C_prob * U(R, r_j2) + (1 - j1_C_prob) * U(S, r_j2)
    EU2_D = j1_C_prob * U(T, r_j2) + (1 - j1_C_prob) * U(P, r_j2)
    
    p1 = exp(EU1_C / mu_j1) / (exp(EU1_C / mu_j1) + exp(EU1_D / mu_j1))
    p2 = exp(EU2_C / mu_j2) / (exp(EU2_C / mu_j2) + exp(EU2_D / mu_j2))
    
    j1_C = rbinom(1, 1, p1)
    j2_C = rbinom(1, 1, p2)
    
    j1_C_prob = 0.9 * j1_C_prob + 0.1 * j1_C
    j2_C_prob = 0.9 * j2_C_prob + 0.1 * j2_C
    
    action_joueur1 = ifelse(j1_C == 1, "C", "D")
    action_joueur2 = ifelse(j2_C == 1, "C", "D")
    
    gain1 = ifelse(action_joueur1=="C" & action_joueur2=="C", R,
                   ifelse(action_joueur1=="C" & action_joueur2=="D", S,
                          ifelse(action_joueur1=="D" & action_joueur2=="C", T, P)))
    
    gain2 = ifelse(action_joueur1=="C" & action_joueur2=="C", R,
                   ifelse(action_joueur1=="C" & action_joueur2=="D", T,
                          ifelse(action_joueur1=="D" & action_joueur2=="C", S, P)))
    
    data_list[[k]] = data.frame(
      pair = p,
      tour = t,
      action_joueur1 = action_joueur1,
      action_joueur2 = action_joueur2,
      gain1 = gain1,
      gain2 = gain2,
      mu1 = mu_j1,
      mu2 = mu_j2,
      r1 = r_j1,
      r2 = r_j2
    )
    
    k = k + 1
  }
}

data = do.call(rbind, data_list)

View(data)


# === 1. Paramètres de la paire 1 ===
pair1 <- subset(data, pair == 1)

mu1_pair1 <- unique(pair1$mu1)[1]  # Premier mu1 de la paire 1
mu2_pair1 <- unique(pair1$mu2)[1]  # Premier mu2 de la paire 1
r1_pair1  <- unique(pair1$r1)[1]   # Premier r1 de la paire 1
r2_pair1  <- unique(pair1$r2)[1]   # Premier r2 de la paire 1

cat("=== Paramètres de la paire 1 ===\n")
cat("Joueur 1 : r1 =", round(r1_pair1, 3), ", mu1 =", round(mu1_pair1, 3), "\n")
cat("Joueur 2 : r2 =", round(r2_pair1, 3), ", mu2 =", round(mu2_pair1, 3), "\n")

# Probabilités observées (sans arrondi)
p_obs <- mean(pair1$action_joueur1 == "C")
q_obs <- mean(pair1$action_joueur2 == "C")

# Affichage sans arrondi
cat("p* (joueur 1) =", p_obs, "\n")
cat("q* (joueur 2) =", q_obs, "\n")

#-------------------------------------------------------------------------------
#               QRE AVEC AVERSION AU RISQUE (CRRA)
#            Effet de l'hétérogénéité de l'aversion au risque
#-------------------------------------------------------------------------------
# paramètres fixes
r1 = 0.2
mu = 0.7

# valeurs de r2 à tester
r2_values = c(0, 0.2, 0.5, 0.8, 0.99)

# stockage des résultats
results_r2 = data.frame(r2 = numeric(), p = numeric(), q = numeric())

for (r2 in r2_values) {
  
  p = 0.5
  q = 0.5
  
  for(i in 1:1000){
    
    EU1_C = q * U(R, r1) + (1 - q) * U(S, r1)
    EU1_D = q * U(T, r1) + (1 - q) * U(P, r1)
    
    EU2_C = p * U(R, r2) + (1 - p) * U(S, r2)
    EU2_D = p * U(T, r2) + (1 - p) * U(P, r2)
    
    p_best = exp(EU1_C / mu) / (exp(EU1_C / mu) + exp(EU1_D / mu))
    q_best = exp(EU2_C / mu) / (exp(EU2_C / mu) + exp(EU2_D / mu))
    
    alpha = 0.9
    p_new = alpha * p_best + (1 - alpha) * p
    q_new = alpha * q_best + (1 - alpha) * q
    
    if(abs(p_new - p) < 1e-8 & abs(q_new - q) < 1e-8) break
    
    p = p_new
    q = q_new
  }
  
  results_r2 = rbind(results_r2, data.frame(r2 = r2, p = p, q = q))
}

print(results_r2)

# graphique
plot(results_r2$r2, results_r2$p, 
     type = "b", 
     pch = 19, 
     col = "blue",
     xlab = "r2 (aversion au risque du joueur 2)",
     ylab = "Probabilité de coopération",
     ylim = c(0, 0.4),
     main = "Effet de l'hétérogénéité des préférences (r1 = 0.2, μ = 0.7)")

lines(results_r2$r2, results_r2$q, 
      type = "b", 
      pch = 17, 
      col = "red")

legend("topright", 
       legend = c("p* (joueur 1)", "q* (joueur 2)"), 
       col = c("blue", "red"), 
       pch = c(19, 17), 
       lty = 1)

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
  cat("Simulation", s, "sur", Nsim, "\n")
  dat_mc <- simulate_data_hetero(mu1_true, mu2_true, r1_true, r2_true, n = n)   
  est <- estimate_model_static_hetero(dat_mc)                                    
  results[s, ] <- est
}

colnames(results) <- c("mu1_hat", "mu2_hat", "r1_hat", "r2_hat")

cat("Moyennes empiriques :\n")
print(colMeans(results))
cat("\nÉcarts-types empiriques :\n")
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

#-------------------------------------------------------------------------------
#                         III Estimation
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         III Estimation - MLE
#-------------------------------------------------------------------------------

set.seed(123)
dat_hetero <- simulate_data_hetero(mu1_true, mu2_true, r1_true, r2_true, n = 2000)


# ------------------------------------------------------------------------------
# loglik_QRE_CRRA : log-vraisemblance TOTALE, ell_N(theta) = sum_{i=1}^N ln f(y_i;theta)
# Cette fonction boucle sur les N observations et SOMME les contributions.
# C'est celle qu'on passe a optim() : optim minimise -ell_N(theta), donc trouve
# theta_hat = argmax_theta ell_N(theta) = le MLE.
# C'est aussi cette fonction qui, via optim(hessian = TRUE), produira la matrice
# opt$hessian utilisee plus bas pour construire H_N (Sample Hessian estimator).
# ------------------------------------------------------------------------------
loglik_QRE_CRRA <- function(params, dat_hetero){
  
  mu1 <- params[1]
  mu2 <- params[2]
  r1  <- params[3]
  r2  <- params[4]
  
  # Sécurité : paramètres invalides
  if(any(params <= 0) || any(is.na(params)) || any(is.infinite(params))){
    return(1e10)
  }
  
  # Sécurité : r ne doit pas être trop proche de 1 (problème numérique)
  if(abs(r1 - 1) < 1e-6 || abs(r2 - 1) < 1e-6){
    return(1e10)
  }
  
  # Croyances empiriques
  p_emp <- mean(dat_hetero$a1 == 1)
  q_emp <- mean(dat_hetero$a2 == 1)
  
  loglik <- 0
  
  for(k in 1:nrow(dat_hetero)){    # <-- boucle sur les N observations : c'est le "sum_i" de ell_N(theta)
    
    a1 <- dat_hetero$a1[k]
    a2 <- dat_hetero$a2[k]
    
    # ---- Joueur 1 ----
    EU1_C <- q_emp * U(R, r1) + (1 - q_emp) * U(S, r1)
    EU1_D <- q_emp * U(T, r1) + (1 - q_emp) * U(P, r1)
    
    # Sécurité : valeurs infinies ou NA
    if(any(is.na(c(EU1_C, EU1_D))) || any(is.infinite(c(EU1_C, EU1_D)))){
      return(1e10)
    }
    
    p1_C <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
    p1_C <- max(min(p1_C, 1 - 1e-10), 1e-10)
    
    # ---- Joueur 2 ----
    EU2_C <- p_emp * U(R, r2) + (1 - p_emp) * U(S, r2)
    EU2_D <- p_emp * U(T, r2) + (1 - p_emp) * U(P, r2)
    
    if(any(is.na(c(EU2_C, EU2_D))) || any(is.infinite(c(EU2_C, EU2_D)))){
      return(1e10)
    }
    
    p2_C <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
    p2_C <- max(min(p2_C, 1 - 1e-10), 1e-10)
    
    # Log-vraisemblance (avec sécurités)
    if(a1 == 1){
      loglik <- loglik + log(max(p1_C, 1e-15))    # <-- accumulation = la somme sum_i ln f(y_i;theta)
    } else {
      loglik <- loglik + log(max(1 - p1_C, 1e-15))
    }
    
    if(a2 == 1){
      loglik <- loglik + log(max(p2_C, 1e-15))
    } else {
      loglik <- loglik + log(max(1 - p2_C, 1e-15))
    }
    
    # Sécurité : si loglik devient infini
    if(is.infinite(loglik) || is.na(loglik)){
      return(1e10)
    }
  }
  
  return(-loglik)   # <-- -ell_N(theta) : c'est CE signe et CETTE somme qu'optim va minimiser
}

#-------------------------------------------------------------------------------
#                         III - 1  Log-vraisemblance INDIVIDUELLE ln f(y_i;theta), i = 1..N (pour les scores)
# renvoie un vecteur de longueur N 
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# loglik_i_QRE_CRRA : contribution individuelle ln f(y_i;theta), SANS sommer.
# Meme calcul que loglik_QRE_CRRA mais on garde chaque i separement (vecteur
# de longueur N au lieu d'un seul scalaire).
# Necessaire pour le bloc "Outer product estimator" du cours : il faut le score
# individuel s_i = d ln f(y_i;theta)/d theta de CHAQUE observation, pas leur
# somme, pour pouvoir construire s_i * s_i' terme a terme.
# ------------------------------------------------------------------------------
loglik_i_QRE_CRRA <- function(params, dat_hetero){
  mu1 <- params[1]; mu2 <- params[2]
  r1  <- params[3]; r2  <- params[4]
  
  p_emp <- mean(dat_hetero$a1 == 1)
  q_emp <- mean(dat_hetero$a2 == 1)
  
  EU1_C <- q_emp * U(R, r1) + (1 - q_emp) * U(S, r1)
  EU1_D <- q_emp * U(T, r1) + (1 - q_emp) * U(P, r1)
  EU2_C <- p_emp * U(R, r2) + (1 - p_emp) * U(S, r2)
  EU2_D <- p_emp * U(T, r2) + (1 - p_emp) * U(P, r2)
  
  p1_C <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
  p2_C <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
  p1_C <- min(max(p1_C, 1e-10), 1 - 1e-10)
  p2_C <- min(max(p2_C, 1e-10), 1 - 1e-10)
  
  ll1 <- ifelse(dat_hetero$a1 == 1, log(p1_C), log(1 - p1_C))
  ll2 <- ifelse(dat_hetero$a2 == 1, log(p2_C), log(1 - p2_C))
  ll1 + ll2   # vecteur de longueur N : ln f(y_i;theta), i = 1..N (PAS de sum() ici)
}

# ------------------------------------------------------------------------------
# compute_variance_estimators : implemente EXACTEMENT les 4 lignes du cours
# (Hansen, section 10.15, Sample Hessian estimator + Outer product estimator),
# avec la notation E_N[.] = (1/N) * sum_i(.).
#
#   H_N(theta_hat)   = - E_N[ d2/dtheta dtheta' ln f(y_i;theta_hat) ]
#   V1               = [H_N(theta_hat)]^{-1}
#
#   I_theta(theta_hat) = E_N[ (d/dtheta ln f) (d/dtheta' ln f) ]
#   V2                 = [I_theta(theta_hat)]^{-1}
#
# NB : optim(hessian = TRUE) renvoie la Hessienne de la fonction MINIMISEE
# (-ell_N), SOMMEE sur les N observations (pas la moyenne) :
#   opt$hessian = - sum_i d2 ln f(y_i;theta_hat)/dtheta dtheta'
# Il faut donc diviser par N pour obtenir H_N(theta_hat) tel que defini
# dans le cours (E_N = moyenne, pas somme) -- c'est ce que fait la ligne
# "H_N <- opt$hessian / N" ci-dessous.
# Meme logique pour I_theta : t(scores) %*% scores = sum_i s_i s_i', donc on
# divise aussi par N.
#
# Enfin, V1 et V2 estiment la variance asymptotique de sqrt(N)*(theta_hat -
# theta_0) (resultat theorique standard du MLE), donc pour obtenir la
# covariance de theta_hat lui-meme (celle qui donne les vrais ecarts-types),
# il faut diviser une SECONDE fois par N : Cov(theta_hat) = V_j / N.
#
#
#     III - 0  Fonctions generiques (utilisees par les DEUX modeles)
# ------------------------------------------------------------------------------
compute_variance_estimators <- function(opt, ll_i_fn, dat, par_names = NULL){
  
  theta_hat <- opt$par
  if(is.null(par_names)) par_names <- names(theta_hat)
  N <- nrow(dat)
  
  ## ---- Sample Hessian estimator (cf. cours) ----
  H_N <- opt$hessian / N            # -E_N[ d2 ln f / dtheta dtheta' ]  (formule du cours, ligne 1)
  V1  <- solve(H_N)                 # V1 = [H_N]^{-1}                   (formule du cours, ligne 2)
  
  ## ---- Outer product estimator / BHHH (cf. cours) ----
  scores  <- numDeriv::jacobian(func = function(par) ll_i_fn(par, dat), x = theta_hat)  # matrice N x k : ligne i = s_i'
  I_theta <- (t(scores) %*% scores) / N   # = E_N[ s_i s_i' ]           (formule du cours, ligne 3)
  V2      <- solve(I_theta)               # V2 = [I_theta]^{-1}         (formule du cours, ligne 4)
  
  dimnames(V1) <- list(par_names, par_names)
  dimnames(V2) <- list(par_names, par_names)
  
  ## ---- Covariance de theta_hat (etape supplementaire, non ecrite dans le
  ## screen mais necessaire pour passer de Var(sqrt(N)(theta_hat-theta_0)) a
  ## Var(theta_hat) : meme logique que Var(X_bar) = Var(X)/N pour une moyenne ----
  Cov1 <- V1 / N
  Cov2 <- V2 / N
  dimnames(Cov1) <- list(par_names, par_names)
  dimnames(Cov2) <- list(par_names, par_names)
  
  se1 <- sqrt(diag(Cov1)); names(se1) <- par_names   # ecarts-types de theta_hat (Hessienne)
  se2 <- sqrt(diag(Cov2)); names(se2) <- par_names   # ecarts-types de theta_hat (OPG)
  
  list(N = N, H_N = H_N, I_theta = I_theta,
       V1 = V1, V2 = V2, Cov1 = Cov1, Cov2 = Cov2,
       se1 = se1, se2 = se2)
}

#-------------------------------------------------------------------------------
#                         III - 1  Estimation MLE
#-------------------------------------------------------------------------------

# optim minimise -ell_N(theta) = -sum_i ln f(y_i;theta) -> theta_hat = argmax ell_N
# hessian = TRUE demande a R de calculer numeriquement, EN theta_hat, la derivee
# seconde de la fonction minimisee -> c'est opt1$hessian, utilise plus bas.
opt1 <- optim(
  par = c(mu1 = 0.5, mu2 = 0.5, r1 = 0.3, r2 = 0.3),
  fn = loglik_QRE_CRRA,
  dat_hetero = dat_hetero,
  method = "L-BFGS-B",
  lower = c(0.1, 0.1, 0.1, 0.1),
  upper = c(2, 2, 1, 1),
  hessian = TRUE
)

opt1

#-------------------------------------------------------------------------------
#                III - 2  Estimation de la variance (apres le MLE)
#-------------------------------------------------------------------------------

# Appel de la fonction generique definie plus haut, sur le resultat de l'optim
# et sur la log-vraisemblance individuelle (necessaire pour les scores s_i).
var1 <- compute_variance_estimators(opt1, loglik_i_QRE_CRRA, dat_hetero,
                                    par_names = c("mu1","mu2","r1","r2"))

cat("N =", var1$N, "\n")   # rappel de la taille d'echantillon utilisee

# V1, V2 : les objets du cours (variance asymptotique de sqrt(N)(theta_hat-theta_0))
print(round(var1$V1, 5))
print(round(var1$V2, 5))

# se_V1, se_V2 : les ecarts-types REELLEMENT utilisables pour theta_hat
# (deja divises par N via Cov1, Cov2 dans la fonction ci-dessus)
print(round(rbind(se_V1 = var1$se1, se_V2 = var1$se2), 5))

# cov2cor s'appelle APRES que var1 existe, et sur Cov1 (la vraie covariance de
# theta_hat, PAS V1 qui est N fois plus grande -- mais la correlation
# implicite est la meme dans les deux cas car cov2cor est invariante d'echelle)
print(round(cov2cor(var1$Cov1), 5))


#-------------------------------------------------------------------------------
#                III  # Modèle contraint : r1 = r2 = 0
#-------------------------------------------------------------------------------


loglik_QRE_CRRA_constrained <- function(params, dat_hetero){
  
  mu1 <- params[1]
  mu2 <- params[2]
  r1  <- 0  # fixé
  r2  <- 0  # fixé
  
  # Sécurités
  if(any(params <= 0) || any(is.na(params)) || any(is.infinite(params))){
    return(1e10)
  }
  
  p_emp <- mean(dat_hetero$a1 == 1)
  q_emp <- mean(dat_hetero$a2 == 1)
  
  loglik <- 0
  
  for(k in 1:nrow(dat_hetero)){
    
    a1 <- dat_hetero$a1[k]; a2 <- dat_hetero$a2[k]
    
    # Joueur 1
    EU1_C <- q_emp * U(R, r1) + (1 - q_emp) * U(S, r1)
    EU1_D <- q_emp * U(T, r1) + (1 - q_emp) * U(P, r1)
    if(any(is.na(c(EU1_C, EU1_D))) || any(is.infinite(c(EU1_C, EU1_D)))){
      return(1e10)
    }
    p1_C <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
    p1_C <- max(min(p1_C, 1 - 1e-10), 1e-10)
    
    # Joueur 2
    EU2_C <- p_emp * U(R, r2) + (1 - p_emp) * U(S, r2)
    EU2_D <- p_emp * U(T, r2) + (1 - p_emp) * U(P, r2)
    if(any(is.na(c(EU2_C, EU2_D))) || any(is.infinite(c(EU2_C, EU2_D)))){
      return(1e10)
    }
    p2_C <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
    p2_C <- max(min(p2_C, 1 - 1e-10), 1e-10)
    
    # Log-vraisemblance
    if(a1 == 1){
      loglik <- loglik + log(max(p1_C, 1e-15))
    } else {
      loglik <- loglik + log(max(1 - p1_C, 1e-15))
    }
    if(a2 == 1){
      loglik <- loglik + log(max(p2_C, 1e-15))
    } else {
      loglik <- loglik + log(max(1 - p2_C, 1e-15))
    }
    
    if(is.infinite(loglik) || is.na(loglik)){
      return(1e10)
    }
  }
  
  return(-loglik)
}

# Estimation du modèle contraint
optim(
  par = c(mu1 = 0.5, mu2 = 0.5),
  fn = loglik_QRE_CRRA_constrained,
  dat_hetero = dat_hetero,
  method = "L-BFGS-B",
  lower = c(0.1, 0.1),
  upper = c(2, 2)
)

#-------------------------------------------------------------------------------
#                         III Estimation - MLE sans fixer q 
#-------------------------------------------------------------------------------


#---- Log-vraisemblance TOTALE avec mise a jour sequentielle des croyances ----
loglik_QRE_CRRA_learning <- function(params, dat_hetero){
  
  mu1 <- params[1]
  mu2 <- params[2]
  r1  <- params[3]
  r2  <- params[4]
  
  # Sécurités
  if(any(params <= 0) || any(is.na(params)) || any(is.infinite(params))){
    return(1e10)
  }
  if(abs(r1 - 1) < 1e-6 || abs(r2 - 1) < 1e-6){
    return(1e10)
  }
  
  # Initialisation des croyances
  p <- 0.5
  q <- 0.5
  alpha <- 0.9  # ← MODIFICATION : apprentissage plus rapide (0.5 au lieu de 0.9)
  
  loglik <- 0
  
  for(t in 1:nrow(dat_hetero)){
    
    a1 <- dat_hetero$a1[t]
    a2 <- dat_hetero$a2[t]
    
    # ---- Joueur 1 ----
    EU1_C <- q * U(R, r1) + (1 - q) * U(S, r1)
    EU1_D <- q * U(T, r1) + (1 - q) * U(P, r1)
    
    if(any(is.na(c(EU1_C, EU1_D))) || any(is.infinite(c(EU1_C, EU1_D)))){
      return(1e10)
    }
    
    p1_C <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
    p1_C <- max(min(p1_C, 1 - 1e-10), 1e-10)
    
    # ---- Joueur 2 ----
    EU2_C <- p * U(R, r2) + (1 - p) * U(S, r2)
    EU2_D <- p * U(T, r2) + (1 - p) * U(P, r2)
    
    if(any(is.na(c(EU2_C, EU2_D))) || any(is.infinite(c(EU2_C, EU2_D)))){
      return(1e10)
    }
    
    p2_C <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
    p2_C <- max(min(p2_C, 1 - 1e-10), 1e-10)
    
    # Log-vraisemblance
    if(a1 == 1){
      loglik <- loglik + log(max(p1_C, 1e-15))
    } else {
      loglik <- loglik + log(max(1 - p1_C, 1e-15))
    }
    
    if(a2 == 1){
      loglik <- loglik + log(max(p2_C, 1e-15))
    } else {
      loglik <- loglik + log(max(1 - p2_C, 1e-15))
    }
    
    if(is.infinite(loglik) || is.na(loglik)){
      return(1e10)
    }
    
    # Mise à jour des croyances
    p <- alpha * p + (1 - alpha) * a1
    q <- alpha * q + (1 - alpha) * a2
  }
  
  return(-loglik)
}

# ---- Log-vraisemblance INDIVIDUELLE (une contribution par periode) ----
loglik_i_QRE_CRRA_learning <- function(params, dat_hetero){
  mu1 <- params[1]; mu2 <- params[2]
  r1  <- params[3]; r2  <- params[4]
  
  n <- nrow(dat_hetero)
  p <- 0.5; q <- 0.5; alpha <- 0.9
  ll_t <- numeric(n)
  
  for(t in 1:n){
    a1 <- dat_hetero$a1[t]; a2 <- dat_hetero$a2[t]
    
    EU1_C <- q * U(R, r1) + (1 - q) * U(S, r1)
    EU1_D <- q * U(T, r1) + (1 - q) * U(P, r1)
    EU2_C <- p * U(R, r2) + (1 - p) * U(S, r2)
    EU2_D <- p * U(T, r2) + (1 - p) * U(P, r2)
    
    p1_C <- exp(EU1_C / mu1) / (exp(EU1_C / mu1) + exp(EU1_D / mu1))
    p2_C <- exp(EU2_C / mu2) / (exp(EU2_C / mu2) + exp(EU2_D / mu2))
    p1_C <- min(max(p1_C, 1e-10), 1 - 1e-10)
    p2_C <- min(max(p2_C, 1e-10), 1 - 1e-10)
    
    ll1_t <- if(a1 == 1) log(p1_C) else log(1 - p1_C)
    ll2_t <- if(a2 == 1) log(p2_C) else log(1 - p2_C)
    ll_t[t] <- ll1_t + ll2_t
    
    p <- alpha * p + (1 - alpha) * a1
    q <- alpha * q + (1 - alpha) * a2
  }
  ll_t
}


# Estimation avec meilleures bornes et initialisation
opt <- optim(
  par = c(mu1 = 0.7, mu2 = 0.7, r1 = 0.3, r2 = 0.3),
  fn = loglik_QRE_CRRA_learning, dat_hetero = dat_hetero,
  method = "L-BFGS-B",
  lower = c(0.1, 0.1, 0.01, 0.01),
  upper = c(3, 3, 0.99, 0.99),      # mêmes bornes pour les deux joueurs
  hessian = TRUE, control = list(maxit = 5000)
)

print(opt$par)
print(opt$value)


#----------------------------------------------------------
# Profil de vraisemblance de $\mu_2$ (autres paramètres ré-optimisés à chaque point, multi-start)
#----------------------------------------------------------
mu2_grid <- seq(0.3, 6, by = 0.3)
profile_ll_true <- numeric(length(mu2_grid))

# point de depart initial
par_start <- c(mu1 = 0.6, r1 = 0.25, r2 = 0.6)

for(i in seq_along(mu2_grid)){
  
  m2_fixed <- mu2_grid[i]
  
  fn_conditional <- function(par_sub){
    full_par <- c(par_sub[1], m2_fixed, par_sub[2], par_sub[3])
    loglik_QRE_CRRA_learning(full_par, dat_hetero)
  }
  
  # plusieurs points de depart pour eviter les minima locaux
  starts <- list(
    par_start,                          # warm start (solution precedente)
    c(mu1 = 0.6, r1 = 0.25, r2 = 0.6),  # point de depart "neutre"
    c(mu1 = 0.6, r1 = 0.25, r2 = 2.5)   # point de depart avec r2 eleve
  )
  
  best <- NULL
  for(s in starts){
    o <- tryCatch(
      optim(par = s, fn = fn_conditional, method = "L-BFGS-B",
            lower = c(0.1, 0.01, 0.01), upper = c(2, 0.99, 5)),
      error = function(e) NULL
    )
    if(!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  
  profile_ll_true[i] <- best$value
  par_start <- best$par   # warm start pour le point suivant
}

plot(mu2_grid, -profile_ll_true, type = "l",
     xlab = "mu2 (fixe)", ylab = "log-vraisemblance profilee",
     main = "Vrai profil de vraisemblance en mu2 (multi-start)")

#----------------------------------------------------------
# Estimation de la variance (memes conventions que le modele statique :
# H_N, I_theta, V1, V2 au sens du cours, puis Cov(theta_hat) = V_j / N)
#----------------------------------------------------------

var2 <- compute_variance_estimators(opt, loglik_i_QRE_CRRA_learning, dat_hetero,
                                    par_names = c("mu1","mu2","r1","r2"))

cat("N =", var2$N, "\n")

cat("\n--- H_N(theta_hat) ---\n")
print(round(var2$H_N, 5))
cat("\n--- I_theta(theta_hat) (OPG) ---\n")
print(round(var2$I_theta, 5))

cat("\n--- V1 = [H_N]^-1 (au sens du cours) ---\n")
print(round(var2$V1, 5))
cat("\n--- V2 = [I_theta]^-1 (au sens du cours) ---\n")
print(round(var2$V2, 5))

cat("\n--- Cov(theta_hat) via Hessienne (= V1 / N) ---\n")
print(round(var2$Cov1, 5))
cat("\n--- Cov(theta_hat) via OPG (= V2 / N) ---\n")
print(round(var2$Cov2, 5))

cat("\n--- Ecarts-types de theta_hat ---\n")
print(round(rbind(se_hessian = var2$se1, se_opg = var2$se2), 5))

cat("\n--- Corrélation implicite (via Cov1) ---\n")
print(round(cov2cor(var2$Cov1), 5))


#-------------------------------------------------------------------------------
#         III Estimation - Equilibre de Nash perturbé (best response)
#                   avec hétérogénéité (r1 et r2)
#-------------------------------------------------------------------------------

# =========================
# Best response déterministe (argmax) - avec r spécifique
# =========================
best_response <- function(q, r){
  
  EU_C <- q * U(R, r) + (1 - q) * U(S, r)
  EU_D <- q * U(T, r) + (1 - q) * U(P, r)
  
  if(EU_C > EU_D){
    return(1)   # Coopérer
  } else {
    return(0)   # Dévier
  }
}

# =========================
# Perturbation des croyances (ε2)
# =========================
beliefs <- function(p, q, sigma_b){
  set.seed(123)
  
  p_pert <- p + rnorm(1, 0, sigma_b)
  q_pert <- q + rnorm(1, 0, sigma_b)
  
  p_pert <- min(max(p_pert, 0), 1)
  q_pert <- min(max(q_pert, 0), 1)
  
  return(list(p = p_pert, q = q_pert))
}

# =========================
# Log-vraisemblance (best-response perturbé) - avec r1 et r2
# =========================
loglik_br_perturbed <- function(params, dat_hetero){
  
  r1 <- params[1]
  r2 <- params[2]
  sigma_a <- params[3]
  sigma_b <- params[4]
  
  ll <- 0
  
  for(i in 1:nrow(dat_hetero)){
    
    # actions observées
    a1 <- dat_hetero$a1[i]
    a2 <- dat_hetero$a2[i]
    
    # croyances empiriques (historique)
    if(i == 1){
      p0 <- 0.5
      q0 <- 0.5
    } else {
      p0 <- mean(dat_hetero$a1[1:(i-1)] == 1, na.rm = TRUE)
      q0 <- mean(dat_hetero$a2[1:(i-1)] == 1, na.rm = TRUE)
    }
    
    # perturbation des croyances
    b <- beliefs(p0, q0, sigma_b)
    p <- b$p
    q <- b$q
    
    # best responses (r1 pour joueur 1, r2 pour joueur 2)
    br1 <- best_response(q, r1)  # joueur 1 utilise r1
    br2 <- best_response(p, r2)  # joueur 2 utilise r2
    
    # erreurs d'exécution
    e1 <- a1 - br1
    e2 <- a2 - br2
    
    # vraisemblance gaussienne
    ll <- ll - (e1^2 + e2^2) / (2 * sigma_a^2) - log(sigma_a)
  }
  
  return(-ll)
}

# =========================
# Estimation
# =========================
set.seed(123)

# Estimation avec contrainte r2 > r1
optim(
  par = c(r1 = 0.3, r2 = 0.6, sigma_a = 0.5, sigma_b = 0.1),
  fn = loglik_br_perturbed,
  dat_hetero = dat_hetero,
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, 0.01, 0.01),
  upper = c(2, 2, 2, 1),
  control = list(maxit = 1000)
)