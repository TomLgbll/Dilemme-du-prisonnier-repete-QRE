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
opt_br <- optim(
  par = c(r1 = 0.3, r2 = 0.6, sigma_a = 0.5, sigma_b = 0.1),
  fn = loglik_br_perturbed,
  dat_hetero = dat_hetero,
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, 0.01, 0.01),
  upper = c(2, 2, 2, 1),
  control = list(maxit = 1000)
)

cat("\n==============================\n")
cat("   RESULTAT BEST RESPONSE\n")
cat("==============================\n\n")

cat("Paramètres estimés :\n")
print(round(opt_br$par, 4))

cat("\nValeur de la log-vraisemblance :",
    round(-opt_br$value, 4), "\n")