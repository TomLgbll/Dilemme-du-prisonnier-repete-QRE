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

cat("\n==============================\n")
cat("        RESULTATS QRE\n")
cat("==============================\n\n")

print(
  transform(
    results,
    p = round(p, 4),
    q = round(q, 4)
  ),
  row.names = FALSE
)

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

cat("\n==============================\n")
cat("     PARAMETRES PAIRE 1\n")
cat("==============================\n\n")
cat("Joueur 1 : r1 =", round(r1_pair1, 3), ", mu1 =", round(mu1_pair1, 3), "\n")
cat("Joueur 2 : r2 =", round(r2_pair1, 3), ", mu2 =", round(mu2_pair1, 3), "\n")


# Probabilités observées (sans arrondi)
p_obs <- mean(pair1$action_joueur1 == "C")
q_obs <- mean(pair1$action_joueur2 == "C")

cat("\n------------------------------\n")
cat(" PROBABILITES DE COOPERATION\n")
cat("------------------------------\n")
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

cat("\n==============================\n")
cat("   EFFET DE L'AVERSION AU RISQUE\n")
cat("==============================\n\n")
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
