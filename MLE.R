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

cat("\n==============================\n")
cat("        RESULTAT MLE STATIQUE\n")
cat("==============================\n\n")
opt1

#-------------------------------------------------------------------------------
#                III - 2  Estimation de la variance (apres le MLE)
#-------------------------------------------------------------------------------

# Appel de la fonction generique definie plus haut, sur le resultat de l'optim
# et sur la log-vraisemblance individuelle (necessaire pour les scores s_i).
var1 <- compute_variance_estimators(opt1, loglik_i_QRE_CRRA, dat_hetero,
                                    par_names = c("mu1","mu2","r1","r2"))
cat("\n==============================\n")
cat("      VARIANCE DES ESTIMATEURS\n")
cat("==============================\n\n")


cat("N =", var1$N, "\n")   # rappel de la taille d'echantillon utilisee

# V1, V2 : les objets du cours (variance asymptotique de sqrt(N)(theta_hat-theta_0))
cat("--- V1 : Hessienne ---\n")
print(round(var1$V1, 5))

cat("\n--- V2 : OPG ---\n")
print(round(var1$V2, 5))

# se_V1, se_V2 : les ecarts-types REELLEMENT utilisables pour theta_hat
# (deja divises par N via Cov1, Cov2 dans la fonction ci-dessus)
cat("\n--- Ecarts-types ---\n")
print(round(rbind(se_V1 = var1$se1, se_V2 = var1$se2), 5))

# cov2cor s'appelle APRES que var1 existe, et sur Cov1 (la vraie covariance de
# theta_hat, PAS V1 qui est N fois plus grande -- mais la correlation
# implicite est la meme dans les deux cas car cov2cor est invariante d'echelle)
cat("\n--- Matrice de corrélation ---\n")
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
cat("\n==============================\n")
cat("      RESULTAT MLE LEARNING\n")
cat("==============================\n\n")

cat("Paramètres estimés :\n")
print(opt$par)

cat("\nValeur de la log-vraisemblance :\n")
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

cat("\n==============================\n")
cat(" VARIANCE MODELE LEARNING\n")
cat("==============================\n\n")

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

