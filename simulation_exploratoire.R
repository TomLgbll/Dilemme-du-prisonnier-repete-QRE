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

cat("\n==============================\n")
cat("   FREQUENCE DE COOPERATION\n")
cat("==============================\n\n")

cat("Joueur 1 :", round(prop_C_j1 * 100, 2), "% de coopération\n")
cat("Joueur 2 :", round(prop_C_j2 * 100, 2), "% de coopération\n\n")


cat("------------------------------\n")
cat("Coopération moyenne par paire\n")
cat("------------------------------\n")

aggregate(
  cbind(
    mean_j1 = as.numeric(data$action_joueur1 == "C"),
    mean_j2 = as.numeric(data$action_joueur2 == "C")
  ),
  by = list(pair = data$pair),
  FUN = mean
)
