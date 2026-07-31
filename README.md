Projet réalisé dans le cadre d'un **stage de Master 1 Statistique**.

## Objectif

L'objectif de ce stage est d'appliquer les outils de l'économétrie des interactions stratégiques et de la théorie des jeux au dilemme du prisonnier répété.

Le projet consiste à :
- Modéliser le comportement d'agents en interaction répétée.
- Simuler des agents hétérogènes caractérisés par des paramètres comportementaux.
- Estimer ces paramètres à partir des données simulées via des méthodes économétriques.

## Méthodologie

Le modèle est progressivement enrichi :
- QRE standard avec un paramètre de rationalité/bruit $\mu$.
- Introduction d'une hétérogénéité entre joueurs.
- Ajout d'une fonction d'utilité CRRA permettant de modéliser l'aversion au risque.
- Estimation du modèle complet avec les paramètres $(\mu_1,\mu_2,r_1,r_2)$.

## Structure du projet

Le dépôt contient deux versions du code :

- **`Stage_code.R`** : version complète regroupant l'ensemble du projet dans un seul script. Elle permet d'exécuter toutes les étapes de manière séquentielle.
- **Version modulaire** : code séparé en plusieurs fichiers afin de faciliter la lecture, la maintenance et la réutilisation des fonctions.

L'organisation est la suivante :


```text
Dilemme-du-prisonnier-repete-QRE
│
├── Stage_final_juillet.R        # Script complet du projet
│
└── Version modulaire
    │
    ├── main.R                   # Script principal
    │
    ├── main_construction.R      # Construction des données
    │   │
    │   ├── simulation_exploratoire.R  # Simulation des agents
    │   └── QRE_modeles.R              # Modèles QRE
    │
    └── main_estimation.R        # Estimation des modèles
        │
        ├── Monte_Carlo.R        # Expériences Monte Carlo
        ├── MLE.R                # Maximum de Vraisemblance
        └── best_response.R      # Modèle best response perturbé
```
