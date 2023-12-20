
#---------------- Simulations Eclatement clusters ----------------

rKbarre_to_K <-function(params, size){
  #---------------------- Variables ----------------------------
  # idx_K : Indices des colonnes pour les clusters de chaque composante, afin de permettre le choix de quel cluster final influencera quels clusters au sein de chaque composantes
  # size : permet de choisir le nombre de tirage fait au sein des groupes pour composer les groupes éclatés.
  # clusters : Lien entre K_barre et K
  
  Q = params$Q
  K_barre = params$K_barre
  K = params$K
  
  #---------------------- Indice composantes ----------------------------
  idx_K <- list()
  for(q in 1:Q){
    idx_K[[q]] <- 1:K[q]
  }
  
  #---------------------- Eclatement des clusters ----------------------------
  clusters <- list()
  for(k in 1:(K_barre-1)){
    tmp2 <- list()
    for(q in 1:Q){
      if(identical(idx_K[[q]], integer(0)) ){next;} #vérification vecter indice pas vide
      
      #------ Partie du code modifiable pour ajuster les simulations ---------------
      
      ## Simulation qui fait le plus de sens :
      # si K multiple de K_barre : permet d'avoir un nombre de clusters qui font de bons groupes
      idx <- unique(sample(idx_K[[q]],size = size,replace = F))
      
      ## Simulation avec + ou - de sens dans la création des clusters :
      #Choisi un nombre size de clusters à chaque fois, avec possibilité d'en tirer les mêmes,
      #idx <- unique(sample(idx_K[[q]],size = size,replace = T))
      
      
      
      ## Simulation chaotique :
      # Nombre de clusters tirés aléatoirement au lieu de size
      #idx <- unique(sample(idx_K[[q]],size = sample(1:K[q],1),replace = T))
      
      #------------------------------------------------------------------
      idx_K[[q]] <- idx_K[[q]][-which(idx_K[[q]] %in% idx)]
      tmp2[[q]] <- idx
    }
    clusters[[k]] <- tmp2
  }
  
  
  #-----------------------------  Dernier cluster -----------------------------------
  
  # On se sert des clusters non utilisés pour former le dernier groupe
  k = K_barre
  tmp2 <- list()
  for(q in 1:Q){
    if(identical(idx_K[[q]], integer(0)) ){next;}
    tmp2[[q]] <-  idx_K[[q]]
  }
  clusters[[k]] <- tmp2
  
  #---------------------- Output ----------------------------
  return(clusters)
}


rZ <- function(params,Z_barre,clusters,sorted_by_Z=F){
  #---------------------- Variables ----------------------------
  # N : Nombre d'observations à générer
  # Q : Nombre de composantes du mélange de vues
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # Z : Variables simulées d'appartenance à des clusters selon la composante du mélange de vue
  # Z_barre : vrai clustering des données.
  
  # sorted_by_Z : variable indiquant l'ordre de classement des variables
  # F -> dans l'ordre de Z_barre
  # T -> dans l'ordre des clusters générés dans chaque Z
  
  N <- params$N
  Q <- params$Q
  K <- params$K
  Z <- list()
  
  #---------------------- Génération des Z ----------------------------
  for(q in 1:Q){
    val <- matrix(0,N,K[q])
    for(i in 1:N){
      if (q > length(clusters[[Z_barre[i]]])){idx <- sample(1:K[q],1)}
      # Erreur cas où on a pas de groupe prédéfini pour le dernier élément des listes
      else if(is.null(clusters[[Z_barre[i]]][[q]])){idx <- sample(1:K[q],1)}
      # Erreur cas où on a pas de groupe prédéfini pour un élément (PAS DERNIER des listes).
      else if(length(clusters[[Z_barre[i]]][[q]]) ==1) {val[i,clusters[[Z_barre[i]]][[q]]] = 1  }
      # Erreur lorsqu'il y a qu'un seul cluster dans un le choix car sample considère ça comme un int
      # donc utilise la fonction sample.int au lieu de sample.
      else {idx <- sample(x=clusters[[Z_barre[i]]][[q]],size=1);val[i,idx] = 1  }
      
    }
    Z[[q]] <- max.col(val)
  }
  
  if(sorted_by_Z){
    for(q in Q:1){
      Z[[q]] <- Z[[q]][order(Z[[1]])]  #sort(max.col(Z2[[q]]))
    }
  }
  
  
  return(Z)
}


rmultimodal <- function(params,size=2,sorted_by_Z=F){
  #---------------------- Variables ----------------------------
  # params$N : Nombre d'observations à générer
  # params$V : Nombre de vues à générer
  # params$Q : Nombre de composantes du mélange de vues
  # params$K_barre : Nombre de vrais clusters
  # params$K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # clusters : Lien entre K_barre et K
  
  
  # params$p_Zbarre : vecteur de probabilité pour les clusters de Z_barre
  # params$p_W: vecteur de probabilité pour les clusters des vues
  
  # clusters : Lien entre K_barre et K
  # Z_barre : vrai clustering des données.
  # Z :Variables simulées d'appartenance à des clusters selon la composante du mélange de vue
  
  
  # size : Nombre d'éléments à mettre pour le lien entre K_barre et K
  # sorted_by_Z : variable indiquant l'ordre de classement des variables
  # F -> dans l'ordre de Z_barre
  # T -> dans l'ordre des clusters générés dans chaque Z
  # B : Matrices d'émissions
  # A : Matrices d'adjacences simulées
  
  #----------------------  K_Barre -> K ----------------------------
  
  clusters <- rKbarre_to_K(params, size=size)
  
  # ----------------------  Z_barre et Z ----------------------------
  Z_barre <- t(rmultinom(n=params$N,size =1,prob = params$p_Zbarre))
  Z_barre <- sort(max.col(Z_barre))
  
  Z <- rZ(params,Z_barre,clusters,sorted_by_Z=sorted_by_Z)
  
  # ----------------------  Vues et mélange de vues ----------------------------
  W <- t(rmultinom(n=params$V,size =1,prob = params$p_W ))
  W <- sort(max.col(W))
  
  
  # ----------------------  Matrices alpha ----------------------------
  
  # A mettre directement dans params ?
  
  B <- list()
  a = 0.99
  b=0.01
  for(q in 1:params$Q) {B[[q]] = matrix(b,params$K[q],params$K[q]); diag(B[[q]]) <- a}
  
  # ----------------------  Matrices d'adjacences ----------------------------
  A <- array(NA,dim=c(params$N,params$N,params$V))
  for(l in 1:params$V){
    for(i in 1:(params$N-1)){
      for(j in (i+1):params$N){
        q = W[l]
        Zi = Z[[q]][i]
        Zj =  Z[[q]][j]
        A[i,j,l] <- A[i,j,l] <- rbinom(1,1,prob=B[[q]][Zi,Zj])
        A[j,i,l] <- A[i,j,l]
      }
    }
    diag(A[,,l]) <- 0
  }
  
  # ----------------------  Output ----------------------------
  params$B <- B
  simulations <- list(A = A, Z = Z, Z_barre = Z_barre, clusters = clusters, W = W )
  output <- list(params=params,simulations=simulations)
  return(output)
}
