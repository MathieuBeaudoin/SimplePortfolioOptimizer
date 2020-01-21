##################################################################################################
#####                                 IMPORTANT DISCLAIMER                                   #####

###       THIS CODE IS NOT INTENDED TO ASSIST IN ANY ACTUAL, REAL-WORLD INVESTMENT DECISIONS.  ###
### I WILL NOT BE HELD RESPONSIBLE FOR ANY LOSSES FLOWING FROM THE USE OR MISUSE OF THIS CODE. ###
##################################################################################################

rm(list=ls())

### HYPERPARAMETRES ###
ITERATIONS <- 1000  # Nombre de gradient descents a effectuer
MAX_ITER <- 100     # Maximum d'iterations pour une gradient descent
INIT_DW <- 1        # Increment/decrement initial = 0.1/exp(INIT_DW)
MAX_DW_LEVELS <- 5  # Plus petit increment/decrement = 0.1/exp(MAX_DW_LEVELS)
MAX_PRECISION <- 0.0001
PENALITE_CONTRAINTE <- 50

### PARAMETRES D'IMPORTATION DE DONNEES ###
# Mode pour calcul de MU, SIGMA et RHO
input_mode <- 2   # 1 pour entree manuelle, 2 pour calcul a partir de donnees brutes d'un CSV
# Parametres pour mode 2
source <- "data.csv"
orientation <- 2  # 1 si chaque colonne est un actif, 2 si chaque rangee est un actif
separateur <- ";" # Separateur utilise par le fichier CSV

### PARAMETRES D'OPTIMISATION ###
target_MU <- 0.10
target_SIGMA <- 0.20
contrainte_poids_neg <- T   # Contrainte: F si on permet des poids negatifs dans le portefeuille
score_type <- 2 # 1 = minimiser difference quadratique
                # 2 = minimiser variance pour rendement donne
                # 3 = maximiser rendement pour variance donnee

### DONNEES ###
getData <- function(mode=1, src="", orient=1, separ=";") {
  # Mode 1: mu, sigma et rho "hard-coded"
  # Mode 2: données brutes importées
    # Orientations (mode 2): 1 -> chaque colonne=1 actif; 2 -> chaque rangee=1 actif
  if(mode==1){
    MU <- c(0.06, 0.1, 0.03)
    SIGMA <- c(0.1, 0.2, 0)
    RHO <- c(   1, -0.2,    0, 
             -0.2,    1,    0, 
                0,    0,    1)
    N <- length(MU)
    RHO <- matrix(RHO, byrow=T, nrow=sqrt(length(RHO)))
  }
  else {
    raw <- as.matrix(read.csv(src, sep=separ, header=F))
    if(orientation==2) raw <- t(raw)
    N <- dim(raw)[2]
    MU <- SIGMA <- rep(0,N)
    RHO <- matrix(rep(0,N*N), nrow=N)
    for(i in 1:N) {
      MU[i] <- mean(raw[,i])
      SIGMA[i] <- sqrt(mean(raw[,i]**2) - MU[i]**2)
    }
    for(i in 1:N){
      for(j in 1:N){
        RHO[i,j] <- (mean(raw[,i] * raw[, j]) - MU[i]*MU[j]) / (SIGMA[i]*SIGMA[j])
      }
    }
  }
  return(list(N, MU, SIGMA, RHO))
}
DATA <- getData(mode=input_mode, src=source, orient=orientation, separ=separateur)

# Morcellement
N <- DATA[[1]]
MU <- DATA[[2]]
SIGMA <- DATA[[3]]
RHO <- DATA[[4]]

# Calcule le rendement espere du portefeuille
ptf_mu <- function(W) {
  return(sum(MU*W))
}

# Calcule la volatilite du portefeuille
ptf_sigma <- function(W) {
  variance <- 0
  for(i in 1:N) {
    for(j in 1:N) {
      variance <- variance + RHO[i,j] * prod(W[c(i,j)]) * prod(SIGMA[c(i,j)])
    }
  }
  return(sqrt(variance))
}

# Initialisation aleatoire, pour eviter de se prendre dans un minimum local
init_weight <- function() {
  amplitudes <- rexp(N)
  centre <- mean(amplitudes)
  if(contrainte_poids_neg) W <- amplitudes / sum(amplitudes)
  else  W <- (amplitudes-centre)/centre + 1/N
  return(W)
}

# Mesure a minimiser
score <- function(W) {
  metric <- 0
  if(score_type==1) {
    metric <- (ptf_mu(W)-target_MU)^2 + (ptf_sigma(W)-target_SIGMA)^2
  }
  else if(score_type==2) {
    if(ptf_mu(W)<target_MU) metric <- PENALITE_CONTRAINTE*(target_MU-ptf_mu(W))
    metric <- metric + ptf_sigma(W)
  }
  else if(score_type==3) {
    if(ptf_sigma(W)>target_SIGMA) metric <- PENALITE_CONTRAINTE*(ptf_sigma(W)-target_SIGMA)
    metric <- metric + 1/ptf_mu(W)
  }
  else metric <- PENALITE_CONTRAINTE
  return(metric)
}

# Maintient les contraintes sum(W)==1 et W>=0 (si applicable)
adjust_weights <- function(W,i,dw) {
  dW <- rep(0,N)
  dW[i] <- dw
  dW[-i] <- -dw * W[-i]/sum(W[-i])
  W <- W + dW
  if(contrainte_poids_neg) {
    if(is.na(sum(W<0)>0)) print(W)
    while(sum(W<0)>0) {
      negatifs <- which(W<0)
      autres <- W[-c(i, negatifs)]
      a_redistribuer <- sum(W[negatifs])
      W[negatifs] <- 0
      autres <- autres + a_redistribuer * autres/sum(autres)
      W[-c(i, negatifs)] <- autres
    }
  }
  return(W)
}

# Approximation numerique des derivees partielles
derivee_part <- function(W,i,dw) {
  if(contrainte_poids_neg){
    if(dw<0) dw <- min(-dw, W[i])
    else dw <- min(dw, 1-W[i])
  }
  return(score(adjust_weights(W,i,dw)) - score(W))
}

# Recherche d'un portefeuille minimisant la distance quadratique
# entre les caracteristiques du portefeuille et les objectifs
gradient_descent <- function() {
  W <- init_weight()
  iter <- 0
  dw_level <- INIT_DW
  learning <- as.data.frame(matrix(c(W, 
                                     0.1/exp(dw_level), 
                                     ptf_mu(W), 
                                     ptf_sigma(W),
                                     score(W)), byrow=T, nrow=1))
  while((abs(ptf_mu(W)-target_MU) > MAX_PRECISION |
         abs(ptf_sigma(W)-target_SIGMA) > MAX_PRECISION) &
        dw_level <= MAX_DW_LEVELS &
        iter <= MAX_ITER) {
    dw <- 0.1/exp(dw_level)
    iter <- iter+1
    derivees <- c()
    for(i in 1:(2*N)) {
      if(i<=N) derivee <- derivee_part(W,i,dw)
      else derivee <- derivee_part(W,i-N,-dw)
      derivees <- c(derivees, derivee)
    }
    meilleur <- which(derivees==min(derivees))
    if(length(meilleur>1)) meilleur <- meilleur[1] # En cas d'egalite
    if(min(derivees)>=0) dw_level <- dw_level+1
    else {
      if(meilleur>N) {
        dw <- -1*dw
        meilleur <- meilleur-N
      }
      W <- adjust_weights(W, meilleur, dw)
      this_iter <- matrix(c(W, dw, ptf_mu(W), ptf_sigma(W), score(W)), 
                          byrow=T, nrow=1)
      learning <- rbind(learning, this_iter)
    }
  }
  colnames(learning) <- c(paste0("w",1:N), "Increment", "Mu", "Sigma", "Score")
  return(list(W, learning, iter))
}

# Repete la "gradient descent" pour valider que les solutions convergent
main_loop <- function() {
  resultats <- as.data.frame(matrix(rep(0, N+4), byrow=T, nrow=1))
  learning_details <- list()
  for(i in 1:ITERATIONS) {
    one_run <- gradient_descent()
    W <- one_run[[1]]
    mu <- ptf_mu(W)
    sigma <- ptf_sigma(W)
    cycles <- one_run[[3]]
    print(paste0("Iteration #", i, " Mu=", mu, " Sigma=", sigma,
                 " Cycles G/D=", cycles))
    this_result <- matrix(c(W, mu, sigma, score(W), cycles), 
                          byrow=T, nrow=1)
    resultats <- rbind(resultats, this_result)
    learning_details[[i]] <- one_run[[2]]
  }
  colnames(resultats) <- c(paste0("w",1:N), "Mu", "Sigma", "Score", "Cycles")
  return(list(resultats[-1,], learning_details))
}

# Execution de l'analyse et presentation des resultats
loop <- main_loop()
analyse <- loop[[1]]
output <- t(analyse[which(analyse$Score==min(analyse$Score)),])
contrainte_poids_neg <- F
loop <- main_loop()
analyse <- loop[[1]]
output <- cbind(output, t(analyse[which(analyse$Score==min(analyse$Score)),]))
output <- output[1:(N+2),]
rownames(output) <- c(paste("Actif", 1:N), "Rendement", "Ecart-type")
colnames(output) <- c("Avec contrainte", "Sans contrainte")
write.csv(output, "ptf.csv")
