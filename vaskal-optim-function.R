######## Vasicek 


Vasicek_ML_Estimation <- function(ZK){
  n = length(ZK)
  delta = 1  # delta
  Sx = sum(ZK[1:n - 1])
  Sy = sum(ZK[2:n])
  Sxx = sum((ZK[1:n - 1])^2)
  Syy = sum((ZK[2:n])^2)
  Sxy = sum((ZK[1:n - 1]) * (ZK[2:n]))
  theta = (Sy * Sxx - Sx * Sxy)/((n - 1) * (Sxx - Sxy) - (Sx^2 - Sx * Sy))
  (Sy*Sxx - Sx*Sxy) / ( n*(Sxx - Sxy) - (Sx^2 - Sx*Sy))
  kappa = -log((Sxy - theta * Sx - theta * Sy + (n - 1) * theta^2)/(Sxx - 2 * theta * Sx + (n - 1) * theta^2))/delta
  a = exp(-kappa * delta)
  sigmah2 = (Syy - 2 * a * Sxy + a^2 * Sxx - 2 * theta * (1 - a) * (Sy - a * Sx) + (n - 1) * theta^2 * (1 - a)^2)/(n - 1)
  sigma = sqrt((sigmah2) * 2 * kappa/(1 - a^2))
  list = list(kappa = kappa, theta = theta, sigma = sigma)
  return(list)
}


Vasicek_Simulation_neu <- function(ZK,T=length(ZK),N){
  
  ## Vasicek Maximum Likelihood Estimation
  n <- length(ZK)
  delta <- 1  # delta
  Sx <- sum(ZK[1:n - 1])
  Sy <- sum(ZK[2:n])
  Sxx = sum((ZK[1:n - 1])^2)
  Syy = sum((ZK[2:n])^2)
  Sxy = sum((ZK[1:n - 1]) * (ZK[2:n]))
  theta = (Sy * Sxx - Sx * Sxy)/((n - 1) * (Sxx - Sxy) - (Sx^2 - Sx * Sy))
  (Sy*Sxx - Sx*Sxy) / ( n*(Sxx - Sxy) - (Sx^2 - Sx*Sy))
  kappa = -log((Sxy - theta * Sx - theta * Sy + (n - 1) * theta^2)/(Sxx - 2 * theta * Sx + (n - 1) * theta^2))/delta
  a = exp(-kappa * delta)
  sigmah2 = (Syy - 2 * a * Sxy + a^2 * Sxx - 2 * theta * (1 - a) * (Sy - a * Sx) + (n - 1) * theta^2 * (1 - a)^2)/(n - 1)
  sigma = sqrt((sigmah2) * 2 * kappa/(1 - a^2))
  
  
  ##### Vasicek Simulation 
  # R_0 = Startwert, 1 Wert Zinskurve ? 
  # T = Länge Zinskurve (YTM) , N = Simulationsanzahl 
  # ZK in Anleihenpreise umwandeln exp(-z(t,T)*(T-t))
  
  Vek_1 <- c(rep(99,T*N))
  R_t <- matrix(Vek_1,nrow = T,ncol=N)
  R_t[1,] <- c(rep(ZK[1],N))
  Vasicek_Bond <- matrix(Vek_1,nrow = T,ncol=N)
  Vasicek_Yield <-  matrix(Vek_1,nrow = T,ncol=N)
  for(k in 1:N){
    for(i in 2:T){ # tau = T-t also Restlaufzeit in k=tau für die for schleife
      R_t[i,k] <- R_t[i-1,k] + kappa*(theta-R_t[i-1,k])+sigma*rnorm(1) # Shortrate
      
    }
    for (i in 1:T) {
      Vasicek_Bond[i,k] <- exp(-R_t[i,k]*1/kappa*(1-exp(-kappa*(k)))+(theta-sigma^2/(2*kappa^2))*(1/kappa*(1-exp(-kappa*(k)))-k)-sigma^2/(4*kappa)*(1/kappa*(1-exp(-kappa*(k))))^2)
      Vasicek_Yield[i,k]<- -log(Vasicek_Bond[i,k]/(i)) /100
    }  
  }
  return(Vasicek_Yield)
}

Vasicek_ZK_Bewertung <- Vasicek_Simulation_neu(ZK=RFR_spot_no_VA,N=10)
head(tests)

plot(RFR_spot_no_VA,ylim=c(-0.005,0.05))
points(Vasicek_ZK_Bewertung[,1],col="blue")
points(Vasicek_ZK_Bewertung[,2],col="red")
points(Vasicek_ZK_Bewertung[,3],col="green")

hist(tests[,2]-tests[,1]) # Histogram Differenz Vasicek Zinskurven Schätzung



#############################################################################

# Kalibrierung KQ 
plot(RFR_spot_no_VA)

RFR_spot_no_VA
Bonds_RFR_spot_no_VA<-c()
for (i in 1:length(RFR_spot_no_VA)) {
  Bonds_RFR_spot_no_VA[i] <- exp(-RFR_spot_no_VA[i]*i)
}
plot(Bonds_RFR_spot_no_VA)  # Bondpreise EIOPA RFR
plot(RFR_spot_no_VA)

write.csv(RFR_spot_no_VA,"RFR.csv")


########################################################################################

### Vasicek Estimation with Kalman Filter
# auf Basis 

# Bundesbank November 2018.csv
y <- read.csv("C:/Users/Administrator/Downloads/Maik-Upwork/RFR.csv",1)
y <- SwapData

### Vasicek model and the Kalman filter

kappa <- c(0.5,0.005)
#c[3]
alpha <- -3 
#c[4]
H <- c(rep(-1,19))
#c[23]
psi <- 0.004

VasKal <- function(param){
  # Initialization
  # Random values -> should theoretically not matter
  
  F_0 <- c(0.05,0.04)
  P_0 <- matrix(c(0.2,0.2,0.2,0.2),2,2)
  tau <- c(1:19)
  tau_n <- length(tau)
  
  # Make Parameters non-negative
  H_vec <- exp(param[4:22])
  alpha_new <- exp(param[3])
  
  kappa <- c(param[1],param[2])
  psi <- param[23]
  
  #2 Thetas, result of theta* = 0, Vasicek model beta = 0  
  
  theta <- psi*alpha_new/kappa
  theta_1 <- psi*alpha_new/kappa[1]
  theta_2 <- psi*alpha_new/kappa[2]
  
  # H_m stands for H Matrix
  H_m <- diag(H_vec,tau_n,tau_n)
  
  # Formulas for A and B 
  A_factor_1 <- c(theta[1]-0.5*(alpha_new^2/kappa[1]^2))
  A_factor_2 <- c(theta[2]-0.5*(alpha_new^2/kappa[2]^2))
  
  B_tau_1 <- (1-exp(-c(kappa[1]*tau)))/kappa[1]
  B_tau_2 <- (1-exp(-c(kappa[2]*tau)))/kappa[2]
  A_tau_1 <- c(A_factor_1*c(B_tau_1-tau)-c((alpha_new^2*B_tau_1^2)/(4*kappa[1])))
  A_tau_2 <- c(A_factor_2*c(B_tau_2-tau)-c((alpha_new^2*B_tau_2^2)/(4*kappa[2])))
  A_tau <- c(A_tau_1+A_tau_2)
  
  # Creating Observation Equation parameters
  A <- (as.matrix((-A_tau)/tau))
  B_tau <- c(c(B_tau_1/tau),c(B_tau_2/tau))
  #tau * factor
  B <- matrix(B_tau,tau_n,2)
  diag_phi <- c(exp(-kappa[1]*1/12),exp(-kappa[2]*1/12))
  phi <- diag(diag_phi)
  
  # Kalman
  # The Kalman filter loop
  
  loglike <- matrix(c(rep(1,202)),202,1)
  set.seed(123)
  for(i in 1:202){
    # Q matrix -> Put in Thesis
    diag_Q_1 <- alpha^2 / kappa[1]*(exp(-kappa[1]*1/12))-exp(-2*kappa[1]*(1/12))*F_0[1]
    diag_Q_2 <- alpha^2 / kappa[2]*(exp(-kappa[2]*1/12))-exp(-2*kappa[2]*(1/12))*F_0[2]
    diag_Q <- c(diag_Q_1,diag_Q_2)
    Q <- diag(diag_Q)
    
    # Kalman Filter Formulas
    
    F_t_cond <- phi%*%F_0
    P_t_cond <- phi%*%P_0 %*% t(phi)+Q
    u <- t(t(y[i,]))-A-B%*%F_t_cond
    v <- B%*%P_t_cond%*%t(B) + H_m
    loglike[i] <- -0.5*(log(det(v))+t(u)%*%inv(v)%*%(u))
    
    if(i==202)
    {break}
    
    K <- P_t_cond%*%t(B)%*%inv(v)
    L <- 1-K %*%B
    F_0 <- F_t_cond+K_%*%u
    P_0 <- L%*%P_t_cond
  }
  ans <<- sum(loglike)
  return(sum(loglike))
}


## Optimizing vaskal function using the same parameters
####################################################################################################

kappa <- c(0.7,0.004)
alpha <- -5 
H <- c(rep(-2,12))
psi <- 0.02

y <- read.csv("C:/Users/Administrator/Downloads/Maik-Upwork/RFR.csv",1)
print(y)
#y <- swapdata

VasKal_optim <- function(param){
  # Initialize the variables 
  # Take some random valuse 
  
  F_0 <- c(0.05,0.04)
  P_0 <- matrix(c(0.2,0.2,0.2,0.2),2,2)
  tau <- c(1:19)
  tau_n <- length(tau)
  
  # Make Parameters non-negative
  H_vec <- exp(param[4:22])
  alpha_new <- exp(param[3])
  
  kappa <- c(param[1],param[2])
  psi <- param[23]
  
  #Thetas, result of theta* = 0, Vasicek model beta = 0  
  
  theta <- psi*alpha_new/kappa
  theta_1 <- psi*alpha_new/kappa[1]
  theta_2 <- psi*alpha_new/kappa[2]
  
  # H_m stands for H Matrix
  H_m <- diag(H_vec,tau_n,tau_n)
  
  # Formulas for A and B 
  A_factor_1 <- c(theta[1]-0.5*(alpha_new^2/kappa[1]^2))
  A_factor_2 <- c(theta[2]-0.5*(alpha_new^2/kappa[2]^2))
  
  B_tau_1 <- (1-exp(-c(kappa[1]*tau)))/kappa[1]
  B_tau_2 <- (1-exp(-c(kappa[2]*tau)))/kappa[2]
  A_tau_1 <- c(A_factor_1*c(B_tau_1-tau)-c((alpha_new^2*B_tau_1^2)/(4*kappa[1])))
  A_tau_2 <- c(A_factor_2*c(B_tau_2-tau)-c((alpha_new^2*B_tau_2^2)/(4*kappa[2])))
  A_tau <- c(A_tau_1+A_tau_2)
  
  # Creating Observation Equation parameters
  A <- (as.matrix((-A_tau)/tau))
  B_tau <- c(c(B_tau_1/tau),c(B_tau_2/tau))
  
  #tau * factor
  B <- matrix(B_tau,tau_n,2)
  diag_phi <- c(exp(-kappa[1]*1/12),exp(-kappa[2]*1/12))
  phi <- diag(diag_phi)
  

  # The Kalman filter loop
  
  loglike <- matrix(c(rep(1,202)),202,1)
  set.seed(123)
  for(i in 1:202){
    diag_Q_1 <- alpha^2 / kappa[1]*(exp(-kappa[1]*1/12))-exp(-2*kappa[1]*(1/12))*F_0[1]
    diag_Q_2 <- alpha^2 / kappa[2]*(exp(-kappa[2]*1/12))-exp(-2*kappa[2]*(1/12))*F_0[2]
    diag_Q <- c(diag_Q_1,diag_Q_2)
    Q <- diag(diag_Q)
    
    # Kalman Filter Formulas
    
    F_t_cond <- phi%*%F_0
    P_t_cond <- phi%*%P_0 %*% t(phi)+Q
    u <- t(t(y[i,]))-A-B%*%F_t_cond
    v <- B%*%P_t_cond%*%t(B) + H_m
    loglike[i] <- -0.5*(log(det(v))+t(u)%*%inv(v)%*%(u))
    
    if(i==202)
    {break}
    
    K <- P_t_cond%*%t(B)%*%inv(v)
    L <- 1-K %*%B
    F_0 <- F_t_cond+K_%*%u
    P_0 <- L%*%P_t_cond
  }
  ans <<- sum(loglike)
  return(sum(loglike))
}

ParameterVector <- c(kappa,alpha,H,psi)
optim_result  <- optim(ParameterVector,VasKal_optim)
