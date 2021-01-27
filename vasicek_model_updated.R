VasicekHelper <- function(r, kappa, theta, sigma, dt = 1/252) {
  term1 <- exp(-1 * kappa * dt)
  term2 <- (sigma^2) * (1 - term1^2) / (2*kappa)
  result <- r*term1 + theta*(1-term1) + sqrt(term2)*rnorm(n=1)
  return(result)
}

VasicekSimulation <- function(N, r0, kappa, theta, sigma, dt = 1/252) {
  short.rates <- rep(0, N)
  short.rates[1] <- r0
  for (i in 2:N) {
    short.rates[i] <- VasicekHelper(short.rates[i - 1], kappa, theta, sigma, dt)
  }
  return(short.rates)
}

VasicekSimulations <- function(M, N, r0, kappa, theta, sigma, dt = 1/252) {
  sim.mat <- matrix(nrow = N, ncol = M)
  for (i in 1:M) {
    sim.mat[, i] <- VasicekSimulation(N, r0, kappa, theta, sigma, dt)
  }
  return(sim.mat)
}

VasicekZeroCouponBondPrice <- function(r0, kappa, theta, sigma, years) {
  b.vas <- (1-exp(-years*kappa)) / kappa
  a.vas <- (theta-sigma^2/(2*kappa^2))*(years-b.vas)+(sigma^2)/(4*kappa)*b.vas^2
  return(exp(-a.vas-b.vas*r0))
}

VasicekYieldCurve <- function(r0, kappa, theta, sigma, max.maturity=10) {
  yields <- rep(0, max.maturity)
  for (y in 1:max.maturity) {
    yields[y] <- -log(VasicekZeroCouponBondPrice(r0, kappa, theta, sigma, y))/y
  }
  return(yields)
}


VasicekCalibration <- function(fred.ticker = 'RFR', dt = 1/252) {
  
  require(quantmod)
  
  data <- getSymbols(fred.ticker, src = 'RFR', auto.assign = FALSE)
  data <- na.omit(data)/100 
  n <- length(data)
  
  Sx <- sum(data[1:(length(data) - 1)])
  Sy <- sum(data[2:length(data)])
  Sxx <- as.numeric(crossprod(data[1:(length(data) - 1)], data[1:(length(data) - 1)]))
  Sxy <- as.numeric(crossprod(data[1:(length(data) - 1)], data[2:length(data)]))
  Syy <- as.numeric(crossprod(data[2:length(data)], data[2:length(data)]))
  
  theta  <- (Sy * Sxx - Sx * Sxy) / (n* (Sxx - Sxy) - (Sx^2 - Sx*Sy) )
  kappa <- -log((Sxy - theta * Sx - theta * Sy + n * theta^2) /   (Sxx - 2 * theta * Sx + n * theta^2)) / dt
  a <- exp(-kappa*dt)
  sigmah2 <- (Syy - 2 * a * Sxy + a^2 * Sxx - 2 * theta * (1-a) * (Sy - a * Sx) + n * theta^2 * (1 - a)^2)/n
  sigma <- sqrt(sigmah2 * 2 * kappa / (1 - a^2))
  
  r0 <- data[length(data)]
  
  return(c(kappa, theta, sigma, r0))
}


years <- 4
N <- years * 252 
t <- (1:N)/252 

calibration <- VasicekCalibration('RFR')
kappa <- calibration[1]
theta <- calibration[2]
sigma <- calibration[3]
r0 <- calibration[4]

set.seed(666)

test <- VasicekSimulation(N, r0, kappa, theta, sigma)
plot(t, test, type = 'l')

M <- 20
test.mat <- VasicekSimulations(M, N, r0, kappa, theta, sigma)

plot(t, test.mat[, 1], type = 'l', main = 'Short Rates', ylab = 'rt', 
     ylim = c(0, max(test.mat)), col = 1)
for (count in 2:ncol(test.mat)) {
  lines(t, test.mat[, count], col = count)
}

expected <- theta + (r0 - theta)*exp(-kappa*t)
stdev <- sqrt( sigma^2 / (2*kappa)*(1 - exp(-2*kappa*t)))
lines(t, expected, lty=2) 
lines(t, expected + 2*stdev, lty=2) 
lines(t, expected - 2*stdev, lty=2)

VasicekZeroCouponBondPrice(r0, kappa, theta, sigma, years) 

yields <- VasicekYieldCurve(r0, kappa, theta, sigma, 10)
plot(1:10, yields, xlab = 'Maturity', type = 'l', ylab = 'Yield', main = 'Yield Curve')

