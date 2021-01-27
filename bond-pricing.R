# To evaluate bond values using R
# Create present value or use some data here 
pv <- 100

# Create r, which is interest rate
r <- 0.1

# Calculate future value after 1st yr
fv1 <- pv * (1 + r)

# Calculate future value after 2nd yr
fv2 <- pv*(1+r)*(1+r)

# Create vector of cash flows
cf <- c(5,5,5,5,105) 

# Convert to data frame
cf <- data.frame(cf)

# Add column t to cf; t indicates the year in which cash flow is received
cf$t <- as.numeric(rownames(cf))

# Set the yield rate to be 0.06
yr <- 0.06

# Calculate pv_factor;this create a column of discount factors
cf$pv_factor <- 1 / (1 + 0.06)^cf$t

# Calculate pv
cf$pv <- cf$cf * cf$pv_factor

# Calculate the bond price, which is the summation of the present values of future cash flows
sum(cf$pv)


# Here to evalue the bond value through a function
# Create function
# p: par value; r: coupon rate; ttm: time to maturity; y: yield
bondprc <- function(p, r, ttm, y) {
  # rep returns a vector with value = p * r and times = ttm -1 
  cf <- c(rep(p * r, ttm - 1), p * (1 + r))
  cf <- data.frame(cf)
  cf$t <- as.numeric(rownames(cf))
  cf$pv_factor <- 1 / (1 + y)^cf$t
  cf$pv <- cf$cf * cf$pv_factor
  sum(cf$pv)
}

# Verify prior result
bondprc(100,0.05,5,0.06)

# Install and Load Quandl package
install.packages("Quandl")
library(Quandl)

# Obtain Moody's Baa index data
baa <- Quandl("MOODY/DBAAYLD")

# Identify 7/24/17 yield
baa_yield <- subset(baa, baa$DATE == "2017-07-24")

# Convert yield to decimals and view, since the data is in percentage
baa_yield <-baa_yield$DBAA *0.01
baa_yield

# Value bond with the following parameters
bondprc(p = 100, r =0.05, ttm =5, y = 0.0429)

# Plotting the Price/Yield Relationship
# Generate prc_yld
prc_yld <- seq(0.02, 0.4, 0.01)

# Convert prc_yld to data frame
prc_yld <- data.frame(prc_yld)

# Calculate bond price given different yields
for (i in 1:nrow(prc_yld)) {
  prc_yld$price[i] <- bondprc(100, 0.10, 20, prc_yld$prc_yld[i])  
}

# Plot P/YTM relationship. It can observe from the graph that bond price and yield has an inverse relationship
plot(prc_yld,prc_yld$price,
     type = "l",
     col = "blue",
     main = "Price/YTM Relationship")

#=========================================================================================
# Plot 10-Year US Treasury Yields
# Load quantmod package
library(quantmod)
library(xts)
library(zoo)
# Obtain Treasury yield data
t10yr <- getSymbols(Symbols = "DGS10", src = "FRED", auto.assign = FALSE)
#t10yr
# Subset data
head(t10yr)
t10yr <- subset(t10yr["2006-01-01/2016-09-30"] )

# Plot yields
plot(x = index(t10yr),
     y = t10yr$DGS10,
     xlab = "Date",
     ylab = "Yield (%)",
     type = "l",
     col = "red",
     main = "10-Year US Treasury Yields")


# =======================================================================================
# Ploting investing grade yield: Aaa - Baa spread. Data source (Moody Yield)
# Examine first and last six elements in spread
spread_aaa <- Quandl("MOODY/DAAAYLD")
spread_aaa <- subset(spread_aaa, spread_aaa$DATE>"1986-01-03" & spread_aaa$DATE <"2017-08-03")
spread_baa <- Quandl("MOODY/DBAAYLD")
spread_baa <- subset(spread_baa, spread_baa$DATE>"1986-01-03" & spread_baa$DATE <"2017-08-03")
spread <- merge(spread_aaa, spread_baa, by="DATE")
head(spread)
tail(spread)

# Calculate spread$diff
spread$diff <- (spread$DBAA - spread$DAAA)*100

# Plot spread
plot(x = spread$DATE,
     y = spread$diff,
     type = "l",
     xlab = "Date",
     ylab = "Spread (bps)",
     col = "red",
     main = "Aaa - Baa Spread")

#========================================================================================
# Find a bond's yield
# Value bond using 5% yield
bondprc(p = 100, r = 0.05, ttm = 5, y = 0.05)

# Value bond using 7% yield
bondprc(p = 100, r = 0.05, ttm = 5, y = 0.07)

# Value bond using 6% yield
bondprc(p = 100, r = 0.05, ttm = 5, y = 0.06)

# Using uniroot function to find YTM
# Create cash flow vector
cf <- c(-95.79, 5,5,5,5,105)

# Create bond valuation function
bval <- function(i, cf,
                 t=seq(along = cf))
  sum(cf / (1 + i)^t)

# Create ytm() function using uniroot
ytm <- function(cf) {
  uniroot(bval, c(0, 1), cf = cf)$root
}

# Use ytm() function to find yield
ytm(cf)


# Calculate the PV01
PV01 <- abs(bondprc(100,0.1,20,0.1001) - bondprc(100,0.1,20,0.1000))
PV01

#================================================================================
# Calculate approximate duration for a bond
# Calculate bond price today
px <- bondprc(p = 100, r = 0.1, ttm = 20, y = 0.1)
px

# Calculate bond price if yields increase by 1%
px_up <- bondprc(p = 100, r = 0.1, ttm = 20, y = 0.11)
px_up

# Calculate bond price if yields decrease by 1%
px_down <- bondprc(p = 100, r = 0.1, ttm = 20, y = 0.09)
px_down

# Calculate approximate duration
duration <- (px_down - px_up) / (2* px * 0.01)
duration


# Estimating effect on bond price using duration
# Estimate percentage change
duration_pct_change <- -duration * (-0.01)
duration_pct_change

# Estimate dollar change
duration_dollar_change <- px*duration_pct_change
duration_dollar_change

#==================================================================================
# Calculate approximate convexity for a bond
# Calculate approximate convexity
convexity <- (px_down + px_up - 2 * px) / (px * (0.01)^2)
convexity

# Estimate percentage change
convexity_pct_change <- 0.5 * convexity * (0.01)^2
convexity_pct_change

# Estimate dollar change
convexity_dollar_change <- convexity_pct_change * px 
convexity_dollar_change

# Estimating the bond price using duration and convexity
# Estimate change in price
price_change <- duration_dollar_change + convexity_dollar_change

# Estimate price
price <- duration_dollar_change + convexity_dollar_change + px

#==================================================================================
# In this comprehensive example, you will value a bond with a $100 par value, 
# 3% coupon rate, and 8 years to maturity. This bond was rated Aaa by Moody's 
# and it was issued on September 30, 2016. You have determined that this bond's 
# yield is comparable to the yield of bonds with a Aaa rating.

# Load Quandl package
library(Quandl)

# Obtain Moody's Aaa yield
aaa <- Quandl("MOODY/DAAAYLD")

# identify yield on September 30, 2016
aaa_yield <- subset(aaa, aaa$DATE == "2016-09-30")

# Convert yield into decimals
aaa_yield <- as.numeric(aaa_yield$VALUE)*0.01
aaa_yield

# Layout the bond's cash flows
cf <- c(3,3,3,3,3,3,3,103)

# Convert to data.frame
cf <- data.frame(cf)

# Add time indicator
cf$t <- seq(1, 8, 1)

# Calculate PV factor
cf$pv_factor <- 1 / (1 + aaa_yield)^cf$t

# Calculate PV
cf$pv <- cf$cf * cf$pv_factor

# Price bond
sum(cf$pv)

# Code cash flow function
alt_cf <- function(r, p, ttm) {
  c(rep(p * r, ttm - 1), p * (1 + r))
}

# Generate cf vector
alt_cf(r = 0.03, p = 100, ttm = 8)

# Calculate bond price when yield increases

px_up <- bondprc(p = 100, r =0.03, ttm = 8, y = aaa_yield + 0.01)

# Calculate bond price when yield decreases
px_down <- bondprc(p = 100, r =0.03, ttm = 8, y = aaa_yield - 0.01)

# Calculate duration
duration <- (px_down - px_up) / (2 * px* 0.01)

# Calculate percentage effect of duration on price
duration_pct_change <- -duration * 0.01
duration_pct_change

# Calculate dollar effect of duration on price
duration_dollar_change <- duration_pct_change * px
duration_dollar_change

# Calculate convexity measure
convexity <- (px_up + px_down - 2*px)/(px * (0.01)^2)

# Calculate percentage effect of convexity on price
convexity_pct_change <- convexity * 0.5* (0.01)^2
convexity_pct_change

# Calculate dollar effect of convexity on price
convexity_dollar_change <- convexity_pct_change * px
convexity_dollar_change

# Estimate price_change
price_change <- duration_dollar_change + convexity_dollar_change
price_change

# Estimate new_price
new_price <- price_change + px
new_price
