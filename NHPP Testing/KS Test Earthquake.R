# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)

# Load data
eq.data <- read.csv("/Users/ryanrodrigue/Downloads/earthquakes.csv", header=TRUE, sep=",")
eq.data$datetime <- as.POSIXct(eq.data$datetime, format="%Y-%m-%d %H:%M:%S")

# Filter to October-December 2024 only (most recent 3 month period)
eq.data <- eq.data[290257:294855, ]

# Computing lag
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))

# Remove first row
eq.data <- eq.data[-1, ]

# Interarrival times (hours)
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600

# Remove immediate aftershocks (within 2 hours)
eq.data <- eq.data[eq.data$elapsed.time > 2, ]

### MODEL lambda(t) ###

### MODEL lambda(t) ###

# Create year-week variable using ISOweek (from lubridate::isoweek)
eq.data$year.week <- paste0(year(eq.data$datetime), "-W", sprintf("%02d", isoweek(eq.data$datetime)))

# Number of earthquakes per week
freq.week <- data.frame(table(eq.data$year.week))
year.week.unique <- freq.week[,1]
neq.week <- freq.week[,2]

# Number of days per week (all are 7, but let's be robust)
ndays.week <- rep(7, length(neq.week))

# Estimate intensity per day
lambda <- neq.week / ndays.week

# Cumulative number of days to each week
median.time <- c()
ndays.total <- c()
median.time[1] <- ndays.week[1] / 2
ndays.total[1] <- ndays.week[1]
for (i in 2:length(ndays.week)) {
  median.time[i] <- ndays.total[i-1] + ndays.week[i]/2
  ndays.total[i] <- ndays.total[i-1] + ndays.week[i]
}
median.time <- as.numeric(median.time)

# Polynomial regression (same)
median.time.re <- median.time / 1000
median.time.sq <- median.time.re^2
median.time.cu <- median.time.re^3

model <- lm(lambda ~ median.time.re + median.time.sq + median.time.cu)
coefs <- coef(model)

# Define estimated lambda function
lambda.fn <- function(t) {
  coefs[1] + coefs[2]*(t/1000) + coefs[3]*(t/1000)^2 + 
    coefs[4]*(t/1000)^3
}

### Time-rescaling ###

# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)

eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))

# Now rescale each event time by integrating lambda from 0 to t
# Ensure t is valid before integrating
rescaled_times <- sapply(eq.data$time_in_days, function(t) {
  if (is.finite(t) && !is.na(t)) {
    return(integrate(lambda.fn, lower = 0, upper = t)$value)
  } else {
    return(NA)  # Return NA if t is not valid
  }
})

# Compute rescaled interarrival times
rescaled_interarrivals <- diff(c(0, rescaled_times))  # Add 0 to start

### KS Goodness of fit via test for Exp(1) ###

ks_result <- ks.test(rescaled_interarrivals, "pexp", 1)

# Show results
print(ks_result)
      
#plotting lambda approximation
plot(median.time, lambda, xlab = "Days since October 1, 2024", ylab = "Earthquakes per Week")

lines(median.time, lambda.fn(median.time), lwd=2, col="blue")

legend("topright", 
       legend = c("Observed Weekly Intensity", "Estimated λ(t)"), 
       col = c("black", "blue"), 
       pch = c(1, NA), 
       lty = c(NA, 1), 
       lwd = c(1, 2),
       pt.cex = 1.2,
       bty = "o")
