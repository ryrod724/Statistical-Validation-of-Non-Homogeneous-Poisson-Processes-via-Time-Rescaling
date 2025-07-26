# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)
library(readr)

# Load data
eq.data <- read_csv("/Users/ryanrodrigue/Downloads/tnUS.csv") %>%
  filter(st %in% c("OK", "TX", "KS", "NE"), yr %in% c(2020))

eq.data$datetime <- as.POSIXct(eq.data$datetime_utc, format = "%Y-%m-%d %H:%M:%S")

# Computing lag
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))

# Remove first row
eq.data <- eq.data[-1, ]

# Interarrival times (hours)
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600

### MODEL lambda(t) ###

# Create year-week variable using ISOweek
eq.data$year.week <- paste0(year(eq.data$datetime), "-W", sprintf("%02d", isoweek(eq.data$datetime)))

# Number of earthquakes per week
freq.week <- data.frame(table(eq.data$year.week))

# Weekly variables
year.week.unique <- freq.week[,1]
neq.week <- freq.week[,2]
ndays.week <- rep(7, length(neq.week))  # Assume each week is 7 days

# Estimate intensity per day
lambda <- neq.week / ndays.week

# Cumulative number of days to each week (for plotting x-axis)
median.time <- c()
ndays.total <- c()
median.time[1] <- ndays.week[1] / 2
ndays.total[1] <- ndays.week[1]
for (i in 2:length(ndays.week)) {
  median.time[i] <- ndays.total[i-1] + ndays.week[i]/2
  ndays.total[i] <- ndays.total[i-1] + ndays.week[i]
}
median.time <- as.numeric(median.time)

# Estimate lambda(t) with Weibull-like function via nonlinear least squares
weibull.fn <- function(t, a, b, c) {
  a * (t / b)^(c - 1) * exp(-(t / b)^c)
}
weibull.model <- nls(
  lambda ~ a * (median.time / b)^(c - 1) * exp(-(median.time / b)^c),
  start = list(a = 1, b = 100, c = 2),
  control = nls.control(maxiter = 500)
)
weibull.coefs <- coef(weibull.model)

# Define new λ(t)
lambda.fn <- function(t) {
  weibull.fn(t, weibull.coefs["a"], weibull.coefs["b"], weibull.coefs["c"])
}

### Time-rescaling ###
# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)
eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units = "days"))

# Rescale each event time by integrating lambda from 0 to t
rescaled_times <- sapply(eq.data$time_in_days, function(t) {
  if (is.finite(t) && !is.na(t)) {
    return(integrate(lambda.fn, lower = 0, upper = t)$value)
  } else {
    return(NA)
  }
})

# Compute rescaled interarrival times
rescaled_interarrivals <- diff(c(0, rescaled_times))  # Add 0 to start

### KS Goodness-of-fit test for Exp(1) ###
ks_result <- ks.test(rescaled_interarrivals, "pexp", 1)
print(ks_result)

### Plotting λ(t) ###
plot(median.time, lambda, xlab = "Days since January 1, 2020", ylab = "Tornadoes per Day", pch = 1)
lines(median.time, lambda.fn(median.time), lwd = 2, col = "blue")
legend("topright",
       legend = c("Observed Daily Intensity (per week)", "Estimated λ(t)"),
       col = c("black", "blue"),
       pch = c(1, NA),
       lty = c(NA, 1),
       lwd = c(1, 2),
       pt.cex = 1.2,
       bty = "o")
