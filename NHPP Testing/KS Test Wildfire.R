# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)

# Load and clean data
eq.data <- read.csv("/Users/ryanrodrigue/Downloads/wildfireDATA.csv", header=TRUE)
eq.data$DISCOVERY_DATE <- as.Date(eq.data$DISCOVERY_DATE, format="%m/%d/%Y")

# Filter for 2018 and California
eq.data <- subset(eq.data, FIRE_YEAR == 2018 & STATE == "CA")

# Drop rows with missing or malformed time
eq.data <- eq.data[grepl("^\\d{4}$", eq.data$DISCOVERY_TIME), ]

# Parse time
eq.data$DISCOVERY_HOUR <- as.numeric(substr(eq.data$DISCOVERY_TIME, 1, 2))
eq.data$DISCOVERY_MIN <- as.numeric(substr(eq.data$DISCOVERY_TIME, 3, 4))
eq.data$datetime <- as.POSIXct(
  paste(eq.data$DISCOVERY_DATE, sprintf("%02d:%02d:00", eq.data$DISCOVERY_HOUR, eq.data$DISCOVERY_MIN)),
  format="%Y-%m-%d %H:%M:%S"
)

# Sort by datetime
eq.data <- eq.data[order(eq.data$datetime), ]

# Computing lag
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))

# Remove first row
eq.data <- eq.data[-1, ]

# Interarrival times (hours)
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600

### MODEL lambda(t) ###

# Create year-week variable using ISOweek
eq.data$year.week <- paste0(year(eq.data$datetime), "-W", sprintf("%02d", isoweek(eq.data$datetime)))

# Remove events from the last (outlier) week
outlier_week <- tail(sort(unique(eq.data$year.week)), 1)
eq.data <- eq.data[eq.data$year.week != outlier_week, ]

# Number of earthquakes per week
freq.week <- data.frame(table(eq.data$year.week))

# Remove last week (outlier) from summary table
freq.week <- freq.week[-nrow(freq.week), ]

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

# Normal distribution fit
normal_model <- nls(lambda ~ a * dnorm((median.time - mu)/sigma),
                    start = list(a = max(lambda), mu = median(median.time), sigma = sd(median.time)))
normal_coefs <- coef(normal_model)

# Define estimated lambda(t) function
lambda.fn <- function(t) {
  normal_coefs["a"] * dnorm((t - normal_coefs["mu"]) / normal_coefs["sigma"])
}

### Time-rescaling ###

# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)
eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))

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
plot(median.time, lambda,
     xlab = "Days since January 1, 2018", 
     ylab = "Wildfires per Day", 
     pch = 1)

lines(median.time, lambda.fn(median.time), lwd=2, col="blue")

legend("topright", 
       legend = c("Observed Daily Intensity (per week)", "Estimated λ(t)"), 
       col = c("black", "blue"), 
       pch = c(1, NA), 
       lty = c(NA, 1), 
       lwd = c(1, 2),
       pt.cex = 1.2,
       bty = "o")

