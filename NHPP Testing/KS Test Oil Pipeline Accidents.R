# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)

# Load data
eq.data <- read.csv("oilpipelineaccidents.csv", header=TRUE, sep=",")
eq.data$datetime <- as.POSIXct(eq.data$Accident.Date.Time, format="%m/%d/%Y %I:%M %p")

# Filter for dates from Jan to June 2016
eq.data <- eq.data %>%
  filter(datetime >= as.POSIXct("2016-01-01", tz = "UTC") & 
         datetime < as.POSIXct("2016-07-01", tz = "UTC"))


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

# Use constant lambda: average intensity
lambda.const <- mean(lambda, na.rm = TRUE)

# Define constant lambda(t)
lambda.fn <- function(t) {
  return(as.numeric(lambda.const))
}

### Time-rescaling ###

# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)
eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units = "days"))

# Rescale each event time by integrating lambda from 0 to t
rescaled_times <- eq.data$time_in_days * lambda.const

# Compute rescaled interarrival times
rescaled_interarrivals <- diff(c(0, rescaled_times))  # Add 0 to start

### KS Goodness-of-fit test for Exp(1) ###
ks_result <- ks.test(rescaled_interarrivals, "pexp", 1)
print(ks_result)

### Plotting λ(t) ###
plot(median.time, lambda,
     xlab = "Days since January 1, 2016", 
     ylab = "Oil Pipeline Accidents per Day", 
     pch = 1)

lines(median.time, rep(lambda.const, length(median.time)), lwd = 2, col = "blue")

legend("topright", 
       legend = c("Observed Daily Intensity (per week)", "Estimated λ(t)"),
       col = c("black", "blue"),
       pch = c(1, NA),
       lty = c(NA, 1),
       lwd = c(1, 2),
       pt.cex = 1.2,
       bty = "o")
