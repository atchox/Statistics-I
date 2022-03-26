data.import <- read.csv(
  "https://query1.finance.yahoo.com/v7/finance/download/%5ENSEI?period1=1604725694&period2=1636261694&interval=1d&events=history&includeAdjustedClose=true"
)
data.filter <- data.import[data.import$Close != "null", c("Date", "Close")]
data.filter[, "Close"] <- sapply(data.filter[, "Close"], as.double)
n <- nrow(data.filter)
r <- log(data.filter$Close[2:n]) - log(data.filter$Close[2:n - 1])
png("Histogram.png")
hist(r, breaks = 25, freq = F)
mean <- mean(r)
sd <- sd(r)
dev.off()
png("NormalPlot.png")
qqnorm(r)
dev.off()
library(timeDate)
observed <- c(
  skewness(r),
  kurtosis(r),
  IQR(r) / sd,
  length(which(r < mean + sd & r > mean - sd)) * 100 / length(r),
  length(which(r < mean + 2 * sd & r > mean - 2 * sd)) * 100 / length(r),
  length(which(r < mean + 3 * sd & r > mean - 3 * sd)) * 100 / length(r)
)
expected <- c(0, 3, 1.3, (pnorm(1:3) - pnorm(-c(1:3)))*100)
summary <- data.frame(observed = observed, expected = expected)
row.names(summary) <- c(
  "Skewness",
  "Kurtosis",
  "IQR:SD",
  "Within 1 SD",
  "Within 2 SD",
  "Within 3 SD"
)
print(summary)