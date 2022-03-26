sample_count <- 100
lower <- c(0)
upper <- c(0)
count <- 0

for (i in 1:sample_count) {
  x <- rnorm(10)
  y <- rnorm(20, 2)
  diff <- mean(x) - mean(y)
  s <- sqrt((9 * sd(x) + 19 * sd(y)) / 28 * (1 / 10 + 1 / 20))
  t <- qt(1 - 0.05 / 2, 28)
  lower[i] <- diff - t * s
  upper[i] <- diff + t * s
  if (lower[i] < -2 && upper[i] > -2){
    count <- count + 1
  }
}

# Number of confidence intervals which include the true parameter
print(count)

png("Plot.png")
plot(lower, ylim = c(-4, 0), type = "n", xlab = "Samples", ylab = "Intervals")
for (i in 1:sample_count) {
  lines(c(i, i), c(lower[i], upper[i]))
}
abline(-2, 0)
dev.off()