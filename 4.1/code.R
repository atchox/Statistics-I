n <- 10000
lambda <- 2
uniform_data <- runif(n)
exp_data <- -log(uniform_data) / lambda
hist_breaks <- seq(min(exp_data), max(exp_data) + 1, length.out = 25)
curve_breaks <- seq(min(exp_data), max(exp_data) + 1, length.out = 50)
png(filename = "histogram.png", width = 700, height = 700)
hist_data <- hist(exp_data,
  hist_breaks,
  freq = F,
  col = "#ebebeb",
  main = "Histogram of Exponential Data",
  xlab = "Data",
  ylim = c(0, 2))
lines(
  x = curve_breaks,
  y = dexp(curve_breaks, lambda),
  col = "#77120f")
dev.off()
rm(list = c("n",
  "lambda",
  "uniform_data",
  "exp_data",
  "hist_breaks",
  "curve_breaks",
  "hist_data")
)