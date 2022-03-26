m_0 <- 75
data <- c(74.5, 75.0, 72.3, 76.0, 75.2, 75.1, 75.3, 74.9, 74.8)
n <- length(data)
m <- mean(data)
s <- sd(data) / sqrt(n)
z <- (m - m_0) / s

# test statistic
print(z)

# p-value
print(pt(z, n - 1))