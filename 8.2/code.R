m_0 <- 11
data <- c(12.5, 13.0, 11.9, 10.2, 13.1, 13.6, 13.8, 14.0)
n <- length(data)
m <- mean(data)
s <- sd(data) / sqrt(n)
z <- (m - m_0) / s

# test statistic
print(z)

# p-level
print(1 - pt(z, n - 1))