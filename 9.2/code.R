s_1 <- 1.87
n_1 <- 10
s_2 <- 1.25
n_2 <- 10

F <- s_1^2/s_2^2

# Test statistic
print(F)

# 5% and 95% quantiles
print(qf(0.975, n_1 - 1, n_2 - 1))
print(qf(0.025, n_1 - 1, n_2 - 1))