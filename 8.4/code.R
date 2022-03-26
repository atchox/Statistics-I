v_square <- 0.04
s_square <- 0.054
n <- 30
chi_square <- (n - 1) * s_square / v_square

# test statistic
print(chi_square)

# quantile of 0.95
print(qchisq(0.95, n - 1))