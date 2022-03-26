x <- c(380, 645, 360, 900, 540, 670, 820, 1050, 760, 800)
y <- c(550, 780, 530, 1200, 620, 800, 910, 1400, 830, 905)

# Scatterplot
png("Scatterplot.png")
plot(y ~ x, main = "Scatterplot of Given Data")
dev.off()

#Linear Model
lm_sales <- lm(y ~ x)
print(lm_sales)

png("ScatterplotWithRegression.png")
plot(y ~ x)
abline(lm_sales$coefficients)
dev.off()

lm_summary <- summary(lm_sales)
# s*s
print(lm_summary$sigma^2)
# SSE
print(lm_summary$sigma^2 * lm_summary$df[2])

# Correlation Coefficient
print(cor(x, y))

# Coefficent of Determination
print(lm_summary$r.squared)

# 95% confidence interval
print(confint(lm_sales, "x"))