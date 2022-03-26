brake1 <- c(43, 44, 57, 41, 47, 54)
brake2 <- c(51, 65, 62, 67, 58, 57)
brake3 <- c(34, 45, 39, 37, 48, 38)

brakes.life <- c(brake1, brake2, brake3)
brakes.num <- rep(1:3, each = 6)
brakes <- data.frame(
  num = as.factor(brakes.num),
  life = brakes.life
)

# ANOVA table
print(summary(aov(life ~ num, data = brakes)))