hyper.data <- dhyper(0:5, 100, 200, 5)
binom.data <- dbinom(0:5, 5, 1 / 3)

png("Barplot.png")
barplot(rbind(hyper.data, binom.data), beside = T, names.arg = 0:5)
dev.off()