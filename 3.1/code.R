library(HSAUR)
library(qcc)
data("Forbes2000", package = "HSAUR")
countries <- levels(Forbes2000$country)
counts <- rep(0, length(countries))
for (i in seq_len(nrow(Forbes2000))) {
  index <- match(Forbes2000$country[i], countries)
  counts[index] <- counts[index] + 1
}
names(counts) <- countries
png(filename = "countryPareto.png", width = 700, height = 800)
par(mar = c(8, 5, 3, 3))
pareto.chart(counts,
cumperc = seq(0, 100, by = 20),
col = colorRampPalette(c("#204c71", "#ffffff"))(length(c)),
main = "Forbes2000 Countries",
xlab = "Countries",
ylab = "",
cex.lab = 1.5,
mgp = c(11, 1, 0)
)
title(ylab = "Number of Companies", cex.lab = 1.25, mgp = c(3.5, 1, 0))
dev.off()
rm(list = c("countries", "counts", "i", "Forbes2000", "index"))