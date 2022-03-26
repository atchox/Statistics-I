airline.import <- read.table(
  "https://dasl.datadescription.com/download/data/3048",
  header = T,
  sep = "\t",
  fill = TRUE
)
airline.filter <- airline.import[, c("Rate.16", "Rate.17")]
png(filename = "boxplot.png", width = 600, height = 600)
par(mar = c(3, 3, 4, 3), bg = "#fcfcfc", cex.axis = 1.25, cex.main = 1.5)
boxplot(
  airline.filter,
  main = "Airline Rates",
  names = c("2016", "2017"),
  boxcol = "#204c71", boxfill = "#e7f7ff",
  medcol = "#204c71",
  whiskcol = "#08263f",
  staplecol = "#08263f",
  outcol = "#204c71"
)
dev.off()
rm(list = c("airline.import", "airline.filter"))