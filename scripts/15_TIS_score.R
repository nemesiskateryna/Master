# 1. Extract TIS score

tis_values <- immune_scores$TIS_mean_zscore

# 2. Summary statistics

cat("\n--- TIS DISTRIBUTION ---\n")

print(
  summary(tis_values)
)

cat(
  "Standard deviation:",
  sd(tis_values),
  "\n"
)

cat(
  "Skewness:",
  mean(
    (tis_values - mean(tis_values))^3
  ) / sd(tis_values)^3,
  "\n"
)

# 3. Histogram

hist(
  tis_values,
  breaks = 30,
  main = "Distribution of TIS scores in TCGA-SKCM",
  xlab = "TIS mean z-score",
  ylab = "Number of samples"
)

# 4. Density plot

plot(
  density(tis_values),
  main = "Density distribution of TIS scores",
  xlab = "TIS mean z-score"
)

# 5. Median

cat(
  "TIS median:",
  median(tis_values),
  "\n"
)
