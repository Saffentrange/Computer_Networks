install.packages(
  c("network", "ergm", "manynet", "patchwork", "latticeExtra", "texreg", "broom")
)

remotes::install_github("stocnet/autograph", ref = "develop")

library(manynet)
library(autograph)
library(network)
library(ergm)
library(ggraph)
library(patchwork)
library(texreg)
library(broom)

f1 <- as.matrix(read.csv("f1.csv", header = FALSE))
f2 <- as.matrix(read.csv("f2.csv", header = FALSE))
f3 <- as.matrix(read.csv("f3.csv", header = FALSE))

demographic <- read.csv("demographic.csv")
logdistance <- as.matrix(read.csv("logdistance.csv", header = FALSE))
alcohol <- read.csv("alcohol.csv")

# out degree range from f1
max_k <- max(rowSums(f1))
k_vals <- 0:max_k

# change statistic function
change_stat <- function(k, a) {
  if (a == 0) {
    return(ifelse(k==0, 1, 0))
  }
  return((1 - exp(-a))^k)
}

# alpha values to test
decay_val_neg <- c(-1, -0.5)
decay_val_pos <- c(0, 0.1, 0.2, 0.3, 0.5, log(2), 1, 2, 3)

par(mfrow = c(1,2))
for (a in decay_val_neg) {
  y_vals <- sapply(k_vals, change_stat, a = a)
  plot(k_vals, y_vals, type = "b", pch = 19, col = "red",
       main = paste("Alpha =", a), xlab = "Out-degree (k)", ylab = "Change Statistic")
  abline(h = 0, lty = 2)
}

par(mfrow = c(1, 1))
colors <- rainbow(length(decay_val_pos))
plot(NULL, xlim = c(0, max_k), ylim = c(0, 1), xlab = "Out-degree (k)",
     ylab = "Change Statistic", main = "GWD Change Statistic (Non-negative Alpha)")

for (i in 1:length(decay_val_pos)) {
  y_vals <- sapply (k_range, delta_gwd, alpha = alphas_pos[i])
  lines(k_range, y_vals, type = "b", pch = 19, col = colors[i])
}

legend("topright", legend = round(decay_val_pos, 2), col = colors, lty = 1, pch = 19, title = "Alpha", cex = 0.8)
