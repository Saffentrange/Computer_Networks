library(manynet)
library(network)
library(ergm)
library(ggplot2)
library(tidyr)

# i
f1_matrix <- as.matrix(read.csv("f1.csv", header = FALSE))
glasgow_net <- as.network(f1_matrix, directed = TRUE)

max_k <- max(rowSums(f1_matrix))

m1 <- ergm(glasgow_net ~ edges)

# ii
m2_neg <- ergm(glasgow_net ~ edges + gwodegree(1.0, fixed = TRUE))
m2_pos <- ergm(glasgow_net ~ edges + gwodegree(1.0, fixed = TRUE))

# iii
sim_out_dist <- function(formula, coefficients, n_sim = 100, max_val) {
  sims <- simulate(formula, coef = coefficients, nsim = n_sim)
  
  dist_list <- lapply(sims, function(x) {
    deg_counts <- table(factor(rowSums(as.matrix(x)), levels = 0:max_val))
    return(as.data.frame(t(as.matrix(deg_counts))))
  })
  
  return(do.call(rbind, dist_list))
}

set.seed(123)

sim_m1_df <- sim_out_dist(glasgow_net ~ edges, 
                          coef(m1), 100, max_k)

sim_m2_neg_df <- sim_out_dist(glasgow_net ~ edges + gwodegree(1.0, fixed = TRUE), 
                              coef(m2_neg), 100, max_k)

sim_m2_pos_df <- sim_out_dist(glasgow_net ~ edges + gwodegree(1.0, fixed = TRUE), 
                              coef(m2_pos), 100, max_k)

process_df <- function(df, label) {
  df %>%
    pivot_longer(cols = everything(), names_to = "degree", values_to = "node_count") %>%
    mutate(model = label, degree = as.numeric(degree))
}

plot_data <- rbind(
  process_df(sim_m1_df, "M1 (Edges Only)"),
  process_df(sim_m2_neg_df, "M2 (GWOD Negative)"),
  process_df(sim_m2_pos_df, "M2 (GWOD Positive)")
)

ggplot(plot_data, aes(x = factor(degree), y = node_count, fill = model)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~model) +
  labs(title = "Out-Degree Distribution Comparison",
       subtitle = "Alpha = 1.0, 100 Simulations (No Edge Constraints)",
       x = "Out-degree (k)",
       y = "Number of Nodes") +
  theme_minimal() +
  theme(legend.position = "none")
