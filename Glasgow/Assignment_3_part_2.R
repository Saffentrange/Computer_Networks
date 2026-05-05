### 3.2.1

# Alcohol as dependent behaviour variable
alcoholbeh <- as_dependent_rsiena(
  as.matrix(alcohol_data[, 1:3]), 
  type = "behavior"
)

# Rebuild data object with alcohol as dependent variable
mydata_coev <- make_data_rsiena(friendship, alcoholbeh, gender, logdistance)

# Check preconditions
write_report(mydata_coev, outputName = "glasgowCoevReport")

# Moran's I for autocorrelation in alcohol behaviour
moran1 <- nacf(f1, alcohol_data[, 1], lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")
moran2 <- nacf(f2, alcohol_data[, 2], lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")
moran3 <- nacf(f3, alcohol_data[, 3], lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")
autocorr <- rbind(moran1, moran2, moran3)
autocorr[, 2]

## Model specification
myeff_coev <- make_specification(mydata_coev)

# --- Selection part: identical to 3.1.2 ---
myeff_coev <- set_effect(myeff_coev, list(recip),       depvar = "friendship")
myeff_coev <- set_effect(myeff_coev, list(transTrip),   depvar = "friendship")
myeff_coev <- set_effect(myeff_coev, list(cycle3),      depvar = "friendship")
myeff_coev <- set_effect(myeff_coev, list(egoX, altX, sameX), depvar = "friendship",
                         covar1 = "gender")

# Hypothesis i: popularity
myeff_coev <- set_effect(myeff_coev, list(inPopSqrt),   depvar = "friendship")

# Hypothesis ii: alcohol homophily (now using the behaviour dependent variable)
myeff_coev <- set_effect(myeff_coev, list(egoX, altX, simX), depvar = "friendship",
                         covar1 = "alcoholbeh")

# Hypothesis iii: geographic proximity
myeff_coev <- set_effect(myeff_coev, list(X), depvar = "friendship",
                         covar1 = "logdistance")

# --- Influence part: effects explaining alcohol behaviour evolution ---

# Basic controls for behaviour (linear and quadratic shape)
# These are included by default in make_specification but listed here for clarity
# They capture the baseline tendency to increase/decrease alcohol consumption

# Hypothesis iv: Popular students tend to increase or maintain alcohol consumption
# indeg: students who receive more friendship nominations (more popular)
# tend to increase their alcohol consumption
myeff_coev <- set_effect(myeff_coev, list(indeg), depvar = "alcoholbeh",
                         covar1 = "friendship")

# Hypothesis v: Students adjust alcohol consumption to that of their friends
# avSim: average similarity effect — ego moves their behaviour toward
# the average alcohol level of their alters (the core social influence effect)
myeff_coev <- set_effect(myeff_coev, list(avSim), depvar = "alcoholbeh",
                         covar1 = "friendship")

myeff_coev

## Model estimation
SAOMAlgorithm_coev <- sienaAlgorithmCreate(seed = 161, n3 = 3000, nsub = 4,
                                           MaxDegree = c(friendship = 6))

model_coev <- siena07(SAOMAlgorithm_coev,
                      data    = mydata_coev,
                      effects = myeff_coev,
                      returnDeps = TRUE,
                      batch = FALSE)

# Convergence check
t_ratios_coev <- apply(model_coev$sf, 2, mean) / apply(model_coev$sf, 2, sd)
t_ratios_coev
overall_coev <- sqrt(t(apply(model_coev$sf, 2, mean)) %*% 
                       solve(cov(model_coev$sf)) %*% 
                       apply(model_coev$sf, 2, mean))
overall_coev

## Goodness of fit
cl <- makeCluster(4)

gof_coev_id <- test_gof(model_coev, verbose = FALSE,
                        varName = "friendship", IndegreeDistribution, cluster = cl)
gof_coev_od <- test_gof(model_coev, verbose = FALSE,
                        varName = "friendship", OutdegreeDistribution, cluster = cl)
gof_coev_tc <- test_gof(model_coev, verbose = FALSE,
                        varName = "friendship", TriadCensus, cluster = cl)
gof_coev_gd <- test_gof(model_coev, verbose = FALSE,
                        varName = "friendship", GeodesicDistribution)
gof_coev_beh <- test_gof(model_coev, verbose = FALSE,
                         varName = "alcoholbeh", BehaviorDistribution, cluster = cl)

stopCluster(cl)

plot(gof_coev_id)
plot(gof_coev_od)
plot(gof_coev_gd)
plot(gof_coev_tc, center = TRUE, scale = TRUE)
plot(gof_coev_beh)

model_coev