---
  title: "Network Modeling: Assignment 1"
output: pdf_document
---
  
  **Group number**: 21\
**Group members**: Eduard Mihai, Hanzhang Wang, Silvan Affentranger, Linus Michel

# Task 1: Tie dependence and CUGs

## 1.1

We compute $p$ by dividing the observed number of ties by the total number of pairs:
  
  $$
  p=\frac{134}{22 \cdot 21}=134/462\approx0.290
$$
  
  ## 1.2
  
  We compute the probabilities by dividing the observed count of each dyad state by the total number of dyads:
  
  $$
  \begin{aligned}
p_M=\frac{M}{M+A+N}=42/231\approx0.182 \\
p_A=\frac{A}{M+A+N}=50/231\approx0.216 \\
p_N=\frac{N}{M+A+N}=139/231\approx0.602
\end{aligned}
$$
  
  ## 1.3
  
  If we assume independence and use the p in (1), we have:
  
  $$
  \begin{aligned}
\text{Both ties exist:}\quad & p_M' = p^2 = (134/462)^2 \approx 0.084 \\
    \text{Exactly one tie exists:}\quad & p_A' = 2 \cdot p \cdot (1-p) = 2 \cdot (134/462) \cdot (1 - 134/462) \approx 0.412 \\
\text{Neither tie exists:}\quad & p_N' = (1-p)^2 = (1 - 134/462)^2 \approx 0.504
  \end{aligned}
$$

## 1.4

It is not reasonable to assume tie independence for this network.

Comparing the observed values computed at (2) with the expected values computed at (3), we observe that the proportion of mutual dyads is much higher than expected under independence, while the proportion of asymmetric dyads is substantially lower. The proportion of null dyads is also higher than predicted.

This indicates that ties are not independent. In particular, there is evidence of positive reciprocity, meaning that the presence of a tie from $i$ to $j$ increases the likelihood of a tie from $j$ to $i$. Therefore, the assumption of tie independence is not appropriate for further analysis of this network.

## 1.5

If every node is forced to have an outdegree of exactly 9, the probability $\bar{p}$ of a tie appearing in an arbitrary pair becomes:

$$
  \bar{p}=\frac{9 \cdot 22/2}{22 \cdot 21/2}=3/7\approx0.429
$$

Assuming tie independence under this framework, the probability of a mutual tie would be:

$$
  \bar{p_M}=\bar{p}^2=(3/7)^2\approx0.184
$$

While this modeled mutual tie probability ($0.184$) is very close to the observed mutual tie probability computed at (2) ($0.182$), this is still not a reasonable model to assume overall. If all 22 nodes have an outdegree of 9, the total number of ties in the network would have to be $22 \cdot 9 = 198$. However, the observed network only has $m = 134$ ties.

## 1.6

```{r message=FALSE}
library(sna)
library(network)
obsMat <- as.matrix(read.csv("./Lintner/10_W1.csv", sep = ";", row.names = 1))

obsMat[is.na(obsMat)] <- 0
stopifnot(sum(is.na(obsMat)) == 0)

set.seed(161) 
cguRec <- cug.test(obsMat, grecip, cmode = "edges", reps = 5000)
cugInd <- cug.test(
  obsMat,
  centralization,
  cmode = "edges",
  FUN.arg = list(FUN = degree, cmode = "indegree"),
  reps = 5000
)
cugTrans <- cug.test(obsMat, gtrans, cmode = "dyad.census", reps = 5000)
```

The null hypothesis is $H_0$: *"The observed reciprocity is consistent with the random graph with the same number of edges"*. The conditional feature is the number of edges. From the below diagram, in the simulations, no random network has equal or more mutual dyads, and the p-value is $0$. Under a significance level of $\alpha=0.05$, we consider this result significant. Namely, reciprocity is much higher than expected by chance, and the null hypothesis is rejected.

```{r}
plot(cguRec)
paste("p-value =", sprintf("%.3f", cguRec$pgteobs))
```

The null hypothesis is $H_0$: *"The observed indegree centralization is consistent with the random graph with the same number of edges"*. The conditional feature is the number of edges. From the below diagram, a large part of random networks has equal or more mutual dyads, and the p-value is $0.720$, above a significance level of $\alpha=0.05$, so we consider this result not significant. We fail to reject the null hypothesis, suggesting that the network does not display a significant level of hierarchical structure.

```{r}
plot(cugInd)
paste("p-value =", sprintf("%.3f", cugInd$pgteobs))
```

The null hypothesis is $H_0$: *"The observed transitivity is consistent with random graphs having the same dyad census"*. The conditional features consist of the number of mutual dyads, asymmetric dyads, and null dyads. From the below diagram, none of the random networks has equal or more transitive triads, namely the p-value is $0$, under a significance level of $\alpha=0.05$. We consider this result significant, and the null hypothesis is rejected.

```{r}
plot(cugTrans)
paste("p-value =", sprintf("%.3f", cugTrans$pgteobs))
```

# Task 2: MR-QAP Regression

## 2.1

We use wave 2 as the dependent variable, and include wave 1 as one explanatory variable to build the QAP regression model.

```{r, warning=FALSE}
library(sna)
library(network)
library(xtable)
library(ggplot2)
```

```{r}
w1 <- as.matrix(read.csv("./Lintner/10_W1.csv", sep = ";", row.names = 1, check.names = FALSE))
w2 <- as.matrix(read.csv("./Lintner/10_W2.csv", sep = ";", row.names = 1, check.names = FALSE))

good_ids <- Reduce(intersect, list(
  rownames(w1)[complete.cases(w1)],
  rownames(w2)[complete.cases(w2)]))
w1_clean <- w1[good_ids, good_ids]
w2_clean <- w2[good_ids, good_ids]

set.seed(161)
permutations <- 5000
nl0 <- netlogit(w2_clean, w1_clean, reps = permutations, nullhyp = "qapy")
nl0$names <- c("intercept", "wave_1")
print(summary(nl0))
```

Since the two-sided p-value is 0, under a significance level of $\alpha = 0.05$, we fail to reject the null hypothesis. Therefore, the friendship nomination in wave 2 is indeed associated with the one in wave 1.

## 2.2

We aim to add more explanatory variables in the QAP regression model we built in (1).

```{r}
attr <- read.csv("./Lintner/attr.csv", sep = ";", stringsAsFactors = FALSE)

good_ids <- Reduce(intersect, list(
  rownames(w1)[complete.cases(w1)],
  rownames(w2)[complete.cases(w2)],
  attr[
    complete.cases(attr[, c("HISEI", "literacy_end")]),
    "studentID"
  ]
))
w1_clean <- w1[good_ids, good_ids]
w2_clean <- w2[good_ids, good_ids]

# We directly extract cleaned data for classroom 10
attr10 <- attr[
  attr$classroomID == 10 &
  attr$studentID %in% good_ids,
]
```

(i) To test whether students with higher literacy scores are less likely to receive friendship nominations, we construct a dyadic matrix in which each entry $(i,j)$ is assigned the literacy score of the receiver $j$. That is, $\text{litreceiver}_{ij} = \text{literacy\_end}_j$.

    This operationalization reflects the hypothesis because friendship nominations are directed ties from sender $i$ to receiver $j$, and the hypothesis concerns characteristics of the receiver. By assigning the receiver's literacy score to each dyad, the model tests whether higher literacy decreases the likelihood of being nominated as a friend.

An alternative specification would assign the sender's literacy score, i.e. $\text{literacy\_end}_i$, to each dyad. However, this would instead test whether students with higher literacy are more likely to send friendship nominations, which corresponds to a different hypothesis. Therefore, the receiver-based operationalization is appropriate for the present analysis.

(ii) To test whether friendship nominations are more likely between students of the same gender, we construct a dyadic indicator matrix that captures gender similarity between pairs of students. Specifically, for each dyad $(i,j)$, we define:

     $$
     \text{sameGender}_{ij} =
      \begin{cases}
        1 & \text{if students } i \text{ and } j \text{ have the same gender}, \\
        0 & \text{otherwise}.
      \end{cases}
     $$

     This operationalization reflects the hypothesis because the effect of interest concerns a dyadic property—whether two individuals share the same gender—rather than an attribute of either the sender or the receiver alone. By assigning a value of $1$ to same-gender pairs and $0$ otherwise, the model directly tests whether same gender increases the likelihood of a friendship nomination.

(iii) To test whether students with higher HISEI values are more likely to send friendship nominations, we construct a dyadic matrix in which each entry $(i,j)$ is assigned the HISEI value of the sender $i$. That is, $\text{HIsender}_{ij} = \text{HISEI}_i$.

      This operationalization reflects the hypothesis because friendship nominations are directed ties from sender $i$ to receiver $j$, and the hypothesis concerns characteristics of the sender. By assigning the sender's HISEI value to each dyad, the model tests whether students with higher socioeconomic status are more likely to initiate friendship nominations.

An alternative specification would assign the receiver's HISEI value, i.e. $\text{HISEI}_j$, to each dyad. However, this would instead test whether students with higher socioeconomic status are more likely to receive nominations, which corresponds to a different hypothesis. Therefore, the sender-based operationalization is appropriate for the present analysis.

```{r}
literacy_end <- attr10$literacy_end
litreceiver <- matrix(literacy_end, length(literacy_end) , length(literacy_end) , byrow = TRUE)

gender <- attr10$gender
sameGender <- outer(gender, gender, "==") * 1

HISEI <- attr10$HISEI
HIsender <- matrix(HISEI, length(HISEI), length(HISEI), byrow = FALSE)
```

With those definitions, we can run a MR-QAP test now:

```{r}
explanatoryMatrix <- list(
  w1_clean = w1_clean,
  litreceiver = litreceiver,
  sameGender = sameGender,
  HIsender = HIsender
)


nl10 <- netlogit(w2_clean, explanatoryMatrix, reps = permutations, nullhyp = "qapspp")
nl10$names <- c("intercept","w1", "litreceiver", 
              "sameGender", "HIsender")
print(summary(nl10))
```

## 2.3

We estimate a QAP logistic regression model in which friendship nominations in wave 2 are explained by friendship ties in wave 1, gender homophily, the receiver's literacy score, and the sender's HISEI value.

The results show that friendship ties in wave 2 are strongly and significantly associated with friendship ties in wave 1. The estimated coefficient is positive and highly significant, with a one-sided (upper-tail) p-value equal to 0. This indicates strong persistence in friendship relations over time: students are much more likely to nominate as friends those whom they previously nominated.

The coefficient for $\text{sameGender}_{ij}$ is also positive and highly significant (upper-tail p-value equal to 0), providing strong evidence of gender homophily. In particular, the odds ratio (6.481) bigger than 1 suggests that friendship nominations are substantially more likely between students of the same gender.

However, the estimated effect of the sender's HISEI is not statistically significant, with an upper-tail p-value of 0.1006. Since this value exceeds the conventional significance level of 0.05, we fail to reject the null hypothesis of no association between HISEI and the likelihood of sending friendship nominations.

Although the estimated coefficient remains positive, indicating that students with higher HISEI values may be somewhat more likely to send nominations, the lack of statistical significance suggests that this effect is weak and not reliably distinguishable from random variation.

Similarly, the coefficient for $\text{litreceiver}_{ij}$ is not statistically significant. The one-sided (lower-tail) p-value is 0.9290, which is far above conventional significance levels. Therefore, we fail to reject the null hypothesis of no association between literacy and receiving nominations. Moreover, the estimated coefficient is positive, which is inconsistent with the hypothesized negative effect. Hence, there is no evidence supporting hypothesis (i) that students with higher literacy scores are less likely to receive friendship nominations.

In summary, the data provide strong support for hypotheses (ii), as well as for the persistence of friendship ties across waves, but do not support hypothesis (i) and (iii).

## 2.4

```{r}
class_ids <- c(1, 2, 3, 4, 5, 9, 10, 12)

results <- list()

for (c in class_ids) {
  cid <- sprintf("%02d", c)
  
  w1 <- as.matrix(read.csv(
    paste0("./Lintner/", cid, "_W1.csv"),
    sep = ";",
    row.names = 1,
    check.names = FALSE
  ))
  
  w2 <- as.matrix(read.csv(
    paste0("./Lintner/", cid, "_W2.csv"),
    sep=";",
    row.names = 1,
    check.names = FALSE
  ))
  
  good_ids <- Reduce(intersect, list(
    rownames(w1)[complete.cases(w1)],
    rownames(w2)[complete.cases(w2)],
    attr[
      complete.cases(attr[, c("HISEI", "literacy_end")]),
      "studentID"
    ]
  ))
  colnames(w1) <- rownames(w1)
  colnames(w2) <- rownames(w2)
  w1_clean <- w1[good_ids, good_ids]
  w2_clean <- w2[good_ids, good_ids]
  
  attr_c <- attr[
    attr$classroomID == c &
      attr$studentID %in% good_ids, ]
  
  
  litreceiver <- matrix(attr_c$literacy_end, nrow = length(good_ids), ncol = length(good_ids), byrow = TRUE)
  sameGender <- outer(attr_c$gender, attr_c$gender, "==") * 1
  HIsender <- matrix(attr_c$HISEI, nrow = length(good_ids), ncol = length(good_ids), byrow = FALSE)
  
  explanatoryMatrix <- list(
    w1_clean = w1_clean,
    litreceiver = litreceiver,
    sameGender = sameGender,
    HIsender = HIsender
  )
  
  set.seed(161) 
  model <- netlogit(w2_clean, explanatoryMatrix, reps = 5000)
  model$names <- c("intercept","w1", "litreceiver", 
                   "sameGender", "HIsender")
  results[[as.character(c)]] <- summary(model)
}
```

## 2.5

We repeat the analysis for all classrooms and compare the results across them. A clear and consistent pattern emerges for prior friendship ties, whose coefficient is positive and highly significant in all classrooms. This indicates strong persistence in friendship networks over time, suggesting that past ties are a robust predictor of future nominations. This effect appears to be highly generalizable across classrooms.

For gender homophily, the coefficient is positive in most classrooms and statistically significant in the majority of cases. This provides evidence that students tend to nominate peers of the same gender. However, the effect is not significant in all classrooms, indicating some variation in the strength of homophily across different social contexts. Thus, gender homophily appears to be a common but not universal feature.

In contrast, the effects of literacy score and HISEI are generally weak and not statistically significant across classrooms. Their coefficients are small and inconsistent in sign and significance, suggesting that these attributes do not play a systematic role in friendship formation. The lack of consistency indicates that these effects are not generalizable.

Differences across classrooms provide important insight into the robustness of the hypotheses. A mechanism that operates consistently across classrooms can be considered generalizable, whereas variation suggests context dependence, they may be sensitive to local social dynamics.

```{r}
table <- do.call(rbind, lapply(names(results), function(cid) {
  s <- results[[cid]]
  
  data.frame(
    Classroom = cid,
    Variable  = s$names,
    Estimate  = s$coefficients,
    p_value   = s$pgreqabs 
  )
}))
table <- table[table$Variable != "intercept", ]

ggplot(table, aes(x = Estimate, y = p_value)) +
  geom_point(aes(color = p_value < 0.05), size = 3) +
  geom_text(aes(label = Classroom), vjust = -1, size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = c("black", "green4"), name = "Significant (p < 0.05)") +
  labs(
    title = "Effect Sizes and Significance Across Classrooms",
    x = "Estimated Coefficient",
    y = "P-Value"
  ) +
  theme(legend.position = "top")
```

A visualization of p-values across classrooms, with a dashed line indicating the conventional significance level (0.05), provides a clear overview of which effects are consistently significant across settings.

A caterpillar plot of the estimated coefficients across classrooms further illustrates these patterns. The effect of prior friendship ties is consistently positive and large in all classrooms, indicating a robust and generalizable mechanism. The gender effect is also generally positive but shows greater variability, suggesting that homophily depends on classroom-specific dynamics. In contrast, the coefficients for literacy and HISEI fluctuate around zero, with no consistent pattern, reinforcing the conclusion that these variables do not have stable effects on tie formation.

# Task 3: Network Autocorrelation Model

## 3.1

```{r message=FALSE}
library(sna)
library(network)
library(numDeriv)
library(ggplot2)

# Read the adjacency matrix from the target classroom and do the checks
friendships_path <- paste0("./Lintner/", 10, "_W2.csv")
df <- read.csv2(friendships_path, header = TRUE, check.names = FALSE)
student_ids <- df[[1]]
friendships <- as.matrix(df[, -1])
stopifnot(identical(as.integer(colnames(friendships)), student_ids))
colnames(friendships) <- student_ids
rownames(friendships) <- student_ids
mode(friendships) <- "integer"

# Filter the attributes to contain only students from the target classroom
attr_all <- read.csv2("./Lintner/attr.csv", header = TRUE, na.strings = c("", "NA"))
attr_classroom <- attr_all[attr_all$classroomID == as.integer(10),]
stopifnot(identical(attr_classroom$studentID, student_ids))
rownames(attr_classroom) <- student_ids

# We remove students with no literacy_end, gender, or HISEI, or missing connections
ids_bad_attr <- student_ids[is.na(attr_classroom$literacy_end)
                            | is.na(attr_classroom$gender)
                            | is.na(attr_classroom$HISEI)]
ids_bad_friendships <- student_ids[apply(is.na(friendships), 1, any)]
ids_bad <- unique(c(ids_bad_attr, ids_bad_friendships))
keep <- !(student_ids %in% ids_bad)
friendships_clean <- friendships[keep, keep]
attr_clean <- attr_classroom[keep, ]
student_ids_clean <- student_ids[keep]

# Convert categorical data to numeric for the gender field and make literacy end numeric
attr_clean[attr_clean$gender == "female",]$gender <- 1
attr_clean[attr_clean$gender == "male",]$gender <- 0
attr_clean$literacy_end <- as.numeric(attr_clean$literacy_end)

# Apply the row‑stochastic normalization to the friendships
friendships_clean_row_stochastic <- sweep(
  friendships_clean,
  1,
  ifelse(rowSums(friendships_clean) == 0, 1, rowSums(friendships_clean)),
  FUN = "/"
)
```

```{r}
# Fit the model with the literacy end as the dependent variable and the friendship
# in wave 2 as the network variable
nam_model_1 <- lnam(
  y  = attr_clean$literacy_end,
  W1 = friendships_clean_row_stochastic
)
names(nam_model_1$rho1) <- "friendship"
summary(nam_model_1)
```

We observe that the network autocorrelation coefficient is $friendship = 0.852$ and highly significant with $p_{friendship} \le 0.001$. So, this initial model suggests a strong positive peer influence, namely a student's literacy score is highly positively correlated with the average literacy score of their friends.

## 3.2

```{r}
# Add gender and HISEI scores to the model as exogenous variables
covars <- model.matrix(~ gender + HISEI, data = attr_clean)
nam_model_2 <- lnam(
  y  = attr_clean$literacy_end,
  x  = covars,
  W1 = friendships_clean_row_stochastic
)
names(nam_model_2$rho1) <- "friendship"
names(nam_model_2$beta) <- c("intercept", "female", "HISEI")
summary(nam_model_2)
```

In our second model, we added the gender and HISEI scores as exogenous variables.

## 3.3

We observe that both covariates are statistically significant predictors of the ending literacy score. Being female has a strong, significant positive effect on literacy scores ($female = 20.883$ and $p_{female} \le 0.01$), and having a higher HISEI score also significantly increases literacy scores ($HISEI = 0.973$ and $p_{HISEI} \le 0.001$). What's surprising is that the network autocorrelation coefficient flips from a positive value computed at (1) to a negative value ($friendship = -0.940$) while remaining moderately significant ($p_{friendship} \le 0.05$).

## 3.4

Looking first at the variables that appear in both analyses, gender is important in both models, but in different ways. In the MR-QAP, the coefficient for $sameGender$ is positive and highly significant, indicating strong gender homophily, namely friendship nominations are much more likely between students of the same gender. In the autocorrelation model, the coefficient for $female$ is also positive and significant, indicating that females have higher ending literacy scores on average. These findings are not contradictory, because they refer to different mechanisms. More specifically, one is a dyadic similarity effect on friendship ties, while the other is a monadic effect on academic performance.

For HISEI, the two models diverge. In the MR-QAP, the sender’s HISEI has a positive but non-significant effect on sending friendship nominations, so there is no strong evidence that students with higher socioeconomic status nominate more friends. However, in the autocorrelation model, HISEI is positive and highly significant, showing that students from higher-status backgrounds tend to achieve higher literacy scores. This suggests that HISEI matters much more for academic outcomes than for friendship behavior in this classroom.

Regarding the friendship, the findings also differ across the two models. In the MR-QAP, friendship in wave 1 strongly predicts friendship in wave 2, which shows substantial stability in the network over time. By contrast, in the autocorrelation model that takes gender and HISEI into account, the network coefficient is negative and significant. This means that students tend to have literacy scores that are lower when their friends’ average literacy is higher. Importantly, this does not contradict the QAP result since tie persistence and network autocorrelation concern different aspects of the data. A stable friendship network does not necessarily imply positive similarity in literacy outcomes.

## 3.5

```{r message=FALSE}
get_classroom_model <- function(classroom_id) {
  # Read the adjacency matrix from the target classroom and do the checks
  friendships_path <- paste0("./Lintner/", classroom_id, "_W2.csv")
  df <- read.csv2(friendships_path, header = TRUE, check.names = FALSE)
  student_ids <- df[[1]]
  friendships <- as.matrix(df[, -1])
  stopifnot(identical(as.integer(colnames(friendships)), student_ids))
  colnames(friendships) <- student_ids
  rownames(friendships) <- student_ids
  mode(friendships) <- "integer"
  
  # Filter the attributes to contain only students from the target classroom
  attr_all <- read.csv2("./Lintner/attr.csv", header = TRUE, na.strings = c("", "NA"))
  attr_classroom <- attr_all[attr_all$classroomID == as.integer(classroom_id),]
  stopifnot(identical(attr_classroom$studentID, student_ids))
  rownames(attr_classroom) <- student_ids
  
  # We remove students with no literacy_end, gender, or HISEI, or missing connections
  ids_bad_attr <- student_ids[is.na(attr_classroom$literacy_end)
                              | is.na(attr_classroom$gender)
                              | is.na(attr_classroom$HISEI)]
  ids_bad_friendships <- student_ids[apply(is.na(friendships), 1, any)]
  ids_bad <- unique(c(ids_bad_attr, ids_bad_friendships))
  keep <- !(student_ids %in% ids_bad)
  friendships_clean <- friendships[keep, keep]
  attr_clean <- attr_classroom[keep, ]
  student_ids_clean <- student_ids[keep]
  
  # Convert categorical data to numeric for the gender field and make
  # literacy end numeric
  attr_clean[attr_clean$gender == "female",]$gender <- 1
  attr_clean[attr_clean$gender == "male",]$gender <- 0
  attr_clean$literacy_end <- as.numeric(attr_clean$literacy_end)
  
  # Apply the row‑stochastic normalization to the friendships
  friendships_clean_row_stochastic <- sweep(
    friendships_clean,
    1,
    ifelse(rowSums(friendships_clean) == 0, 1, rowSums(friendships_clean)),
    FUN = "/"
  )
  
  # Fit the model
  covars <- model.matrix(~ gender + HISEI, data = attr_clean)
  nam_model <- lnam(
    y  = attr_clean$literacy_end,
    x  = covars,
    W1 = friendships_clean_row_stochastic
  )
  names(nam_model$rho1) <- "friendship"
  names(nam_model$beta) <- c("intercept", "female", "HISEI")
  return (nam_model)
}
```

```{r}
classrooms <- c("01", "02", "03", "04", "05", "09", "10", "12")
all_models <- list()

for (classroom_id in classrooms) {
  all_models[[classroom_id]] = get_classroom_model(classroom_id)  
}
stopifnot(all_models[["10"]]$rho1 == nam_model_2$rho1)
stopifnot(all_models[["10"]]$beta == nam_model_2$beta)
```

## 3.6

```{r}
# Extract data from all classroom in a dataframe
results_list <- lapply(names(all_models), function(class_id) {
  model <- all_models[[class_id]]
  data.frame(
    Classroom = class_id,
    Variable = c("1. Friendship", "2. Gender", "3. HISEI"),
    Estimate = as.numeric(c(model$rho1, model$beta[2], model$beta[3])),
    SE = as.numeric(c(model$rho1.se, model$beta.se[2], model$beta.se[3]))
  )
})
results_long <- do.call(rbind, results_list)

# Calculate the p-values for all variables and classes
z_scores <- results_long$Estimate / results_long$SE
results_long$p_value <- 2 * (1 - pnorm(abs(z_scores)))

# Plot the values using a scatterplot per variable
ggplot(results_long, aes(x = Estimate, y = p_value)) +
  geom_point(aes(color = p_value < 0.05), size = 3) +
  geom_text(aes(label = Classroom), vjust = -1, size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = c("black", "green4"), name = "Significant (p < 0.05)") +
  labs(
    title = "Effect Sizes and Significance Across Classrooms",
    x = "Estimated Coefficient",
    y = "P-Value"
  ) +
  theme(legend.position = "top")
```

Based on the plotted coefficients and p-values across all eight classrooms, the network autocorrelation hypothesis is not generalizable.

First, the network peer effect (Friendship) is highly inconsistent. It is only statistically significant in classroom 10 (where the effect is strongly negative), while the other seven classrooms show no significant peer influence at all, with estimates bouncing on both sides of zero. Second, even the individual attributes (Gender and HISEI) lack universal significance. Although their effect estimates are generally positive across the plots, they only reach statistical significance in classroom 10 (Gender and HISEI) and classroom 01 (HISEI).

As a conclusion, these differences suggest that structural network effects and the influences of the selected individual attributes on literacy are heavily dependent on the context. The lack of generalizability implies that other factors need to be taken into account when considering the academic performance.
