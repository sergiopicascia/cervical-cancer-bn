# SERGIO PICASCIA 943865 - Probabilistic Modeling

library(dplyr)
library(bnlearn)
library(Rgraphviz)
library(gRain)
library(gRbase)

### DATA MANIPULATION ###

# Import data
path <- 'risk_factors_cervical_cancer.csv'
data <- read.csv(path, na='?')
summary(data)

# Remove redundant columns and ones with all zeros 
data <- subset(data, select = -c(STDs.cervical.condylomatosis, STDs.AIDS, STDs..Number.of.diagnosis, 
                                 STDs..Time.since.last.diagnosis, Smokes..years., Smokes..packs.year.,
                                 Hormonal.Contraceptives..years., IUD..years., STDs..Time.since.first.diagnosis,
                                 STDs..number., STDs, STDs.condylomatosis, Dx, STDs.HPV))

# New feature 'Cervical.Cancer': 1 if at least one test is positive, 0 otherwise
data$Cervical.Cancer <- with(data, ifelse((Hinselmann+Schiller+Citology+Biopsy) >= 1, 1, 0))

# Remove rows with NAs
df <- na.omit(data)

# Plotting continuous vars
hist(data$Age)
hist(data$Number.of.sexual.partners)
hist(data$First.sexual.intercourse)
hist(data$Num.of.pregnancies)
pairs(data[, 0:4])

# Discretize continuous vars
df <- df %>% mutate(Age = case_when(Age < 20 ~ '<20',
                                    Age >= 20 & Age < 30 ~ '20-29',
                                    Age >= 30 ~ '30+'),
                    Number.of.sexual.partners = case_when(Number.of.sexual.partners == 1 ~ '1',
                                                          Number.of.sexual.partners == 2 ~ '2',
                                                          Number.of.sexual.partners >= 3 ~ '3+'),
                    First.sexual.intercourse = case_when(First.sexual.intercourse < 16 ~ '<16',
                                                         First.sexual.intercourse >= 16 & 
                                                          First.sexual.intercourse < 19 ~ '16-18',
                                                         First.sexual.intercourse >= 19 ~ '19+'),
                    Num.of.pregnancies = case_when(Num.of.pregnancies == 0 ~ '0',
                                                   Num.of.pregnancies == 1 ~ '1',
                                                   Num.of.pregnancies == 2 ~ '2',
                                                   Num.of.pregnancies >= 3 ~ '3+'))

# Convert variables to factor
df[colnames(df)] <- lapply(df[colnames(df)], factor)

### BAYESIAN NETWORK ###
df <- subset(df, select = -c(Hinselmann, Schiller, Citology, Biopsy)) # Considering only Cervical.Cancer

# Blacklist and whitelist
cols <- colnames(df)
bl1 <- data.frame(from = cols[-grep('Age', cols)], # Prevent parents of Age
                  to   = c('Age'))
bl2 <- data.frame(from = c('Cervical.Cancer'), # Prevent children of Cervical.Cancer
                  to   = cols[-grep('Cervical.Cancer', cols)])

bl <- rbind(bl1, bl2)

wl <- data.frame(from = c('Dx.HPV', 'STDs.HIV', 'Smokes', 'Hormonal.Contraceptives', 'Num.of.pregnancies',
                          'Number.of.sexual.partners', 'STDs.genital.herpes', 'Age', 'IUD', 'First.sexual.intercourse',
                          'Dx.CIN'),
                 to   = c('Cervical.Cancer'))

# Run structure learning algorithms and sum the adjacency matrices
sl_algos <- c(pc.stable, gs, iamb, hc, tabu, rsmax2, mmhc)
models <- list()
adj_mat <- matrix(0L, nrow = 19, ncol = 19)

for (a in sl_algos) {
 model = a(df, blacklist = bl, whitelist = wl)
 models <- append(models, list(model))
 adj_mat <- adj_mat + amat(model)
 graphviz.plot(model, shape = 'ellipse', layout = 'fdp')
}

# Retrieve the most frequent edges
adj_mat[adj_mat < 2] <- 0L
adj_mat[adj_mat >= 2] <- 1L

# Build the BN
model <- empty.graph(cols)
amat(model) <- adj_mat
model <- pdag2dag(model, ordering = cols)
model

# Parameter learning
fitted_model <- bn.fit(model, df, method = 'bayes')

# Plot of the graph
graphviz.plot(model, shape = 'ellipse', layout = 'fdp', 
              highlight = list(nodes = 'Cervical.Cancer', arcs = incoming.arcs(model, 'Cervical.Cancer'),
                               lty = 5, fill = 'darkturquoise', col = 'turquoise4'))

# Plotting marginal probabilities
graphviz.chart(fitted_model, layout = 'fdp', type = 'barprob', scale = c(2, 2), bar.col = "darkturquoise",
               strip.bg = "darkturquoise")

# Conditional probability distributions
for (col in cols) {
 distr <- cpdist(fitted_model, col, (Cervical.Cancer == '1'))
 n <- nrow(distr)
 print(col)
 print(table(distr)/n)
}

# Maximum a posteriori query
gr.fit <- as.grain(fitted_model)
cerv.canc <- setEvidence(object = gr.fit, nodes = 'Cervical.Cancer', states = '1')
joint.post <- querygrain(cerv.canc, type = 'joint')

map <- function(joint) {
 ind_max <- which(sapply(joint, function(v) isTRUE(all.equal(max(joint), v))))
 ind <- arrayInd(ind_max, .dim = dim(joint))
 state <- mapply('[', dimnames(joint), ind)
 prob <- joint[ind_max]
 list(state=state, prob=prob)
}

map(joint.post)