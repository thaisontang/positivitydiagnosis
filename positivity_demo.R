# Libraries
library(nnet); library(tidyverse)
library(tableone); library(survey)

# Generate Toy Dataset
cat.levels <- LETTERS[1:5]
n <- 464
num.expos <- 5

logit <- function(p) {log(p/(1-p))}
expit <- function(q) {1/(1 + exp(-q))}
loglog <- function(p) {log(-log(p))}
expexp <- function(q) {exp(-exp(q))}

# Simulating Patient Covariates
x1 <- rnorm(n)                                                    
x2 <- factor(sample(cat.levels, n,                            
                    replace = T, prob = c(0.35, 0.25, 0.2, 0.1, 0.1)))   
x2 <- relevel(x2, ref = cat.levels[1])

# Simulating Exposure Assignment Probabilities
cat.assign <- function(num) { num:1/(sum(num:1)) }
beta.cat.assign <- function(num) {
  sets <- 0:(num - 1)
  step <- 1.25
  t(sapply(sets, function(x) { step*seq(x, -x, length.out = 4) }))
}
z.param = list(`int` = c(1, 3, 4, 3)/1.75, 
               `beta1` = c(2.25, 2, 1, 6),
               `beta2` = beta.cat.assign(5)
)
y.param = list(`beta1` = 2,
               `beta2` = seq(0, 2, length.out = 5),
               `z` =  c(1.5, 1, 0.5, -0.5),
               `global.int` = 0
)
eta.z <- matrix(z.param$`int`, n, num.expos-1, byrow = TRUE) + 
                   matrix(z.param$`beta1`, n, num.expos-1, byrow = TRUE) * matrix(x1, n, num.expos-1, byrow = FALSE) +
                   z.param$`beta2`[as.numeric(x2),]

prob.z <- exp(eta.z)
prob.z <- cbind( prob.z/(1 + matrix(rowSums(prob.z), n, num.expos-1, byrow = FALSE)), 1/(1 + rowSums(prob.z)))
# Population Average Assignment Probability
colMeans(prob.z)

# Simulating Hospital Assignment
z <- sapply(1:n, function(x) { 
  sample(1:num.expos, 1, prob = prob.z[x,])
} )


# Global Intercept Term
eta <- matrix(y.param$global.int, n, num.expos, byrow = T)

# Hospital Effect
eta <- eta + (matrix(y.param$z, n, num.expos, byrow = T))

# # Hospital-Covariate Interactions (with Continuous Covariate)
# eta <- eta + ind.hosp.int * matrix(y.param.int$`hosp.cov`, n, num.expos, byrow = T) * x2

# Full Covariate Effect
eta <- eta +
    matrix(y.param$beta1, nrow = n, ncol = num.expos, byrow = T) * x1 +
      y.param$`beta2`[as.numeric(x2)]


# Covariate-Covariate Interactions  
# eta <- eta + ind.cov.int * matrix(y.param.int$`cov.cov`, n, num.expos, byrow = F) * x1 * x2

prob.y <- expit(eta)

# Potential Binary Outcomes
yi.all <- matrix(  rbinom(n*num.expos, 1, as.numeric(prob.y)),
                   nrow = n, ncol = num.expos, byrow = FALSE)

# Observed Surival Outcome
yi <- yi.all[cbind(1:n,z)]

# Observed Data
dataset <- data.frame(`ID` = 1:n,
                      `Outcome` = yi,
                      `Drug` = factor(paste0("Drug ", z)),
                      `Continuous` = x1,
                      `Categorical` = x2)

# Exposure Model (Propensity Score)
nnet.mod <- multinom(Drug ~ Continuous + Categorical, data = dataset)

# Assignment Probabilities (Predicted and Observed)
ps <- fitted(nnet.mod)

# IP Weights
ipw <- 1/ps[cbind(1:n, dataset$Drug)]
qs <- quantile(ipw, prob = 0.99)
ipw.trunc <- ifelse(ipw > qs, qs, ipw)
plot(ipw, ipw.trunc)
  
# Overlap Boxplots of Assignment Probabilities for Patients under different Drug Plans being Assigned Drug 3
dx <- 3
boxplot(ps[,dx] ~ dataset$Drug, horizontal=TRUE, xlab=bquote('P(' ~ Z[i]==.(i) ~ '|' ~ x[i] ~ ')'), 
        cex.axis=0.75, cex.lab=0.75, ylab = bquote("Drug ID"))
abline(v=c(0.0,0.5,1.0), lty='dotted')

# All Drug Plans
levels <- levels(dataset$Drug)
op <- par(las=1, mar=c(3.5,3.5,0.1,0.1), oma=0.2 * c(1,1,1,1), mfrow=c(2,3), mgp=c(2,0.75,0))
for (dx in 1:length(levels)) {
  boxplot(ps[,dx] ~ dataset$Drug, horizontal=TRUE, xlab=bquote('P(' ~ Z[i]==.(dx) ~ '|' ~ x[i] ~ ')'), 
          cex=0.5, cex.axis=0.5, cex.lab=0.75)
  abline(v=c(0.0,0.5,1.0), lty='dotted')
  
}
par(op)

# Bhattacharrya Coefficients
grid <- seq(0,1,by=0.01)
bhatta <- matrix(NA, length(levels), length(levels))
rownames(bhatta) <- 1:length(levels)
colnames(bhatta) <- 1:length(levels)
for (i in 1:length(levels)) {
  for (j in 1:length(levels)) {
    index <- ps[dataset$Drug == levels[i],levels[i]]
    ref <- ps[dataset$Drug == levels[j],levels[i]]
    pindex <- rep(0.0, length(grid))
    pref <- rep(0.0, length(grid))
    tab <- table(findInterval(index, grid))/length(index)
    pindex[as.numeric(names(tab))] <- tab
    tab <- table(findInterval(ref, grid))/length(ref)
    pref[as.numeric(names(tab))] <- tab
    bhatta[i,j] <- sum(sqrt(pindex * pref))
  }
}

b.dist <- bhatta            # Bhatt Coefficient
b.dist.values <- as.tibble(b.dist) %>% gather(`Drug 1`,  `B Dist`, factor_key = T) %>%
  mutate(`Drug 2` = factor(rep(1:length(levels), length = ncol(ps)^2)))

plot.sig <- ggplot(b.dist.values, aes(`Drug 1`, `Drug 2`)) +
  geom_tile(aes(fill = `B Dist`), colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "white") + theme_bw() +
  ggtitle("Bhattacharyya Coefficients between Distributions of\nIndex Drugs & Reference Drugs") +
  labs(x = "Reference Drug", y = "Index Drug") +
  theme(legend.title = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 11, margin = margin(t = 18)),
        axis.title.y = element_text(size = 11, margin = margin(r = 18)),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11,
                                  margin = ggplot2::margin(b = 12)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_colourbar(
    frame.colour = "black", label.vjust = 1, title = "Bhattacharyya\nCoefficient\n")) 
print(plot.sig)


# Standardized Mean Differences
(smd.table <- CreateTableOne(vars = c("Continuous", "Categorical"), data = dataset, strata = "Drug",
                            factorVars = "Categorical"))
smd.values <- as.tibble(cbind(`Comparison` = rownames(t(ExtractSmd(smd.table))[-1,]), t(ExtractSmd(smd.table))[-1,])) %>%
  mutate(`Drug 1` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[1])),
         `Drug 2` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[3]))) %>%
  gather(Variable, SMD, Continuous:Categorical, factor_key = T) %>%
  rename(`Covariates` = Variable) %>%
  mutate(`SMD` = as.numeric(SMD))


# Weighted Dataset
wgt.data <- svydesign(ids = ~ 1, data = dataset, weights = ~ ipw.trunc)
(smd.wgt.table <- svyCreateTableOne(vars = c("Continuous", "Categorical"), data = wgt.data, strata = "Drug",
                                    factorVars = "Categorical"))


wgt.smd.values <- as.tibble(cbind(`Comparison` = rownames(t(ExtractSmd(smd.wgt.table))[-1,]),
                                  t(ExtractSmd(smd.wgt.table))[-1,])) %>%
  mutate(`Drug 1` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[1])),
         `Drug 2` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[3]))) %>%
  gather(Variable, SMD, `Continuous`:`Categorical`, factor_key = T) %>%
  rename(`Covariates` = Variable) %>%
  mutate(`SMD` = as.numeric(SMD))



# Combined SMD Results
smd.all.a <- smd.values %>% mutate(`State` = "Pre",
                                   `Category` = paste0("Pre-", `Covariates`))
smd.all.b <- wgt.smd.values %>% mutate(`State` = "Post",
                                       `Category` = paste0("Pre-", `Covariates`))
smd.all <- as.tibble(rbind(smd.all.a, smd.all.b))
smd.all$State <- factor(smd.all$State)
smd.all$State <- relevel(smd.all$State, ref = "Pre")
smd.all$Category <- factor(smd.all$Category)


# Plot 1: Before and After SMD Boxplot

smd.both.vplot <- function(data, title = NULL) {
  w = 13.5
  h = 9.75
  
  v.plot <- ggplot(data, aes(x = `Covariates`, y = SMD)) +
    geom_boxplot(aes(fill = `State`)) +
    
    coord_flip() + theme_bw() + scale_x_discrete(limits = rev(levels(data$`Covariates`))) +
    labs(x = "Covariates", y = "Pairwise Standardized Mean Difference") +
    scale_fill_manual("Weighting", labels = c("Unweighted", "IPW"), values = c("#56B4E9", "#F0E442")) +
    theme(axis.text.x = element_text(size = 10, margin = margin(t = 10)),
          axis.text.y = element_text(size = 10, margin = margin(r = 10)),
          axis.title.x = element_text(size = 12, margin = margin(t = 24)),
          axis.title.y = element_text(size = 12, margin = margin(r = 24)),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 24)),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10),
          
          # Specifies Top Right Corner (1,1 is coordinate for proportionate location on graph)
          legend.position = c(1,1),            
          legend.justification = c(1,1),
          
          # White Fill, Black Border
          legend.background = element_rect(color = "black")) +
    ggtitle(title) +
    geom_hline(linetype = 2, col = "black", yintercept = 0.2) +
    geom_hline(linetype = 2, col = "black", yintercept = 0.5) +
    geom_hline(linetype = 2, col = "black", yintercept = 0.8)
  print(v.plot)
}


smd.both.vplot(smd.all, title = "Pairwise Standardized Mean Difference Distribution
               of Covariates Before and After IP Weighting")
