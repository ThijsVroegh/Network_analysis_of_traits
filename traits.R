## For a more detailed description of this analysis, see the corresponding blog
## on my website >> https://tvroegh.netlify.app/ <<
## -----------------------------------------------------------------------------
library(tidyverse)
library(relaimpo)
library(matrixcalc)
library(EFAtools)
library(psychonetrics)
library(semPlot)
library(qgraph)
library(corrplot)
library(lattice)
library(networktools)
library(bnlearn)
library(Rgraphviz)
library(pcalg)
library(BayesNetBP)
require(reshape2)

set.seed(123)

## Raad factor correlation matrix and interfactor correlation matrix as sources
## for the data
## -----------------------------------------------------------------------------
L <- matrix(c(
0.075, 0.357, 0.33,
0.294, 0.316, 0.298,
-0.017, 0.488, -0.079,
-0.101, 0.79, -0.045,
0.763, 0.028, 0.17,
0.773, -0.114, 0.197,
0.76, 0.073, 0.045,
0.806, -0.05, -0.307,
0.719, 0.005, -0.018,
0.417, 0.214, -0.089,
0.687, -0.152, 0.14,
-0.044, -0.002, 0.575,
0.23, -0.065, 0.404,
0.684, -0.009, -0.067,
-0.1, -0.274, 0.382
),nrow=15,ncol=3,byrow=TRUE)

colnames(L) <- paste0("F", 1:3)
rownames(L) <- paste0("V", 1:15)

# Labels
varLabs <- c(
"Thrill and adventure seeking",
"Experience seeking",
"Boredom susceptibility", 
"Disinhibition",
"TAS",
"Peak experiences",
"Dissociated experiences",
"Openness to inner experiences", 
"Belief in supernatural",
"Emotional extraversion", 
"Intrinsic arousal", 
"Ego strength",
"Intellectual control",
"Cognitive regression", 
"Social desirability")

## Interfactor correlation matrix
## -----------------------------------------------------------------------------
F <- matrix(c(
  1, 0.225, 0.203,
  0.225, 1, -0.079,
  0.203, -0.079, 1
),nrow= 3, ncol = 3, byrow = TRUE)

colnames(F) <- rownames(F) <- paste0("F", 1:3)

## to obtain the implied underlying correlation matrix,
## multiply the matrices as follows: 
## -----------------------------------------------------------------------------
R <- L %*% F %*% t(L)
diag(R) <- 1

## -----------------------------------------------------------------------------
# Short Labels
colnames(R) <- rownames(R) <- c("AdvSe", "ExpSe", "BorSu", "Disin", "TAS","PeakE","Disso","OpeIn","BelSu","EmoEx","InAro","EgoSt","IntCo","CognR","SociD")

corrplot(R, method = "color", addCoef.col = "black",
         number.cex = .5, cl.pos = "n", diag=T,insig = "blank")

## -----------------------------------------------------------------------------
# number of factors
nfac_all <- N_FACTORS(R, N = 523, method = "ML")

nfac_all

## -----------------------------------------------------------------------------
x <- PARALLEL(R, N = 523, method = "ML")
plot(x)

## -----------------------------------------------------------------------------
spss_ml <- EFA(R, n_factors = 3, N = 523, type = "SPSS", method = "ML", rotation = "quartimin")
spss_ml

## -----------------------------------------------------------------------------
calc_matrix <- round(matrix(spss_ml$rot_loadings,
                            nrow = 15, 
                            ncol = 3)
                            ,2)

## Plot correlation matrices
par(mfrow=c(1,3))
corrplot(calc_matrix, method = "color", addCoef.col = "black",
         number.cex = .9, cl.pos = "n", diag=T, insig = "blank")
title(main = "Estimated matrix")
corrplot(abs(L-calc_matrix), method = "color", addCoef.col = "black",
         number.cex = .9, cl.pos = "n", diag=T, insig = "blank")
title(main="Absolute Difference")
corrplot(L, method = "color", addCoef.col = "black",
         number.cex = .9, cl.pos = "n", diag=T, insig = "blank")
title(main="Original factor loadings")

## -----------------------------------------------------------------------------
is.positive.definite(round(R, 2))

## -----------------------------------------------------------------------------
minvalue <- 50
minseed <- 9000

for(i in 35000:40000) {
  set.seed(i)
  M <- chol(R)
  r <- t(M) %*% matrix(rnorm(15*523), nrow = 15, ncol = 523)
  r <- t(r)
  
  rdata <-  as.data.frame(r)
  corrdata <- round(cor(rdata),3)
  
  diff <- abs(R-corrdata)
  diff2 <- diff > 0.05
  
  sumdiff <- sum(diff2)
  
  # updating minimal value
  minvalue <- ifelse ((sumdiff < minvalue), sumdiff, minvalue) 
  minseed <- ifelse ((sumdiff <= minvalue), i, minseed) 
}

## -----------------------------------------------------------------------------
minvalue
minseed

## -----------------------------------------------------------------------------
set.seed(minseed)

M <- chol(R)
nvars <-dim(M)[1]

# number of observations to simulate
nobs <- 523

# Random variables that follow the R correlation matrix
r <- t(M) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
r <- t(r)
rdata <-  as.data.frame(r)
rrdata <- round(cor(rdata),3)

diff <- abs(R-rrdata)
diff2 <- diff > 0.05
sum(diff2)

## -----------------------------------------------------------------------------
corrplot(abs(R-rrdata), method = "color", addCoef.col = "black",
         number.cex = .6, cl.pos = "n", diag=T, insig = "blank")
title(main="Absolute Difference")

## test for overlapping/ redundant variables
## -----------------------------------------------------------------------------
badpairs <- goldbricker(rdata, p = 0.05, method = "hittner2003", threshold = 0.20, corMin = 0.5, progressbar = FALSE)
badpairs

## -----------------------------------------------------------------------------
R13 <- R %>% as.data.frame() %>% 
            dplyr::select(-PeakE,-Disso) %>%
            slice(-6,-7) %>% 
            as.matrix()

## -----------------------------------------------------------------------------
graph <- EBICglasso(R13, n = 523, gamma = 0.5)
g <- qgraph(graph, layout = "spring", title = "GGM with EBICglasso", details = TRUE)

## -----------------------------------------------------------------------------
# plot centrality plot
centralityPlot(EBIC = graph, scale = "z-scores", include =c("ExpectedInfluence"))


## -----------------------------------------------------------------------------
## Adapted from code on relative importance network from Oisín Ryan et al., 
## as part of article (see https://github.com/ryanoisin/SEset)

rownames(R13) <- colnames(R13)
names.sub <- colnames(R13)

# relmat will contain the weights matrix of the relative importance network
relmat <- matrix(0, 13, 13)

# r2vec will contain the predictability scores
r2vec <- numeric(13)

# function needed to reorder the matrix according to "names" 
reorder_mat <- function(matrix, names) { 
  if (is.null(dimnames(matrix))) { 
    stop("Error: matrix must have dimension names") 
  } 
  if (!all(rownames(matrix) %in% names)) { 
    stop("Error: dimnames(matrix) does not match names") 
  } 
  matrix[names,names] 
} 

## -----------------------------------------------------------------------------
for(i in 1:13) 
  {
  # reorder variables so each, in turn, is first in the matrix
  temp <- reorder_mat(R13, c(names.sub[i], names.sub[seq(1:13)[-i]]))
  
  # then estimate relative importance for each variable
  relmat[seq(1:13)[-i], i] <- calc.relimp(temp,type="lmg",rela=TRUE)$lmg
  
  # and the predictability of each variable
  r2vec[i] <- calc.relimp(temp)$R2
  }
relimpnet_c <- as.table(apply(relmat, 2, as.numeric))

# for clarity purposes, we delete all edges below 0.05
relimpnet_censored <- ifelse(relimpnet_c < 0.05 ,0, relimpnet_c)


## -----------------------------------------------------------------------------
relimpnet_plot <- qgraph(relimpnet_censored, labels = names.sub, title= "censored graph", repulsion = .8, asize = 5, pie = r2vec)


## -----------------------------------------------------------------------------
centralityPlot(relimpnet_plot)

## -----------------------------------------------------------------------------
rdata2 <- rdata %>% dplyr::select(-PeakE,-Disso)

## -----------------------------------------------------------------------------
set.seed(38022)

bootnet <- boot.strength(rdata2,
                         R = 100, 
                         algorithm = "hc", 
                         debug = FALSE,
                         algorithm.args = list(restart = 5, perturb = 10))


## -----------------------------------------------------------------------------
avgnet <- averaged.network(bootnet)
avgnet
avgnet$learning

#Compute the score of the Bayesian network
score(avgnet,data = rdata2) # -9025

thresh <- avgnet$learning$args[[1]]
thresh  # optimal significance threshold = 0.5

astr <- arc.strength(avgnet, rdata2, "bic-g")   ## compute edge strengths
astr

nrow(bootnet[with(bootnet, strength > 0.51 & direction > 0.50), ])


## -----------------------------------------------------------------------------
strength.plot(avgnet, astr, main = "Averaged network", shape = "ellipse",threshold = thresh,layout ="dot",
highlight = list(nodes = c("Disin","EgoSt","TAS"),
col = "tomato", fill = "orange"))


## -----------------------------------------------------------------------------
centr_bn <- qgraph(avgnet,DoNotPlot=T)
centralityPlot(centr_bn, theme_bw = T, scale = "z-scores") 


## -----------------------------------------------------------------------------
fitted.bn <- bnlearn::bn.fit(avgnet, rdata2, method = 'mle')
#print(fitted.bn)


## -----------------------------------------------------------------------------
rank <- bootnet[bootnet$from == 'TAS' & bootnet$strength > 0.5, c('from', 'to', 'strength')]
print(rank[order(-rank$strength), ])


## -----------------------------------------------------------------------------
rank <- bootnet[bootnet$from == 'TAS' & bootnet$direction > 0.5, c('from', 'to', 'direction')]
print(rank[order(-rank$direction), ])


## -----------------------------------------------------------------------------
predcor <-structure(numeric(13), names = c("AdvSe", "ExpSe", "BorSu", "Disin","TAS","OpeIn","BelSu","EmoEx","InAro","EgoSt","IntCo","CognR","SociD"))
                    
for (var in names((predcor))) {
  xval = bn.cv(avgnet, data=rdata2, loss = "cor", method = "k-fold",
               loss.args = list(target = var), runs = 10)
  predcor[var] = mean(sapply(xval, function(x) attr(x, "mean")))
  }

round(predcor, digits = 3)

summary(predcor)


## -----------------------------------------------------------------------------
boottab <- bootnet[bootnet$strength > thresh & bootnet$direction > 0.5, ]  ## edges in avgnet
nrow(boottab)

astr3 <- boottab    ## table with direction probabilities
astr3$strength <- astr3$direction  ## use the direction probabilities for edge width


## -----------------------------------------------------------------------------
# Histogram of directionality as part of the avgnet network 
ggplot(data = astr3, aes(x = paste0(from,"->",to), y = direction)) +   
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = 0.50, linetype = "dashed", color = "red") +
  
    coord_flip() +
   
    theme(axis.text.x  = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
  
  theme_minimal() +

  labs(title= "direction probabilities of edges in network")


## -----------------------------------------------------------------------------
par(mfrow = c(1,1))
cpdag <- cpdag(avgnet)
nrow(undirected.arcs(cpdag(cpdag)))
undirected.arcs(cpdag)


## -----------------------------------------------------------------------------
graphviz.plot(cpdag,
              highlight = list(
              nodes = descendants(cpdag, "TAS"), col = "tomato", fill = "orange",
              arcs = undirected.arcs(cpdag), lwd=3))


## -----------------------------------------------------------------------------
isValidGraph(t(amat(avgnet)), type = "dag",verbose=TRUE)  
isValidGraph(t(amat(cpdag)), type = "cpdag",verbose=TRUE) 
isValidGraph(t(amat(cpdag)), type = "pdag",verbose=TRUE)


## -----------------------------------------------------------------------------
plotbn <- bn_to_graphNEL(avgnet)

# all variables are continuous
node.class <-  c("AdvSe"=FALSE, "ExpSe"=FALSE, "BorSu"=FALSE,
                 "Disin"=FALSE, "TAS"  =FALSE, "OpeIn"=FALSE,
                 "BelSu"=FALSE, "EmoEx"=FALSE, "InAro"=FALSE,
                 "EgoSt"=FALSE, "IntCo"=FALSE, "CognR"=FALSE, "SociD"=FALSE)

tree.init.p <- Initializer(dag=plotbn, data=rdata2,
                           node.class = node.class,
                           propagate = TRUE)

tree.init.p@propagated
tree.post <- AbsorbEvidence(tree.init.p, vars=c("TAS"), values=list(3)) # #remember, we work with standardized data 

#check variables and values
tree.post@absorbed.variables
tree.post@absorbed.values

div <- PlotCGBN(tree.init.p, tree.post, fontsize = 24)

## -----------------------------------------------------------------------------
varm <- "EmoEx"
marg.2 <- Marginals(tree.post, varm)
marg.1 <- Marginals(tree.init.p, varm)
mgns <- list(marg.1$marginals[[1]], marg.2$marginals[[1]])
names(mgns) <- c(varm, varm)
types <- c(marg.1$types[1], marg.2$types[1])
marg <- list(marginals=mgns, types=types)
SummaryMarginals(marg.1)
##              Mean        SD n
## EmoEx -0.06827712 0.9850553 1
SummaryMarginals(marg.2)
##           Mean        SD n
## EmoEx 1.012537 0.9161935 1
PlotMarginals(marg, groups=c("before", "after"))

## -----------------------------------------------------------------------------
detach(package:bnlearn, unload=TRUE)

## code adapted from Kan, de Jonge, van der Maas, Levine & Epskamp (2020).
## see https://github.com/KJKan/mcfarland/blob/master/Tutorial.md

yvars    <- colnames(R13)
ny       <- length(yvars)
n_sample <- 523

# latent constructs to be measured (etas)
latents  <- c("Absorption permissiveness",
              "Sensation seeking",
              "Social desirability")

# short names
lvars <- c(
"AP", # absorption permissiveness
"SS", # sensation seeking
"SD"  # Social desirability
)

ne <- length(lvars)

## -----------------------------------------------------------------------------
lambda <- matrix( c (
#P V W
0, 1, 0, # Thrill and adventure seeking
0, 1, 0, # Experience seeking
0, 1, 0, # Boredom susceptibility
0, 1, 0, # Disinhibition
1, 0, 0, # TAS
#1, 0, 0, # Peak experiences
#1, 0, 0, # Dissociated experience
1, 0, 0, # Openness to inner experiences
1, 0, 0, # Belief in supernatural
1, 0, 0, # Emotional extraversion
1, 0, 0, # Intrinsic arousal
0, 0, 1, # Ego strength
0, 1, 1, # Intellectual control
1, 0, 0, # Cognitive regression
0, 0, 1  # Social desirability
),
ncol = ne,
byrow = TRUE,
dimnames = list( yvars, lvars )
)


## -----------------------------------------------------------------------------
lambda_measurement <- lambda
lambda_measurement[ 1, 3 ] <- 1 
lambda_measurement[ 2, 1 ] <- 1 
lambda_measurement[ 2, 3 ] <- 1 
lambda_measurement[ 6, 3 ] <- 1
lambda_measurement[ 8, 2 ] <- 1
lambda_measurement[ 13, 2 ] <- 1

lambda_measurement

## -----------------------------------------------------------------------------
measurementModel <- lvm( covs = ( n_sample - 1 )/n_sample * R13,
    lambda = lambda_measurement,
    nobs = n_sample,
    identification = "variance",
    latents = latents)

measurementModel <- measurementModel %>% runmodel 
measurementModel %>% fit
#measurementModel %>% parameters
#measurementModel %>% MIs

## -----------------------------------------------------------------------------
results_measurementModel <- measurementModel %>% runmodel


## -----------------------------------------------------------------------------
lambda_g               <- cbind( lambda_measurement, g = 0 )
beta_g                 <- cbind( matrix( 0, ne + 1, ne + 1 ) ) 
beta_g[ 1:ne, ne + 1 ] <- 1

gModel    <- lvm( covs = ( n_sample - 1 )/n_sample * R13, 
                  lambda = lambda_g, 
                  beta = beta_g,
                  sigma_zeta = 'empty',
                  nobs = n_sample,
                  identification = "variance" )

results_gModel <- gModel %>% runmodel


## -----------------------------------------------------------------------------
lambda_bifactor    <- cbind( lambda, g = 1 )

bifactorModel    <- lvm( covs = ( n_sample - 1 )/n_sample *R13, 
                         lambda = lambda_bifactor, 
                         sigma_zeta = 'empty',
                         nobs = n_sample,
                         identification = "variance" )

results_bifactorModel <- bifactorModel %>% runmodel

## -----------------------------------------------------------------------------
saturatedModel <- ggm( covs = ( n_sample - 1 )/n_sample*R13,
                          omega = "Full",
                          nobs = n_sample )

# Use 'runmodel' to compute parameters and fit measures
saturatedModel <- saturatedModel %>% runmodel

#This initial model is saturated (df=0)
#saturatedModel

# Check parameters
#saturatedModel %>% parameters

## -----------------------------------------------------------------------------
prunedModel <- saturatedModel %>%
               prune( alpha = 0.01, 
                      recursive = TRUE, 
                      adjust = "fdr")

# Check modification indices
prunedModel %>% MIs

## -----------------------------------------------------------------------------
# Stepup estimation
prunedModel_stepup <- prunedModel %>% 
                      stepup()
prunedModel_stepup %>% fit # Model fit

## -----------------------------------------------------------------------------
# Modelsearch
prunedModel_modelsearch  <- prunedModel %>% 
      modelsearch()
prunedModel_modelsearch %>% fit # Model fit

## -----------------------------------------------------------------------------
compare(saturated = saturatedModel, 
           pruned = prunedModel,
           stepup = prunedModel_stepup,
           search = prunedModel_modelsearch)

## -----------------------------------------------------------------------------
net1 <- getmatrix(prunedModel_modelsearch, "omega")

rownames(net1) <- colnames(net1) <- colnames(rdata2)

qgraph(net1, 
       layout = "spring", 
       title = "Pruned model with modelsearch",
       palette = "colorblind",
       groups = list( "Absorption permissiveness" = which( lambda[ , 1 ] == 1 ),
                      "Sensation seeking"         = which( lambda[ , 2 ] == 1 ),
                      "Social desirability"       = which( lambda[ , 3 ] == 1 ) 
                    ))

## -----------------------------------------------------------------------------
compare( saturated       = saturatedModel,
         ggm_modelsearch = prunedModel_modelsearch,
         measurement     = results_measurementModel,
         bifactor_model  = results_bifactorModel,
         gmodel          = results_gModel)
