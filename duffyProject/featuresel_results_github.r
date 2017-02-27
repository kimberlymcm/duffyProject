# Kimberly McManus
# September 2015
# Scripts for ABC analysis of Duffy locus
# Step 1: model selection for initial frequency via ABC (Approximate Bayesian Computation)
# Step 2: Infers posterior for selection coefficient

# Set current directory
setwd("/Users/kimberlymcmanus/Documents/Duffy/scripts/ABC")
source('featuresel_utils_github.r')

# Initial Frequency Model Selection Step

# Read in simulation files of 5 initial frequencies
de.novo <- readInFile('data/simon_result_20150915_LWK_0_1mill.uniq.txt.200000')
point.one.percent <- readInFile('data/simon_result_20150923_LWK_0.001_200k_fixed.txt')
one.percent <- readInFile('data/simon_result_20150915_0.01_1mill_uniq.txt.200000')
ten.percent <- readInFile('data/simon_result_20150915_LWK_0.1_1mill.uniq.txt.200000')
twenty.five.percent <- readInFile('data/simon_result_20150915_LWK_0.25_1mill.uniq.txt.200000')

# Combine data
max <- 100000 # Ensure all models have the same number of simulations
de.novo <- de.novo[1:max,]
point.one.percent <- point.one.percent[1:max,]
one.percent <- one.percent[1:max,]
ten.percent <- ten.percent[1:max,]
twenty.five.percent <- twenty.five.percent[1:max,]
all.simulations <- as.data.frame(rbind(de.novo, point.one.percent, one.percent, 
                                       ten.percent, twenty.five.percent))
# Create labels for the data frame of all simulations
model.tags <- c(rep("de.novo",NROW(de.novo)), 
                rep("point.one.percent", NROW(point.one.percent)), 
                rep("one.percent", NROW(one.percent)), 
                rep("ten.percent", NROW(ten.percent)), 
                rep("twenty.five.percent", NROW(twenty.five.percent)))
model.tags.factors <- factor(model.tags, c("de.novo", "point.one.percent", "one.percent", 
                                           "ten.percent", "twenty.five.percent"))

# Make sure this package isn't loaded because it clashes with caret
detach("package:mixOmics", unload=TRUE) 
library('caret')
normalize.data <- preProcess(all.simulations[,2:NCOL(all.simulations)], 
                             method=c("center", "scale")) 
normalize.sims <- predict(normalize.data, all.simulations[,2:NCOL(all.simulations)])
# Convert the real stats
normalize.real.stats <- predict(normalize.data, t(as.data.frame(real.stats))) 
names(normalize.real.stats) <- names(normalize.sims)

# Run plsda on features
detach("package:caret", unload=TRUE)
library('mixOmics')
transformed.sims <- plsda(X=normalize.sims, Y=model.tags, ncomp=5) 
real.transformed.data <- predict(object=transformed.sims, newdata=normalize.real.stats)
# The 5 components from PLSDA
colnames(real.transformed.data$variates) <-  c('V1', 'V2', 'V3', 'V4', 'V5') 

library('abc')
library('grDevices')
transformed.sims <- as.data.frame(transformed.sims$variates$X)
colnames(transformed.sims) <-  c('V1', 'V2', 'V3', 'V4', 'V5')

# Model Selection
model.selection.result <- postpr(target=real.transformed.data$variates, 
                                 index=model.tags, sumstat=transformed.sims, 
                                 tol=0.01, method="mnlogistic")
summary(model.selection.result)


# Step 2: Inference of selection coefficient from best fit model from step 1

# Read in data
fileName <- 'data/simon_result_20150915_LWK_0.001_1mill_uniq.txt'
data <- readInFile(fileName)

# Summary statistics we are using for this part (chosen via information gain)
stats.to.use <- c("seg", "fixed", "X2", "H1", "H2", "H12", "uniq", "ehh.avg", 
                  "finalFreq", "sf")

# Make sure this package isn't loaded because it clashes with caret
detach("package:mixOmics", unload=TRUE) 
library('caret')
transformation <- preProcess(data[ ,stats.to.use], method=c("center", "scale", "pca"))
sims.pca <- predict(transformation, data[ ,stats.to.use])
real.data.pca <- predict(transformation, rbind(real.stats[stats.to.use], 
                                               real.stats[stats.to.use]))

# Get ABC result
library('abc')
result <- abc(target=real.data.pca[1,], param=log10(data$s), sumstat=sims.pca, 
           tol=0.01, method="loclinear")

# Plot prior & posterior densities
graph.labels = data.frame(value=c(result$adj.values, log10(data$s)), 
                          m=c(rep("Posterior", length(res$adj.values)),
                              rep("Prior", length(data$s))))
result.plot <- ggplot(graph.labels, aes(x=value)) 
result.plot + geom_density(aes(fill=m), size=1, alpha=0.3) + 
  scale_fill_manual(values=c("#0072B2","#D55E00")) +
  theme(text=element_text(size=30), legend.title=element_blank()) + 
  labs(x="Selection Coefficient", y="Density") +
  scale_x_continuous(breaks=c(-3,-2,-1), 
                     labels=c(parse(text='10^-3'), parse(text='10^-2'), 
                                    parse(text='10^-1')))