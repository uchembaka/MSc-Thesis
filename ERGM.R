########## Libraries ##########

library(corrplot)
library(igraph)
library(statnet)
library(ergm.count)
library(intergraph)
library(RColorBrewer)
library(kableExtra)
library(GGally)
library(ggnet)
library(stargazer)
library(coda)

###############################





################### Load Data sets ####################

load("Net_M_Feb.RData")
load("Net_E_Feb.RData")
cols <- brewer.pal(8, "Set2")
cols.2 <- brewer.pal(9, "Set1")
Net.M <- Net_M_Feb
Net.E <- Net_E_Feb

######################################################


########################## Morning Network ##################################

##### ----- complete network

### Plot
net.1 <- Net.M$weighted_net
plot(asIgraph(net.1), vertex.color = cols[2], vertex.size = 8, vertex.label = NA)


### Edge weight distribution
set.vertex.attribute(net.1, "weight", get.vertex.attribute(net.1, "weight")/100)
par(mfrow=c(1, 2))
hist(get.edge.value(net.1, "edge_weight"), 
     main = "", col = cols.2[2], xlab = "Edge weights")


set.edge.value(net.1, "edge_weight", scale(sqrt(1-(get.edge.value(net.1, "edge_weight"))))[,1])
hist(get.edge.value(net.1, "edge_weight"), 
     main = "", col = cols.2[2], xlab = "Transformed Edge weights")

mtext("Edge Weight Distribution", line = -1.5, outer = T)
par(mfrow=c(1, 1))


### Model Estimation
model1.1 <- ergm(net.1~nodecov("age")+nodecov("weight")+nodefactor("parity")+nodefactor("farm")+nodefactor("treatment")+nodefactor("subtreatment"),
                 control = control.ergm(seed = 19201613, MCMC.samplesize = 5000), response ="edge_weight" , reference =~ StdNormal)
summary(model1.1)

model1.2 <- ergm(net.1~nodecov("age")+nodecov("weight")+nodefactor("parity")+nodefactor("farm")+nodefactor("treatment"),
                 control = control.ergm(seed = 19201613, MCMC.samplesize = 5000), response ="edge_weight" , reference =~ StdNormal)
summary(model1.2)

model1.3 <- ergm(net.1~nodecov("age")+nodefactor("parity")+nodefactor("farm")+nodefactor("treatment"),
                 control = control.ergm(seed = 19201613, MCMC.samplesize = 5000), response ="edge_weight" , reference =~ StdNormal)
summary(model1.3)


### cross-tab betweem farm and subtreatment
xtabs(~Net.M$Net_COV$farm+Net.M$Net_COV$treatment)
xtabs(~Net.M$Net_COV$farm+Net.M$Net_COV$subtreatment)

### aic and bic plot
model1_aic <- sapply(1:3, function(x) summary(eval(parse(text = paste0("model1.",x))))$aic)
model1_bic <- sapply(1:3, function(x) summary(eval(parse(text = paste0("model1.",x))))$bic)

par(mfrow = c(1,2))
plot(model1_aic, type = "b", col=cols.2[1], pch=19,
     main = "AIC", xlab = "Models", ylab = "AIC", xaxt="n")
axis(side = 1, at = seq(1, 3, by = 1))
plot(model1_bic, type = "b", col=cols.2[1], pch=19,
     main = "BIC", xlab = "Models", ylab = "BIC", xaxt="n")
axis(side = 1, at = seq(1, 3, by = 1))


### MCMC diagnostics
mcmc.diagnostics(model1.3, vars.per.page = 6)


### Odds ratio and confidence interval
ci <- exp(cbind(OR = coef(model1.3), confint(model1.3)))

##### ----- unweighted network 

### Model Estimation : Dyad Independent models
net.2 <- Net.M$unweighted_net
model2.1 <- ergm(net.2~edges,
                 control = control.ergm(seed = 19201613))
summary(model2.1)

model2.2 <- ergm(net.2~edges+nodecov("age")+nodecov("weight")+nodefactor("parity")+nodefactor("farm")+nodematch("farm")+nodefactor("treatment"),
                 control = control.ergm(seed = 19201613))
summary(model2.2)

model2.3 <- ergm(net.2~edges+nodecov("age")+nodecov("weight")+nodefactor("parity")+nodefactor("farm")+nodematch("farm"),
                 control = control.ergm(seed = 19201613))
summary(model2.3)

model2.4 <- ergm(net.2~edges+nodecov("age")+nodecov("weight")+nodefactor("farm")+nodematch("farm"),
                 control = control.ergm(seed = 19201613))
summary(model2.4)

model2.5 <- ergm(net.2~edges+nodecov("weight")+nodefactor("farm")+nodematch("farm"),
                 control = control.ergm(seed = 19201613))
summary(model2.5)


### aic, bic and R-squared plot
model2_aic <- sapply(1:5, function(x) summary(eval(parse(text = paste0("model2.",x))))$aic)
model2_bic <- sapply(1:5, function(x) summary(eval(parse(text = paste0("model2.",x))))$bic)

par(mfrow = c(1,3))
plot(model2_aic, type = "b", col=cols.2[1], pch=19,
     main = "AIC and BIC", xlab = "Models", ylab = "AIC and BIC", xaxt="n")
points(model2_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 5, by = 1))
legend("topright", legend = c("AIC", "BIC"), fill = c(cols.2[1], cols.2[2]))

model2_null.dev <- summary(model2.1)$devtable[1]
model2_res.dev <- sapply(1:5, function(x) summary(eval(parse(text = paste0("model2.",x))))$devtable[2])

plot(model2_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 5, by = 1))

plot(1-model2_res.dev/model2_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 5, by = 1))

model2_null.dev <- summary(model2.1)$devtable[1]
model2_res.dev <- sapply(1:5, function(x) summary(eval(parse(text = paste0("model2.",x))))$devtable[2])


par(mfrow = c(1,2))
plot(model2_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 5, by = 1))

plot(1-model2_res.dev/model2_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 5, by = 1))



### Model Estimation : Dyad dependent models
model3.1 <- ergm(net.2~edges+gwesp(0.1, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize =10000, MCMC.interval = 5000, MCMLE.check.degeneracy = T))
summary(model3.1)

model3.2 <- ergm(net.2~edges+nodecov("weight")+nodefactor("farm")+nodematch("farm")+gwesp(0.1, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize =10000, MCMC.interval = 5000, MCMLE.check.degeneracy = T))
summary(model3.2)

model3.3 <- ergm(net.2~edges+nodecov("weight")+nodematch("farm")+gwesp(0.1, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize =10000, MCMC.interval = 5000, MCMLE.check.degeneracy = T))
summary(model3.3)

model3.4 <- ergm(net.2~edges+nodecov("weight")+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize =10000, MCMC.interval = 5000, MCMLE.check.degeneracy = T))
summary(model3.4)

model3.5 <- ergm(net.2~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize =10000, MCMC.interval = 5000, MCMLE.check.degeneracy = T))


### aic, bic and R-squared plot
model3_aic <- sapply(c(1,2,4), function(x) summary(eval(parse(text = paste0("model3.",x))))$aic)
model3_bic <- sapply(c(1,2,4), function(x) summary(eval(parse(text = paste0("model3.",x))))$bic)

par(mfrow = c(1,3))
plot(model3_aic, type = "b", col=cols.2[1], pch=19,
     main = "AIC and BIC", xlab = "Models", ylab = "AIC", xaxt="n")
points(model3_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 5, by = 1))
model3_null.dev <- summary(model3.1)$devtable[1]
model3_res.dev <- sapply(c(1,2,4), function(x) summary(eval(parse(text = paste0("model3.",x))))$devtable[2])

plot(model3_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 5, by = 1))

plot(1-model3_res.dev/model3_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 5, by = 1))

### mcmc.diagnostics plot
mcmc.diagnostics(model3.4)


### GOF
gof.3 <- gof(model3.1, control = control.gof.ergm(seed = 19201613))
gof.4 <- gof(model3.2, control = control.gof.ergm(seed = 19201613))
gof.5 <- gof(model3.3, control = control.gof.ergm(seed = 19201613))
gof.6 <- gof(model3.4, control = control.gof.ergm(seed = 19201613))

### GOF Plots
par(mfrow = c(3,4))
plot(gof.3)
plot(gof.4)
plot(gof.6)


### Odds ratio and confidence interval
ci_2 <- exp(cbind(OR = coef(model3.4), confint(model3.4)))

############################################################################





########################## Evening Network ##################################

##### ----- complete network

net.3 <- Net.E$weighted_net

### Edge weight distribution
set.vertex.attribute(net.3, "weight", get.vertex.attribute(net.3, "weight")/100)

par(mfrow=c(1, 2))
hist(get.edge.value(net.3, "edge_weight"), 
     main = "", col = cols.2[2], xlab = "Edge weights")

set.edge.value(net.3, "edge_weight", scale(sqrt(1-(get.edge.value(net.3, "edge_weight"))))[,1])

hist(get.edge.value(net.3, "edge_weight"), 
     main = "", col = cols.2[2], xlab = "Transformed Edge weights")

mtext("Edge Weight Distribution", line = -1.5, outer = T)
par(mfrow=c(1, 1))


### Model Estimation
model4.1 <- ergm(net.3~nodecov("age")+nodecov("weight")+nodefactor("parity")+nodefactor("farm")+nodefactor("treatment"),
                 control = control.ergm(seed = 19201613, MCMC.samplesize = 5000), response ="edge_weight" , reference =~ StdNormal)
summary(model4.1)

model4.2 <- ergm(net.3~nodecov("age")+nodefactor("parity")+nodefactor("farm")+nodefactor("treatment"),
                 control = control.ergm(seed = 19201613, MCMC.samplesize = 5000), response ="edge_weight" , reference =~ StdNormal)
summary(model4.2)



##### ----- unweighted network


net.4 <- Net.E$unweighted_net

### Model Estimation
model5.1 <- ergm(net.4~edges,
                 control = control.ergm(seed = 19201613))
summary(model5.1)

model5.2 <- ergm(net.4~edges+nodecov("age")+nodefactor("parity")+nodefactor("farm")+nodematch("farm")+nodefactor("treatment"),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50))
summary(model5.2)

model5.3 <- ergm(net.4~edges+nodecov("age")+nodefactor("parity")+nodefactor("farm")+nodematch("farm")+nodefactor("treatment")+gwdsp(.1, T)+gwesp(.1, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50))
summary(model5.3)

model5.4 <- ergm(net.4~edges+nodecov("age")+nodefactor("parity")+nodefactor("farm")+nodematch("farm")+gwdsp(.1, T)+gwesp(.1, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize = 2000, MCMLE.check.degeneracy = T))
summary(model5.4)

model5.5 <- ergm(net.4 ~ edges + nodecov("age") + nodefactor("parity") + nodefactor("farm") + nodematch("farm") + gwdsp(0.25, T) + gwesp(0.25, T),
                 control = control.ergm(seed = 19201613,MCMLE.maxit = 50, MCMC.samplesize = 2000, MCMLE.check.degeneracy = T))
summary(model5.5)

model5.6 <- ergm(net.4 ~ edges + nodecov("age") + nodefactor("parity") + nodefactor("farm") + nodematch("farm") + gwdsp(0.5, T) + gwesp(0.5, T),
                 control = control.ergm(seed = 19201613, MCMLE.maxit = 50, MCMC.samplesize = 2000, MCMLE.check.degeneracy = T))
summary(model5.6)


### aic, bic and R-squared plot
model5_aic <- sapply(1:6, function(x) summary(eval(parse(text = paste0("model5.",x))))$aic)
model5_bic <- sapply(1:6, function(x) summary(eval(parse(text = paste0("model5.",x))))$bic)

par(mfrow = c(1,3))
plot(model5_aic, type = "b", col=cols.2[1], pch=19,
     main = "AIC and BIC", xlab = "Models", ylab = "AIC", xaxt="n")
points(model5_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 6, by = 1))
model5_null.dev <- summary(model5.2)$devtable[1]
model5_res.dev <- sapply(1:6, function(x) summary(eval(parse(text = paste0("model5.",x))))$devtable[2])

plot(model5_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 6, by = 1))

plot(1-model5_res.dev/model5_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 6, by = 1))

### MCMC.diagnostics
mcmc.diagnostics(model5.4)


### GOF
gof.7 <- gof(model5.1, control = control.gof.ergm(seed = 19201613))
gof.8 <- gof(model5.4, control = control.gof.ergm(seed = 19201613))
gof.9 <- gof(model5.6, control = control.gof.ergm(seed = 19201613))
gof.10 <- gof(model5.5, control = control.gof.ergm(seed = 19201613))


### GOF Plot
par(mfrow = c(2,4))
plot(gof.7)
plot(gof.8)

############################################################################
