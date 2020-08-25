########## Libraries ##########

library(igraph)
library(statnet)
library(ergm.count)
library(intergraph)
library(RColorBrewer)
#library(kableExtra)
#library(stargazer)
library(ndtv)
library(htmlwidgets)

###############################




################### Load Data sets ####################

load("Net_M_Feb.RData")
load("Net_M_FebnMar.RData")
load("Net_M_FebnApr.RData")
load("Net_E_Feb.RData")
load("Net_E_FebnMar.RData")
load("Net_E_FebnApr.RData")

cols <- brewer.pal(8, "Set2")
cols.2 <- brewer.pal(9, "Set1")

######################################################




########################## Morning Network ##################################

### Make Dynamic Network 

tmp <- Net_M_Feb$unweighted_net
set.network.attribute(tmp, "vertex.pid", "vertex.names")

tmp2 <- Net_M_FebnMar$unweighted_net
set.network.attribute(tmp2, "vertex.pid", "vertex.names")

tmp3 <- Net_M_FebnApr$unweighted_net
set.network.attribute(tmp3, "vertex.pid", "vertex.names")

temp.net <- networkDynamic(network.list = list(tmp, tmp2, tmp3), vertex.pid = get.network.attribute(tmp2, "vertex.pid"), create.TEAs = TRUE, vertex.TEA.names = "farm")


### Plot of each network
g1 <- asIgraph(Net_M_Feb$unweighted_net)
g2 <- asIgraph(Net_M_FebnMar$unweighted_net)
g3 <- asIgraph(Net_M_FebnApr$unweighted_net)

par(mfrow=c(1,3))

plot(g1, vertex.color=V(g1)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms : Feb")
legend("topleft", title = "Farms", legend = levels(as.factor(Net_M_Feb$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

plot(g2, vertex.color=V(g2)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms : Mar")
legend("topleft", title = "Farms", legend = levels(as.factor(Net_M_FebnMar$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

plot(g3, vertex.color=V(g3)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms : Apr")
legend("topleft", title = "Farms", legend = levels(as.factor(Net_M_FebnApr$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])


### Degree distribution
par(mfrow=c(1,3))
hist(igraph::degree(g1), col =cols.2[2], main="February : Degree Distribution", 
     xlab = "Degree")

hist(igraph::degree(g2), col =cols.2[2], main="March : Degree Distribution", 
     xlab = "Degree")

hist(igraph::degree(g3), col =cols.2[2], main="April : Degree Distribution", 
     xlab = "Degree")


### Density and Transitivity
data.frame("Network" = c("Feb", "Mar", "Apr"), 
           "Density" = sapply(1:3, function(x) edge_density(eval(parse(text = paste0("g",x)))))
           , "Transitivity" = sapply(1:3, function(x) transitivity(eval(parse(text = paste0("g",x))))))


#### Assortativity of largest components
lgc <- decompose.graph(g1)[[1]]
lgc2 <- decompose.graph(g2)[[4]]
lgc3 <- decompose.graph(g3)[[1]]

assort.df<- cbind("Farm" = assortativity.nominal(lgc, V(lgc)$farm),
                  "Parity" = assortativity.nominal(lgc, V(lgc)$parity))

assort.df2<- cbind("Farm" = assortativity.nominal(lgc2, V(lgc2)$farm),
                   "Parity" = assortativity.nominal(lgc2, V(lgc2)$parity))

assort.df3<- cbind("Farm" = assortativity.nominal(lgc3, V(lgc3)$farm),
                   "Parity" = assortativity.nominal(lgc3, V(lgc3)$parity))

### Visualize change 
render.d3movie(temp.net, output.mode = 'htmlWidget')


### number of edges and triangle in each network
tErgmStats(temp.net, "~edges+triangle")




## Models
temp.fit1 <- stergm(temp.net,
                    formation = ~edges,
                    dissolution = ~edges,
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))

temp.fit2 <- stergm(temp.net,
                    formation = ~edges+nodefactor("farm")+nodematch("farm"),
                    dissolution = ~edges+nodefactor("farm")+nodematch("farm"),
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))

temp.fit3 <- stergm(temp.net,
                    formation = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T),
                    dissolution = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T),
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))

temp.fit4 <- stergm(temp.net,
                    formation = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T)+gwdsp(0.1, T),
                    dissolution = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T)+gwdsp(0.1,T),
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))


### AIC, BIC and R-SQD
form1_aic <- sapply(1:4, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$formation$aic)
form1_bic <- sapply(1:4, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$formation$bic)

diss1_aic <- sapply(1:4, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$dissolution$aic)
diss1_bic <- sapply(1:4, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$dissolution$bic)

par(mfrow = c(2,3))
plot(form1_aic, type = "b", col=cols.2[1], pch=19,
     main = "Formation Models: AIC and BIC", xlab = "Models", ylab = "AIC and BIC", xaxt="n")
points(form1_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 4, by = 1))

form1_null.dev <- summary(temp.fit)$formation$devtable[1]
form1_res.dev <- sapply(1:4, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$formation$devtable[2])

plot(form1_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Formation Models: Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 4, by = 1))

plot(1-form1_res.dev/form1_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "Formation Models:  R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 4, by = 1))

plot(diss1_aic , type = "b", col=cols.2[1], pch=19,
     main = "Dissolution Models:  AIC and BIC", xlab = "Models", ylab = "AIC and BIC", xaxt="n")
points(diss1_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 4, by = 1))

diss1_null.dev <- summary(temp.fit)$dissolution$devtable[1]
diss1_res.dev <- sapply(1:4, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$dissolution$devtable[2])


plot(diss1_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Dissolution Models: Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 4, by = 1))

plot(1-diss1_res.dev/diss1_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "Dissolution Models:  R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 4, by = 1))


### MCMC Diagnostics
mcmc.diagnostics(temp.fit4)

### GOF 
gof1 <- gof(temp.fit)
gof2 <- gof(temp.fit2)
gof3 <- gof(temp.fit3)
gof4 <- gof(temp.fit4)


### GOF Plots
par(mfrow=c(2,4))
plot(gof1, main = NULL)

par(mfrow=c(2,4))
plot(gof4, main = NULL)

########################################################################################




########################## Evening  Network ##################################

## Make Dynamic Network 
tmp4 <- Net_E_Feb$unweighted_net
set.network.attribute(tmp4, "vertex.pid", "vertex.names")

tmp5 <- Net_E_FebnMar$unweighted_net
set.network.attribute(tmp5, "vertex.pid", "vertex.names")

tmp6 <- Net_E_FebnApr$unweighted_net
set.network.attribute(tmp6, "vertex.pid", "vertex.names")

temp.net2 <- networkDynamic(network.list = list(tmp4, tmp5, tmp6), vertex.pid = get.network.attribute(tmp5, "vertex.pid"), create.TEAs = TRUE, vertex.TEA.names = "farm")


### Plot of each network
g1 <- asIgraph(Net_E_Feb$unweighted_net)
g2 <- asIgraph(Net_E_FebnMar$unweighted_net)
g3 <- asIgraph(Net_E_FebnApr$unweighted_net)

par(mfrow=c(1,3))

plot(g1, vertex.color=V(g1)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms")
legend("topleft", title = "Farms", legend = levels(as.factor(Net_M_Feb$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

plot(g2, vertex.color=V(g2)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms")
legend("topleft", title = "Farms", legend = levels(as.factor(Net_M_FebnMar$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

plot(g3, vertex.color=V(g3)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms")
legend("topleft", title = "Farms", legend = levels(as.factor(Net_M_FebnApr$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])


### Degree distribution
par(mfrow=c(1,3))
hist(igraph::degree(g1), col =cols.2[2], main="February : Degree Distribution", 
     xlab = "Degree")

hist(igraph::degree(g2), col =cols.2[2], main="March : Degree Distribution", 
     xlab = "Degree")

hist(igraph::degree(g3), col =cols.2[2], main="April : Degree Distribution", 
     xlab = "Degree")


### Density and Transitivity
data.frame("Network" = c("Feb", "Mar", "Apr"), 
           "Density" = sapply(1:3, function(x) edge_density(eval(parse(text = paste0("g",x)))))
           , "Transitivity" = sapply(1:3, function(x) transitivity(eval(parse(text = paste0("g",x))))))



### Visualize network change over time
render.d3movie(temp.net2, output.mode = 'htmlWidget')


### number of edges and triangle in each network
tErgmStats(temp.net2, "~edges+triangle")


## Models
temp.fit5 <- stergm(temp.net2,
                    formation = ~edges,
                    dissolution = ~edges,
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))

temp.fit6 <- stergm(temp.net2,
                    formation = ~edges+nodefactor("farm")+nodematch("farm"),
                    dissolution = ~edges+nodefactor("farm")+nodematch("farm"),
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))

temp.fit7 <- stergm(temp.net2,
                    formation = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T),
                    dissolution = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T),
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))

temp.fit8 <- stergm(temp.net2,
                    formation = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T)+gwdsp(0.1, T),
                    dissolution = ~edges+nodefactor("farm")+nodematch("farm")+gwesp(0.15, T)+gwdsp(0.1,T),
                    estimate = "CMLE",
                    times = 0:2,
                    control = control.stergm(seed = 19201613))


### AIC, BIC and R-SQD
form2_aic <- sapply(5:8, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$formation$aic)
form2_bic <- sapply(5:8, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$formation$bic)

diss2_aic <- sapply(5:8, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$dissolution$aic)
diss2_bic <- sapply(5:8, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$dissolution$bic)

par(mfrow = c(2,3))
plot(form2_aic, type = "b", col=cols.2[1], pch=19,
     main = "Formation Models: AIC and BIC", xlab = "Models", ylab = "AIC and BIC", xaxt="n")
points(form2_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 4, by = 1))

form2_null.dev <- summary(temp.fit5)$formation$devtable[1]
form2_res.dev <- sapply(5:8, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$formation$devtable[2])

plot(form2_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Formation Models: Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 4, by = 1))

plot(1-form2_res.dev/form2_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "Formation Models:  R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 4, by = 1))

plot(diss2_aic , type = "b", col=cols.2[1], pch=19,
     main = "Dissolution Models:  AIC and BIC", xlab = "Models", ylab = "AIC and BIC", xaxt="n")
points(diss2_bic, type = "b", col=cols.2[2], pch=19)
axis(side = 1, at = seq(1, 4, by = 1))

diss2_null.dev <- summary(temp.fit5)$dissolution$devtable[1]
diss2_res.dev <- sapply(5:8, function(x) summary(eval(parse(text = paste0("temp.fit",x))))$dissolution$devtable[2])

plot(diss2_res.dev, type = "b", col=cols.2[1], pch=19,
     main = "Dissolution Models: Residual Deviance", xlab = "Models", ylab = "Residual Deviance", xaxt="n")

axis(side = 1, at = seq(1, 4, by = 1))

plot(1-diss2_res.dev/diss2_null.dev, type = "b", col=cols.2[1], pch=19,
     main = "Dissolution Models:  R-Squared", xlab = "Models", ylab = "R-Squared", xaxt="n")
axis(side = 1, at = seq(1, 4, by = 1))


### GOF
gof5 <- gof(temp.fit5)
gof6 <- gof(temp.fit6)
gof7 <- gof(temp.fit7)
gof8 <- gof(temp.fit8)


### GOF - plots
par(mfrow=c(2,4))
plot(gof5, main = NULL)

par(mfrow=c(2,4))
plot(gof7, main = NULL)

########################################################################################