########## Libraries ##########

library(igraph)
library(blockmodels)
library(igraph)
library(intergraph)
library(RColorBrewer)
library(kableExtra)
library(GGally)
library(ggnet)
library(stargazer)
library(corrplot)
library(sna)
library(e1071)

###############################




################### Load Data sets ####################

load("Net_M_Feb.RData")
load("Net_E_Feb.RData")
cols <- brewer.pal(8, "Set2")
cols.2 <- brewer.pal(9, "Set1")
Net.M <- Net_M_Feb
Net.E <- Net_E_Feb

######################################################




######################## Milk yield Network ########################

### Network plot
net.2 <- Net.M$unweighted_net
net.4 <- Net.E$unweighted_net


par(mfrow=c(1,2))
plot(asIgraph(net.2), vertex.color = cols.2[2], vertex.size = 8, vertex.label = NA, main="Morning Network")
plot(asIgraph(net.4), vertex.color = cols.2[2], vertex.size = 8, vertex.label = NA, main="Evening Network")


### Degree distribution
g1 <- asIgraph(net.2)
g2 <- asIgraph(net.4)

edge_density(g1)
edge_density(g2)


par(mfrow=c(1,2))
hist(igraph::degree(g1), probability = T, col =cols.2[2], main="Morning : Degree Distribution", 
     xlab = "Degree")
lines(density(igraph::degree(g1)))
hist(igraph::degree(g2), probability = T, col =cols.2[2], main="Evening : Degree Distribution", 
     xlab = "Degree")
lines(density(igraph::degree(g2)), lwd=2)

table(igraph::degree(g1))
table(igraph::degree(g2))


### Transitivity 
transitivity(g1)
transitivity(g2)


### Add age and BCS cat to network covariate
Net.M$Net_COV <- transform(Net.M$Net_COV, Age_cat = cut(Net.M$Net_COV$Age, breaks = c(0,1.99,2.99,3.99,4.99, max(Net.M$Net_COV$Age,na.rm = T)), labels = 1:5))
Net.M$Net_COV <- transform(Net.M$Net_COV, BCS_cat = cut(Net.M$Net_COV$BCS, breaks = c(2.49,2.74,2.99,3.24, 3.50), labels = 1:4))


set.vertex.attribute(net.2, "age_cat", as.character(Net.M$Net_COV$Age_cat))
set.vertex.attribute(net.2, "BCS_cat", as.character(Net.M$Net_COV$BCS_cat))

Net.E$Net_COV <- transform(Net.E$Net_COV, Age_cat = cut(Net.E$Net_COV$Age, breaks = c(0,1.99,2.99,3.99,4.99, max(Net.E$Net_COV$Age,na.rm = T)), labels = 1:5))
Net.E$Net_COV <- transform(Net.E$Net_COV, BCS_cat = cut(Net.E$Net_COV$BCS, breaks = c(2.49,2.74,2.99,3.24, 3.50), labels = 1:4))


set.vertex.attribute(net.4, "age_cat", as.character(Net.E$Net_COV$Age_cat))
set.vertex.attribute(net.4, "BCS_cat", as.character(Net.E$Net_COV$BCS_cat))


### largest component plots by categories (Morning)
set.seed(19201613)
g1 <- asIgraph(net.2)
style <- layout.fruchterman.reingold(g1)

par(mfrow = c(3,2))
plot(g1, vertex.color=V(g1)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms", layout=style)
legend("topleft", title = "Farms", legend = levels(as.factor(Net.M$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

plot(g1, vertex.color=brewer.pal(9, "Set1")[as.integer(V(g1)$parity)], vertex.label = NA, vertex.size = 5, asp =0, main = "Parity", layout=style)
legend("topleft", title = "Parity", legend = levels(as.factor(Net.M$Net_COV$parity)), fill =brewer.pal(9, "Set1"))

plot(g1, vertex.color=V(g1)$BCS_cat, vertex.label = NA, vertex.size = 5, asp =0, main = "BCS Categories", layout=style)
legend("topleft", title = "BCS categories", legend = levels(as.factor(Net.M$Net_COV$BCS_cat)), fill = categorical_pal(8)[1:5])

plot(g1, vertex.color=V(g1)$age_cat, vertex.label = NA, vertex.size = 5, asp =0, main = "Age Categories", layout=style)
legend("topleft", title = "Age categories", legend = levels(as.factor(Net.M$Net_COV$Age_cat)), fill = categorical_pal(8)[1:6])

plot(g1, vertex.color=as.integer(V(g1)$treatment), vertex.label = NA, vertex.size = 5, asp =0, main = "Treatments", layout=style)
legend("topleft", title = "Treatment", legend = 1:5, fill = categorical_pal(8)[1:5])

plot(g1, vertex.color=V(g1)$subtreatment, vertex.label = NA, vertex.size = 5, asp =0, main = "Subtreatments", layout=style)
legend("topleft", title = "Subtreatment", legend = levels(as.factor(Net.M$Net_COV$subtreatment)), fill = c("white", categorical_pal(8)[c(1:4)]))



### largest component plot by categories (Evening)
set.seed(19201613)
g2 <- asIgraph(net.4)
style <- layout.fruchterman.reingold(g2)


par(mfrow = c(3,2))
plot(g2, vertex.color=V(g2)$farm, vertex.label = NA, vertex.size = 5, asp =0, main = "Farms", layout=style)
legend("topleft", title = "Farms", legend = levels(as.factor(Net.E$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

plot(g2, vertex.color=brewer.pal(9, "Set1")[as.integer(V(g2)$parity)], vertex.label = NA, vertex.size = 5, asp =0, main = "Parity", layout=style)
legend("topleft", title = "Parity", legend = levels(as.factor(Net.E$Net_COV$parity)), fill =brewer.pal(9, "Set1"))

plot(g2, vertex.color=V(g2)$BCS_cat, vertex.label = NA, vertex.size = 5, asp =0, main = "BCS Categories", layout=style)
legend("topleft", title = "BCS categories", legend = levels(as.factor(Net.E$Net_COV$BCS_cat)), fill = categorical_pal(8)[1:5])

plot(g2, vertex.color=V(g2)$age_cat, vertex.label = NA, vertex.size = 5, asp =0, main = "Age Categories", layout=style)
legend("topleft", title = "Age categories", legend = levels(as.factor(Net.E$Net_COV$Age_cat)), fill = categorical_pal(8)[1:6])

plot(g2, vertex.color=as.integer(V(g2)$treatment), vertex.label = NA, vertex.size = 5, asp =0, main = "Treatments", layout=style)
legend("topleft", title = "Treatment", legend = 1:5, fill = categorical_pal(8)[1:5])

plot(g2, vertex.color=V(g2)$subtreatment, vertex.label = NA, vertex.size = 5, asp =0, main = "Subtreatments", layout=style)
legend("topleft", title = "Subtreatment", legend = levels(as.factor(Net.E$Net_COV$subtreatment)), fill = c("white", categorical_pal(8)[c(1:4)]))


### Assortativity Mixing
lgc <- decompose.graph(g1)[[1]]
lgc2 <- decompose.graph(g2)[[1]]
assort.df<- cbind("Farm" = assortativity.nominal(lgc, V(lgc)$farm),
                  "Parity" = assortativity.nominal(lgc, V(lgc)$parity),
                  "Age Categories" = assortativity.nominal(lgc, V(lgc)$age_cat),
                  "BCS Categories" = assortativity.nominal(lgc, V(lgc)$BCS_cat))

assort.df2<- cbind("Farm" = assortativity.nominal(lgc2, V(lgc2)$farm),
                   "Parity" = assortativity.nominal(lgc2, V(lgc2)$parity),
                   "Age Categories" = assortativity.nominal(lgc2, V(lgc2)$age_cat),
                   "BCS Categories" = assortativity.nominal(lgc2, V(lgc2)$BCS_cat))

####################################################################



######################## Community Detection Algorithms ########################

##### ----- Girvan-Newman -----
set.seed(19201613)
lgc <- decompose.graph(g1)[[1]]
cluster <- cluster_edge_betweenness(lgc)

lgc2 <- decompose.graph(g2)[[1]]
cluster2 <- cluster_edge_betweenness(lgc2)

par(mfrow=c(1,2))
plot(as.dendrogram(cluster), main="Morning Network")
plot(as.dendrogram(cluster2), main="Evening Network")

memberships <- cut_at(cluster, 4)
memberships2 <- cut_at(cluster2, 4)

# Morning
style <- layout.fruchterman.reingold(lgc)

par(mfrow=c(1,2))
plot(lgc, vertex.color=V(lgc)$farm, vertex.label = NA, vertex.size = 5, main = "original", layout=style)
legend("topleft", title = "Farms", legend = levels(as.factor(Net.M$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])
plot(lgc, vertex.label = NA, vertex.size = 5, vertex.color = memberships,  main="Detected", layout=style)

table(V(lgc)$farm, memberships)
classAgreement(table(V(lgc)$farm, memberships))$crand


# Evening
style <- layout.fruchterman.reingold(lgc2)

par(mfrow=c(1,2))
plot(lgc2, vertex.color=V(lgc2)$farm, vertex.label = NA, vertex.size = 5, main = "original", layout=style)
legend("topleft", title = "Farms", legend = levels(as.factor(Net.E$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])
plot(lgc2, vertex.label = NA, vertex.size = 5, vertex.color = memberships2,  main="Detected", layout=style)

table(V(lgc2)$farm, memberships2)
classAgreement(table(V(lgc2)$farm, memberships2))$crand


##### ----- Heirarchical Clustering -----
set.seed(19201613)
lgc <- decompose.graph(g1)[[1]]
Y <- as.matrix(get.adjacency(lgc)) 
dissimilarity <- dist(Y)
# cluster <- hclust(dissimilarity, method = "single") 
cluster <- hclust(dissimilarity, method = "complete") 
# cluster <- hclust(dissimilarity, method = "average") 

lgc2 <- decompose.graph(g2)[[1]]
Y2 <- as.matrix(get.adjacency(lgc2)) 
dissimilarity2 <- dist(Y2)
cluster2 <- hclust(dissimilarity2, method = "complete") 


par(mfrow=c(1,2))
plot(cluster, labels = FALSE)
abline(h=5.3, col=cols[2])
plot(cluster2, labels = FALSE)
abline(h=6.5, col=cols[2])


memberships <- cutree(cluster, k = 4)
memberships2 <- cutree(cluster2, k = 3)

style <- layout.fruchterman.reingold(lgc)

par(mfrow=c(1,2))
plot(lgc, vertex.color=V(lgc)$farm, vertex.label = NA, vertex.size = 5, main = "Original", layout=style)
plot(lgc, vertex.label = NA, vertex.size = 5, vertex.color = memberships, layout=style, main="Detected")

table(V(lgc)$farm, memberships)
classAgreement(table(V(lgc)$farm, memberships))$crand

style <- layout.fruchterman.reingold(lgc2)

par(mfrow=c(1,2))
plot(lgc2, vertex.color=V(lgc2)$farm, vertex.label = NA, vertex.size = 5, main = "Original", layout=style)
plot(lgc2, vertex.label = NA, vertex.size = 5, vertex.color = memberships2, layout=style, main="Detected")

table(V(lgc2)$farm, memberships2)
classAgreement(table(V(lgc2)$farm, memberships2))$crand


##### ----- Louvain Algorithm -----
set.seed(19201613)
cluster <- cluster_louvain(lgc)
memberships <- membership(cluster)
style <- layout.fruchterman.reingold(lgc)

par(mfrow=c(1,2))
plot(lgc, vertex.color=V(lgc)$farm, vertex.label = NA, vertex.size = 5, main = "Original", layout=style)
plot(lgc, vertex.label = NA, vertex.size = 5, vertex.color = memberships, layout=style, main="Detected")

table(V(lgc)$farm, memberships)
classAgreement(table(V(lgc)$farm, memberships))$crand

set.seed(19201613)
cluster2 <- cluster_louvain(lgc2)
memberships2 <- membership(cluster2)
style <- layout.fruchterman.reingold(lgc2)

par(mfrow=c(1,2))
plot(lgc2, vertex.color=V(lgc2)$farm, vertex.label = NA, vertex.size = 5, main = "Original", layout=style)
plot(lgc2, vertex.label = NA, vertex.size = 5, vertex.color = memberships2, layout=style, main="Detected")

table(V(lgc2)$farm, memberships2)
classAgreement(table(V(lgc2)$farm, memberships2))$crand

################################################################################





#################################### Stochastic Block Model ####################################
g1 <- asIgraph(net.2)
g2 <- asIgraph(net.4)

par(mfrow=c(1,2))
image(Net.M$unweighted_matrix)
image(Net.E$unweighted_matrix)

set.seed(19201613)
farms <- V(g1)$farm
memberships <- make_clusters(g1, as.numeric(farms)) # create a communities object using farm
style <- layout.fruchterman.reingold(g1)


farms2 <- V(g2)$farm
memberships2 <- make_clusters(g2, as.numeric(farms)) # create a communities object using farm
style <- layout.fruchterman.reingold(g2)


set.seed(19201613)
sbm <- BM_bernoulli(membership_type = "SBM_sym",
                    adj = Net.M$unweighted_matrix, 
                    verbosity = 1, 
                    plotting="",
                    explore_min = 1,
                    explore_max = 9) 

sbm$estimate()
soft_clustering <- sbm$memberships[[6]]$Z 
hard_clustering <- apply(soft_clustering, 1, which.max) 


sbm2 <- BM_bernoulli(membership_type = "SBM_sym",
                     adj = Net.E$unweighted_matrix, 
                     verbosity = 1, 
                     plotting="",
                     explore_min = 1,
                     explore_max = 9) 

sbm2$estimate()
soft_clustering2 <- sbm2$memberships[[6]]$Z 
hard_clustering2 <- apply(soft_clustering2, 1, which.max) 

## ICL Plot
par(mfrow=c(1,2))
plot(sbm$ICL, pch=19, type="b", col=cols.2[1])
plot(sbm2$ICL, pch=19, type="b", col=cols.2[1])

## Block composition plot
par(mfrow=c(1,2))
sbm$memberships[[6]]$plot()
sbm2$memberships[[6]]$plot()

## Connection Probability 
par(mfrow=c(1,2))
sbm$plot_parameters(6)
sbm2$plot_parameters(6) 

set.seed(19201613)

memberships <- make_clusters(g1, hard_clustering) # create the new communities object according to partition found
memberships2 <- make_clusters(g2, hard_clustering2) # create the new communities object according to partition found

## SBM vs Farm (Morning)
par(mfrow=c(1,2))
style <- layout.fruchterman.reingold(g1)
plot(memberships, mark.groups = NULL, edge.color = NULL, g1, vertex.label = NA, vertex.size = 5, layout = style, main = "SBM")
legend("bottomleft", title = "Cluster", legend = levels(as.factor(memberships$membership)), fill = categorical_pal(8))

plot(g1, vertex.color=V(g1)$farm, vertex.label = NA, vertex.size = 5, main = "Farms", layout=style)
legend("bottomleft", title = "Farms", legend = levels(as.factor(Net.M$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

## SBM vs Farm (Evening)
par(mfrow=c(1,2))
style <- layout.fruchterman.reingold(g2)
plot(memberships2, mark.groups = NULL, edge.color = NULL, g2, vertex.label = NA, vertex.size = 5, layout = style, main = "SBM")
legend("bottomleft", title = "Cluster", legend = levels(as.factor(memberships2$membership)), fill = categorical_pal(8))

plot(g2, vertex.color=V(g2)$farm, vertex.label = NA, vertex.size = 5, main = "Farms", layout=style)
legend("bottomleft", title = "Farms", legend = levels(as.factor(Net.M$Net_COV$farm)), fill = categorical_pal(8)[c(1,2,3,6)])

table(V(g1)$farm, memberships$membership)
classAgreement(table(V(g1)$farm, memberships$membership))$crand

table(V(g2)$farm, memberships2$membership)
classAgreement(table(V(g2)$farm, memberships2$membership))$crand

################################################################################################