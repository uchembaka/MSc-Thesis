########## Libraries ##########

library(birk)
library(DescTools)
library(lubridate)
library(RColorBrewer)
library(slickR)
library(htmlwidgets)
library(corrplot)

###############################



############# Load Data sets ##############

load("Net_M_Feb.RData")
cols <- brewer.pal(8, "Set2")
cols.2 <- brewer.pal(9, "Set1")
Net <- Net_M_Feb

##########################################




#################### Plots and Tables ####################

### unweighted network adjacency matrix, correlation plots
corrplot(Net$unweighted_matrix, cl.pos = "n", tl.pos = "n")
tmp.lo <- Net$Loess_data
tmp2 <- c(1,9,6,32)
tmp <- c(3,7,39,15)

par(mfrow=c(2, 2))
for (i in 1:4) {
  plot(seq(0, 30, length=20), tmp.lo[,tmp2[i]], type = "b", col=cols.2[1], pch = 19, 
       ylim = c(0, 25), ylab = "Yield", xlab = "Time",
       main = paste0("Correlation : ",round(cor(tmp.lo[,tmp[i]], tmp.lo[,tmp2[i]]),2)))
  points(seq(0, 30, length=20), tmp.lo[,tmp[i]], col = cols.2[2], pch = 19)
  lines(seq(0, 30, length=20), tmp.lo[,tmp[i]], col = cols.2[2])
  legend("bottomright", legend = c( paste0("Vertex ",tmp[i]),  paste0("Vertex ",tmp2[i])), fill =c(cols.2[2], cols.2[1]))
}
par(mfrow=c(1, 1))


### SCC
scc_lab <- c("scc<=200", "scc_201-300", "scc_301-400", "scc_>400")

SCC_RF_cat <- cut(Net$Net_COV$SCC_RF, breaks = c(0, 200, 300, 400, max(Net$Net_COV$SCC_RF, na.rm = T)), labels = scc_lab)
SCC_LF_cat <- cut(Net$Net_COV$SCC_LF, breaks = c(0, 200, 300, 400, max(Net$Net_COV$SCC_LF, na.rm = T)), labels = scc_lab)
SCC_RH_cat <- cut(Net$Net_COV$SCC_RH, breaks = c(0, 200, 300, 400, max(Net$Net_COV$SCC_RH, na.rm = T)), labels = scc_lab)
SCC_LH_cat <- cut(Net$Net_COV$SCC_LH, breaks = c(0, 200, 300, 400, max(Net$Net_COV$SCC_LH, na.rm = T)), labels = scc_lab)

Net_Cov <- transform(Net$Net_COV, SCC_RF_cat=SCC_RF_cat, SCC_LF_cat=SCC_LF_cat, SCC_RH_cat = SCC_RH_cat, SCC_LH_cat=SCC_LH_cat)

teat <- c("Right-Front", "Left-Front", "Right-Hind", "Left-Hind")

par(mfrow=c(2, 2))
for (i in 1:4) {
  p1 <- barplot(table(Net_Cov[,25+i]),
                xlab = "SCC Group",
                ylab = "Frequency",
                main = paste0(teat[i]),
                ylim = c(0, 150),
                col = cols)
  text(p1,table(Net_Cov[,25+i]), labels =table(Net_Cov[,25+i]), pos = 3)
}

par(mfrow=c(1, 1))


### Farm
p2 <- barplot(table(Net$Net_COV$farm),
              xlab = "Farm",
              ylab = "Frequency",
              main = "",
              ylim = c(0, 45),
              col = cols)
text(p2,table(Net$Net_COV$farm), labels =table(Net$Net_COV$farm), pos = 3)
text(p2,table(Net$Net_COV$farm), labels =paste0(round(prop.table(table(Net$Net_COV$farm))*100), "%"), pos = 1)


### Age, BCS and weight Summary 
kable(t(round(do.call(cbind, lapply(Net$Net_COV[c("Age", "BCS", "weight")], summary)),2)), "latex")


### Age plot
Net_Cov <- transform(Net_Cov, Age_cat = cut(Net_Cov$Age, breaks = c(0,1.99,2.99,3.99,4.99, max(Net_Cov$Age,na.rm = T)), labels = 1:5))
p3 <- barplot(table(Net_Cov$Age_cat),
              xlab = "Age group",
              ylab = "Frequency",
              main = "",
              ylim = c(0, 50),
              col = cols)
text(p3,table(Net_Cov$Age_cat), labels =table(Net_Cov$Age_cat), pos = 3)
text(p3,table(Net_Cov$Age_cat), labels =paste0(round(prop.table(table(Net_Cov$Age_cat))*100), "%"), pos = 1)


### BCS Plot
Net_Cov <- transform(Net_Cov, BCS_cat = cut(Net_Cov$BCS, breaks = c(2.5,2.75,3.0,3.25, 3.50), labels = 1:4))
p4 <- barplot(table(Net_Cov$BCS_cat),
              xlab = "BCS",
              ylab = "Frequency",
              main = "",
              ylim = c(0, 70),
              col = cols)
text(p4,table(Net_Cov$BCS_cat), labels =table(Net_Cov$BCS_cat), pos = 3)
text(p4,table(Net_Cov$BCS_cat), labels =paste0(round(prop.table(table(Net_Cov$BCS_cat))*100), "%"), pos = 1)


### Weight Plot
par(mfrow = c(1,2))
hist(Net_Cov$weight/100, probability = T, col = cols.2[2],
     xlab = "weight", main = "")
lines(density(Net_Cov$weight/100), col = cols.2[1], lwd=2)

boxplot(Net_Cov$weight/100, pch = 19)
mtext("Weight Distribution", line = -1.5, outer = T)
par(mfrow = c(1,1))


### Parity
p6 <- barplot(table(Net$Net_COV$parity),
              xlab = "Parity",
              ylab = "Frequency",
              main = "Parity Distribution",
              ylim = c(0, 50),
              col = cols)
text(p6,table(Net$Net_COV$parity), labels =table(Net$Net_COV$parity), pos = 3)


### affected quarter
table(sapply(1:nrow(Net$Net_COV), function(x)ifelse(Net$Net_COV$qtr_B[x]=="None", 0, 1)))


### Treatment
par(mfrow=c(2,1))
p5 <- barplot(table(Net$Net_COV$treatment),
              xlab = "Treatment",
              ylab = "Frequency",
              main = "",
              ylim = c(0, 50),
              col = cols)
text(p5,table(Net$Net_COV$treatment), labels =table(Net$Net_COV$treatment), pos = 3)


### Subtreatment
p7 <- barplot(table(Net$Net_COV$subtreatment),
              xlab = "Subtreatment",
              ylab = "Frequency",
              main = "",
              ylim = c(0, 60),
              col = cols)
text(p7,table(Net$Net_COV$subtreatment), labels =table(Net$Net_COV$subtreatment), pos = 3)
par(mfrow=c(1,1))

#########################################################