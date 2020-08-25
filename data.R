########## Libraries ##########

library(birk)
library(DescTools)
library(lubridate)
library(RColorBrewer)
library(slickR)
library(htmlwidgets)
library(corrplot)
library(igraph)
library(statnet)
library(ergm.count)
library(intergraph)
library(kableExtra)
library(GGally)
library(ggnet)

###############################



################### Load Data sets ####################

# Time series data
Net_TS_Feb <- read.csv("Network_TS_Feb.csv")
Net_TS_Mar <- read.csv("Network_TS_March.csv")
Net_TS_Apr <- read.csv("Network_TS_April.csv")
Net_TS_May <- read.csv("Network_TS_May.csv")
Net_TS_Jun <- read.csv("Network_TS_June.csv")
Net_TS_Jul <- read.csv("Network_TS_July.csv")

# Covariate data
Net_COV_Feb <- read.csv("Network_COV_Feb.csv")
Net_COV_Mar <- read.csv("Network_COV_March.csv")
Net_COV_Apr <- read.csv("Network_COV_April.csv")
Net_COV_May <- read.csv("Network_COV_May.csv")
Net_COV_Jun <- read.csv("Network_COV_June.csv")
Net_COV_Jul <- read.csv("Network_COV_July.csv")

######################################################



########################## Data processing functions ##################################

## Remove Cows with NA in covariates of interest
NAs <- function(data){
  return(which(is.na(data$weight) | is.na(data$BCS) | is.na(data$treatment) | is.na(data$subtreatment) 
               | is.na(data$SCC_RF) | is.na(data$SCC_LF) | is.na(data$SCC_RH) | is.na(data$SCC_LH)
               | data$treatment == "" | data$subtreatment == ""))
}


## not-in operator
"%!in%" <-Negate("%in%")

## Remove duplicates 
rm.dup <- function(data){
  if(any(duplicated(data$Cow))){
    return (data[-which(duplicated(data$Cow)),])
  }else{
    return (data)
  }
}


## Update TS dataset to include only cows present in updated COV data set
update.TS <- function(data.TS, data.COV){
  COV_id <- which(unique(data.TS$tb_num) %in% data.COV$Cow) # cows present in both TS and COV
  TS_ind <- which(data.TS$tb_num %in% unique(data.TS$tb_num)[COV_id]) # indexes of cows in TS present in COV
  data.TS <- data.TS[TS_ind,]
  return(data.TS)
}


## Split ```Net_cov``` and ```Net_TS``` into morning and evening milking
milking.time <- function(data.TS){
  M_M_ind <- which(data.TS$milking_time == 1) # M_M_ind : morrning milking indexes
  E_M_ind <- which(data.TS$milking_time == 2) # E_M_ind : evening milking indexes
  return(list("M" = data.TS[M_M_ind,], "E" = data.TS[E_M_ind,]))
}


## Update COV dataset to include only cows present in updated TS data set
update.COV <- function(data.TS, data.COV){
  return(data.COV[which(data.COV$Cow %in% unique(data.TS$tb_num)),])
}

## Remove time series obs < 21
TS.21 <- function(data.TS){
  tab <- as.data.frame(table(data.TS$tb_num))
  B21 <- tab[which(tab$Freq < 21), 1]
  B21_ind <- which(data.TS$tb_num %in% B21)
  if(identical(B21_ind, integer(0))){
    return(data.TS)
  }else{
    return(data.TS[-B21_ind,])
  }
}

## Clear old factor levels from categorical variables
un.factor <- function(data.COV){
  data.COV$Cow <- as.character(data.COV$Cow)
  data.COV$treatment <- as.character(data.COV$treatment)
  data.COV$subtreatment <- as.character(data.COV$subtreatment)
  data.COV$parity <- as.integer(data.COV$parity)
  data.COV$farm <- as.integer(data.COV$farm)
  
  return(data.COV)
}

## Update factor levels
re.factor <- function(data.COV){
  data.COV <- un.factor(data.COV)
  data.COV$Cow <- as.factor(data.COV$Cow)
  data.COV$treatment <- as.factor(data.COV$treatment)
  data.COV$subtreatment <- as.factor(data.COV$subtreatment)
  data.COV$parity <- as.factor(data.COV$parity)
  data.COV$farm <- as.factor(data.COV$farm)
  
  return(data.COV)
}


## TS id factor levels
re.factor.TS <- function(data.TS){
  data.TS$tb_num <- as.character(data.TS$tb_num)
  data.TS$tb_num <- as.factor(data.TS$tb_num)
  return(data.TS)
}

## Time series plots of yields for each cow in ```Net_TS_M/E``` (check for outliers)
ts.plots.carousel <- function(data.TS, data.COV){
  ifelse(!dir.exists("tmp_plots"), dir.create("tmp_plots"), F)
  
  unlink("tmp_plots/*")
  
  cols <- brewer.pal(12, "Paired")
  
  data.TS$tb_num <- as.character(data.TS$tb_num)
  
  i <- 0
  while (i < nrow(data.COV)) {
    
    png(file=paste0("tmp_plots/slide",i,".png"), width = 600)
    par(mfrow = c(3,3))
    tmp1 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+1]) ,]
    tmp2 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+2]) ,]
    tmp3 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+3]) ,]
    tmp4 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+4]) ,]
    tmp5 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+5]) ,]
    tmp6 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+6]) ,]
    tmp7 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+7]) ,]
    tmp8 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+8]) ,]
    tmp9 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+9]) ,]
    
    for (j in 1:9) {
      if(identical(eval(parse(text = paste0("tmp",j)))$tb_num, character(0))){
        break
      }
      
      plot(eval(parse(text = paste0("tmp",j)))$yield, type = "b", pch = 19, col=cols[2],
           main = paste0(i+j," : ",data.COV$Cow[i+j]),
           ylab = "Yield",
           xlab = "Time")
      
    }
    
    dev.off()
    i <- i+9
  }
  
  # Plot carousel
  return(slickR(list.files("tmp_plots",full.names = TRUE,pattern = 'png'), height = "100%", width = "100%")+
           settings(dots = T))
}


## Remove outliers 
rm.outliers <- function(data.TS, data.COV, indexes, obs){
  tmp <- c()
  
  for (i in 1:length(indexes)) {
    cond <- data.TS$tb_num == data.COV$Cow[indexes[i]]
    cond2 <- eval(parse(text = paste0("data.TS$yield ",obs[i])))
    tmp <- c(tmp,which(cond & cond2))
  }
  
  return(data.TS[-tmp,])
}


## Plot loess vs original values
loess.ts.carousel <- function(data.TS, data.COV, span = 0.30){
  ifelse(!dir.exists("tmp_plots"), dir.create("tmp_plots"), F)
  
  unlink("tmp_plots/*")
  
  cols <- brewer.pal(12, "Paired")
  
  data.TS$tb_num <- as.character(data.TS$tb_num)
  data.TS$milk_date <- as.Date(data.TS$milk_date)
  
  i <- 0
  tmp_pred <- c()
  
  while (i < nrow(data.COV)) {
    
    png(file=paste0("tmp_plots/slide",i,".png"), width = 600)
    par(mfrow = c(3,3))
    tmp1 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+1]) ,]
    tmp2 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+2]) ,]
    tmp3 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+3]) ,]
    tmp4 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+4]) ,]
    tmp5 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+5]) ,]
    tmp6 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+6]) ,]
    tmp7 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+7]) ,]
    tmp8 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+8]) ,]
    tmp9 <- data.TS[which(data.TS$tb_num %in% data.COV$Cow[i+9]) ,]
    
    for (j in 1:9) {
      tmp <- eval(parse(text = paste0("tmp",j)))
      
      if(identical(tmp$tb_num, character(0))){
        break
      }
      time <- as.numeric(tmp$milk_date - tmp$milk_date[1])
      seq <- seq(min(time), max(time), length=20)
      yield.lo <- loess(tmp$yield ~ time,span = span)
      tmp_pred <- predict(yield.lo,data.frame(time = seq), se = TRUE)
      
      plot(tmp$yield, type = "b", pch = 19, col=cols[2],
           main = paste0(i+j," : ",data.COV$Cow[i+j]),
           ylab = "Yield",
           xlab = "Time")
      
      points(seq, tmp_pred$fit,  pch = 19, col=cols[6])
      lines(seq, tmp_pred$fit, pch = 19, col=cols[6])
      legend("bottomright", legend = c("Raw data",  "Smooth trend"), fill =c(cols[2], cols[6]))
    }
    
    dev.off()
    i <- i+9
  }
  
  # Plot carousel
  return(slickR(list.files("tmp_plots",full.names = TRUE,pattern = 'png'), height = "100%", width = "100%")+
           settings(dots = T))
}


## Generate network objects
Gen.Net <- function(data.TS, data.COV, span = 0.30,threshold = 0.8){
  data.TS$milk_date <- as.Date(data.TS$milk_date)
  
  tmp_mat <- matrix(NA, nrow = 20, ncol = length(unique(data.TS$tb_num)))
  tmp_pred <- c()
  p=1
  
  for (i in 1:nrow(data.COV)) {
    
    tmp_dat <- data.TS[which(data.TS$tb_num == data.COV$Cow[i]),]
    
    if(nrow(tmp_dat)>0){
      time <- as.numeric(tmp_dat$milk_date-tmp_dat$milk_date[1])
      seq <- seq(min(time), max(time), length=20)
      yield.lo <- loess(tmp_dat$yield ~ time, span = span)
      tmp_pred <- predict(yield.lo,data.frame(time = seq), se = TRUE)
      tmp_mat[,p] <- tmp_pred$fit
      p=p+1
    }
    
  }
  
  
  adjmat_weighted <- cor(tmp_mat)
  diag(adjmat_weighted) <- 0
  
  adjmat_unweighted <- round(adjmat_weighted,1)
  adjmat_unweighted[adjmat_unweighted < threshold] <- 0
  adjmat_unweighted <- round(adjmat_unweighted)
  
  weighted_net <- as.network(adjmat_weighted, loops = F, directed = F, matrix.type = "adjacency")
  network.vertex.names(weighted_net) <- as.character(data.COV$Cow) # set vertex names
  set.vertex.attribute(weighted_net, c("farm", "age", "treatment", "subtreatment", "qtr_B", "num_of_bugs", "type_bugs", "weight", "BCS", "parity"),
                       c(data.COV[, c("farm", "Age", "treatment", "subtreatment", "qtr_B", "num_of_bugs", "type_bugs", "weight", "BCS", "parity")]))
  edge_weight <- as.vector(adjmat_weighted)
  set.edge.value(weighted_net, "edge_weight", edge_weight)
  
  unweighted_net <- as.network(adjmat_unweighted, loops = F, directed = F, matrix.type = "adjacency")
  network.vertex.names(unweighted_net) <- as.character(data.COV$Cow) # set vertex names
  set.vertex.attribute(unweighted_net, c("farm", "age", "treatment", "subtreatment", "qtr_B", "num_of_bugs", "type_bugs", "weight", "BCS", "parity"),
                       c(data.COV[, c("farm", "Age", "treatment", "subtreatment", "qtr_B", "num_of_bugs", "type_bugs", "weight", "BCS", "parity")]))
  
  return(list("weighted_Matrix" = adjmat_weighted, 
              "unweighted_matrix" = adjmat_unweighted, 
              "unweighted_net" = unweighted_net, 
              "weighted_net" = weighted_net, "Net_COV" = data.COV, "Loess_data" = tmp_mat))
}
```




## Separate cows in previous month from cows in current month 
remove.prev.month <- function(newdataTS, newdataCOV, prevdataCOV){
  newdataCOV <- rm.dup(newdataCOV)
  ind <- which(newdataCOV$Cow %in% prevdataCOV$Cow)
  tmp.COV <- newdataCOV[-ind,]
  tmp.TS <- update.TS(newdataTS, tmp.COV)
  return(list("TS" = tmp.TS,"COV" = tmp.COV))  
}

## Get full data on cows from previous month
prev.month <- function(newdataTS, newdataCOV, prevdataCOV){
  newdataCOV <- rm.dup(newdataCOV)
  ind <- which(newdataCOV$Cow %in% prevdataCOV$Cow)
  tmp.COV <- newdataCOV[ind,]
  tmp.TS <- update.TS(newdataTS, tmp.COV)
  print(nrow(tmp.TS))
  tmp.TS.M <- TS.21(milking.time(tmp.TS)$M)
  tmp.TS.E <- TS.21(milking.time(tmp.TS)$E)
  tmp.COV.M <- update.COV(tmp.TS.M, tmp.COV)
  tmp.COV.E <- update.COV(tmp.TS.E, tmp.COV)
  return(list("TS.M" = tmp.TS.M, "TS.E" = tmp.TS.E, "COV.M" = tmp.COV.M, "COV.E" = tmp.COV.E))
}


## Update values of variables missing in current month with values from previous month
replace.NAs <- function(newdataCOV, prevdatCOV){
  
  newdataCOV <- newdataCOV[which(newdataCOV$Cow %in% prevdatCOV$Cow), ]
  prevdatCOV <- prevdatCOV[which(prevdatCOV$Cow %in% newdataCOV$Cow), ]
  
  
  newdataCOV <- un.factor(newdataCOV)
  prevdatCOV <- un.factor(prevdatCOV)
  
  newdataCOV <- newdataCOV[order(newdataCOV$Cow),]
  prevdatCOV <- prevdatCOV[order(prevdatCOV$Cow),]
  
  for (i in 1:nrow(newdataCOV)) {
    
    
    if(as.character(newdataCOV$Cow[i]) == as.character(prevdatCOV$Cow[i])){
      newdataCOV$SCC_RF[i] <- ifelse((is.na(newdataCOV$SCC_RF[i]) | newdataCOV$SCC_RF[i] == ""), prevdatCOV$SCC_RF[i], newdataCOV$SCC_RF[i])
      newdataCOV$SCC_LF[i] <- ifelse((is.na(newdataCOV$SCC_LF[i]) | newdataCOV$SCC_LF[i] == ""), prevdatCOV$SCC_LF[i], newdataCOV$SCC_LF[i])
      newdataCOV$SCC_RH[i] <- ifelse((is.na(newdataCOV$SCC_RH[i]) | newdataCOV$SCC_RH[i] == ""), prevdatCOV$SCC_RH[i], newdataCOV$SCC_RH[i])
      newdataCOV$SCC_LH[i] <- ifelse((is.na(newdataCOV$SCC_LH[i]) | newdataCOV$SCC_LH[i] == ""), prevdatCOV$SCC_LH[i], newdataCOV$SCC_LH[i])
      
      newdataCOV$farm[i] <- ifelse((is.na(newdataCOV$farm[i]) | newdataCOV$farm[i] == ""), prevdatCOV$farm[i], newdataCOV$farm[i])
      
      newdataCOV$treatment[i] <- ifelse((is.na(newdataCOV$treatment[i]) | newdataCOV$treatment[i] == ""), prevdatCOV$treatment[i], newdataCOV$treatment[i])
      
      newdataCOV$subtreatment[i] <- ifelse((is.na(newdataCOV$subtreatment[i]) | newdataCOV$subtreatment[i] == ""), prevdatCOV$subtreatment[i], newdataCOV$subtreatment[i])
      
      #newdataCOV$qtr_B[i] <- ifelse((is.na(newdataCOV$qtr_B[i]) | newdataCOV$qtr_B[i] == ""), prevdatCOV$qtr_B[i], newdataCOV$qtr_B[i])
      
      #newdataCOV$num_of_bugs[i] <- ifelse((is.na(newdataCOV$num_of_bugs[i]) | newdataCOV$num_of_bugs[i] == ""), prevdatCOV$num_of_bugs[i], newdataCOV$num_of_bugs[i])
      
      #newdataCOV$type_bugs[i] <- ifelse((is.na(newdataCOV$type_bugs[i]) | newdataCOV$type_bugs[i] == ""), prevdatCOV$type_bugs[i], newdataCOV$type_bugs[i])
      
      newdataCOV$weight[i] <- ifelse((is.na(newdataCOV$weight[i]) | newdataCOV$weight[i] == ""), prevdatCOV$weight[i], newdataCOV$weight[i])
      
      newdataCOV$BCS[i] <- ifelse((is.na(newdataCOV$BCS[i]) | newdataCOV$BCS[i] == ""), prevdatCOV$BCS[i], newdataCOV$BCS[i])
    }else{
      stop("Cow id should be same at all iteration")
    }
  }
  
  return(newdataCOV)
}

########################################################################################




#################### February Network ####################

Net_COV_Feb.2 <- Net_COV_Feb[-NAs(Net_COV_Feb),]
Net_COV_Feb.2 <- rm.dup(Net_COV_Feb.2)
Net_TS_Feb.2 <- update.TS(Net_TS_Feb, Net_COV_Feb.2)

Net_TS_Feb.M <- TS.21(milking.time(Net_TS_Feb.2)$M)
Net_TS_Feb.E <- TS.21(milking.time(Net_TS_Feb.2)$E)

Net_COV_Feb.M <- update.COV(Net_TS_Feb.M, Net_COV_Feb.2)
Net_COV_Feb.E <- update.COV(Net_TS_Feb.E, Net_COV_Feb.2)

# morning
Net_COV_Feb.M <- re.factor(Net_COV_Feb.M)
Net_TS_Feb.M <- re.factor.TS(Net_TS_Feb.M)
ts.plots.carousel(Net_TS_Feb.M, Net_COV_Feb.M)
Net_TS_Feb.M <- rm.outliers(Net_TS_Feb.M, Net_COV_Feb.M, c(104, 105, 115, 92, 45), c("<5", ">16", ">30", ">25", "<2"))
loess.ts.carousel(Net_TS_Feb.M, Net_COV_Feb.M)

Net_M_Feb <- Gen.Net(Net_TS_Feb.M, Net_COV_Feb.M)

save(Net_M_Feb, file = "Net_M_Feb.RData")

# Evening
Net_COV_Feb.E <- re.factor(Net_COV_Feb.E)
Net_TS_Feb.E <- re.factor.TS(Net_TS_Feb.E)
ts.plots.carousel(Net_TS_Feb.E, Net_COV_Feb.E)
Net_TS_Feb.M <- rm.outliers(Net_TS_Feb.E, Net_COV_Feb.E, c(119, 26, 63, 103, 108), c("<4", ">15", ">12", "<2", ">14"))
loess.ts.carousel(Net_TS_Feb.E, Net_COV_Feb.E)


Net_E_Feb <- Gen.Net(Net_TS_Feb.E, Net_COV_Feb.E)
save(Net_E_Feb, file = "Net_E_Feb.RData")

###########################################################




#################### March Network (Only cows present in February) ####################


feb <- prev.month(Net_TS_Mar,Net_COV_Mar, unique(rbind(Net_COV_Feb.M, Net_COV_Feb.E)))

# update missing values with values from feb
feb_COV_M <- replace.NAs(feb$COV.M, unique(rbind(Net_COV_Feb.M, Net_COV_Feb.E)))
feb_COV_E <- replace.NAs(feb$COV.E, unique(rbind(Net_COV_Feb.M, Net_COV_Feb.E)))

# morning
feb_TS_M <- feb$TS.M
feb_COV_M <- re.factor(feb_COV_M)
feb_TS_M <- update.TS(feb_TS_M, feb_COV_M)
feb_TS_M <- re.factor.TS(feb_TS_M)
ts.plots.carousel(feb_TS_M, feb_COV_M)
feb_TS_M  <- rm.outliers(feb_TS_M, feb_COV_M, c(8, 131, 45, 61, 57, 65, 83, 92, 104), c("<9", "<9", ">25", ">21", "<15", "<12", ">25", "<13", "<1"))
loess.ts.carousel(feb_TS_M , feb_COV_M)
Net_M_FebnMar <- Gen.Net(feb_TS_M , feb_COV_M)
save(Net_M_FebnMar, file = "Net_M_FebnMar.RData")

# evening
feb_TS_E <- feb$TS.E

feb_COV_E <- re.factor(feb_COV_E)
feb_TS_E <- update.TS(feb_TS_E, feb_COV_E)
feb_TS_E <- re.factor.TS(feb_TS_E)
ts.plots.carousel(feb_TS_E, feb_COV_E)
feb_TS_E  <- rm.outliers(feb_TS_E, feb_COV_E, c(66), c(">15"))
loess.ts.carousel(feb_TS_E , feb_COV_E)
Net_E_FebnMar <- Gen.Net(feb_TS_E , feb_COV_E)
save(Net_E_FebnMar, file = "Net_E_FebnMar.RData")


#######################################################################################




#################### April Network (Only cows present in February) ####################


feb <- prev.month(Net_TS_Apr,Net_COV_Apr, unique(rbind(Net_COV_Feb.M, Net_COV_Feb.E)))

# update missing values with values from feb
feb_COV_M <- replace.NAs(feb$COV.M, unique(rbind(Net_COV_Feb.M, Net_COV_Feb.E)))
feb_COV_E <- replace.NAs(feb$COV.E, unique(rbind(Net_COV_Feb.M, Net_COV_Feb.E)))

feb_TS_M <- feb$TS.M

feb_COV_M <- re.factor(feb_COV_M)
feb_TS_M <- update.TS(feb_TS_M, feb_COV_M)
feb_TS_M <- re.factor.TS(feb_TS_M)
ts.plots.carousel(feb_TS_M, feb_COV_M)
feb_TS_M  <- rm.outliers(feb_TS_M, feb_COV_M, c(1,5,112,117,119,123,128,44,58,83,89), c("<4", ">16", ">14", ">25", "<5", "<10", "<5", ">30", "<10", "<15", "<5"))
loess.ts.carousel(feb_TS_M , feb_COV_M)
Net_M_FebnApr <- Gen.Net(feb_TS_M , feb_COV_M)
save(Net_M_FebnApr, file = "Net_M_FebnApr.RData")

# evening
feb_TS_E <- feb$TS.E

feb_COV_E <- re.factor(feb_COV_E)
feb_TS_E <- update.TS(feb_TS_E, feb_COV_E)
feb_TS_E <- re.factor.TS(feb_TS_E)
ts.plots.carousel(feb_TS_E, feb_COV_E)
feb_TS_E  <- rm.outliers(feb_TS_E, feb_COV_E, c(116,122,27,22,40,43,50,60,71,84,87,83,15), c("<2", ">15", "<2", ">12", "<4", ">14", "<2", "<2",">20", ">11", ">20", "<2", ">16"))
loess.ts.carousel(feb_TS_E , feb_COV_E)
Net_E_FebnApr <- Gen.Net(feb_TS_E , feb_COV_E)
save(Net_E_FebnApr, file = "Net_E_FebnApr.RData")


#######################################################################################
