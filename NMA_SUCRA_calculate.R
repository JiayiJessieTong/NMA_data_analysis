# This code is used to calculate SUCRA (area under the cumulative probability line)
# last update: 06/01/2019

library(mvtnorm)
library(ggplot2)
library(reshape)
library(scales)
library(dplyr)
library(plyr)
library(ggpubr)
load("/Users/Jessie/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis/Results.RData")
source("/Users/Jessie/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis")

# calculate SUCRA function
cal_sucra <- function(datm, l){
  sucra = c()
  for (i in c(1:length(l))){
    drug = l[i]
    dt = datm[which(datm$Treatments == drug),]
    dt$cum_prob = cumsum(dt$Probability)/max(cumsum(dt$Probability))
    sucra[i] <- sum(dt$cum_prob)/length(dt$cum_prob)
  }
  return(sucra)
}

result_matrix_eff = matrix(0, ncol = 11, nrow = 100)
result_matrix_fifty = matrix(0, ncol = 11, nrow = 100)
result_matrix_safety = matrix(0, ncol = 11, nrow = 100)

for (rep in c(1:100)){
  print(paste("The loop number is", rep, "out of 100"))
  
  # this step is the boostrap setp 
  # Aim of this step is to get the point estimates of sucra for several time and then calcualte the mean
  # dataout: raw data
  # number of unique studies in dataout
  num_study = length(unique(dataout$ID)) 
  
  # sample with replacement 
  sample_study = sample(unique(dataout$ID), num_study, replace = T)
  
  # number of unique study in sample_study
  # length(unique(sample_study)) 
  
  # get the new dataset
  newdata = dataout[which(dataout$ID == sample_study[1]),]
  for (i in c(2:length(sample_study))){
    tmp = dataout[which(dataout$ID == sample_study[i]),]
    newdata = rbind(newdata, tmp)
  }
  
  #check uniqueness of t1 and t2
  print(paste("t1:", length(unique(newdata$t1))))
  print(paste("t2:", length(unique(newdata$t2))))
  
  #Run analysis
  out = CLNMA.equal.tau.fullout(newdata)
  
  mu1 = out[[1]]
  mu2 = out[[2]]
  tau1 = out[[3]]
  tau2 = out[[4]]
  Varmatrix =out[[5]]
  
  #Confidence interval
  mu1u = mu1+1.96*sqrt(diag(Varmatrix)[1:21])
  mu1l = mu1-1.96*sqrt(diag(Varmatrix)[1:21])
  
  mu2u = mu2+1.96*sqrt(diag(Varmatrix)[23:43])
  mu2l = mu2-1.96*sqrt(diag(Varmatrix)[23:43])
  
  
  ##########################################################################
  ########################## six drugs #####################################
  ##########################################################################
  #SUCRA plot
  Nsim = 10000
  drug_id_list = c(4, 8, 9, 10, 15, 17)
  ########################### Efficacy only ################################
  m = mu1[drug_id_list]
  Varm = Varmatrix[drug_id_list,drug_id_list]
  y = rmvnorm(Nsim,m,Varm)
  R1 = apply(y,1,function(x){order(x,decreasing = T)})
  get.count = function(x){
    ct = rep(0,6)
    t = table(x)
    ct[as.numeric(rownames(t))] = t
    return(ct)
  }
  C1 = apply(R1,1,get.count)
  C1 = data.frame(C1)
  colnames(C1) = as.character(1:6)
  
  #Add an id variable for the filled regions
  datm <- melt(cbind(C1, ind = rownames(C1)), id.vars = c('ind'))
  colnames(datm) = c("Treatments", "rank","Probability")
  datm$Treatments = as.character(datm$Treatments)
  l1 = c("citalopram","escitalopram", "fluoxetine","fluvoxamine", "paroxetine", "sertraline")
  datm$Treatments =mapvalues(datm$Treatment, from =1:6, to=l1)
  
  # calculate sucra
  result_6_eff <- cal_sucra(datm, l1)
  
  ########################### Efficacy 50% +Safety 50% ################################
  m = mu1[drug_id_list]-mu2[drug_id_list]
  Varm = Varmatrix[drug_id_list,drug_id_list]+Varmatrix[drug_id_list+22,drug_id_list+22]-
    Varmatrix[drug_id_list,drug_id_list+22]-t(Varmatrix[drug_id_list,drug_id_list+22])
  y = rmvnorm(Nsim,m,Varm)
  R2 = apply(y,1,function(x){order(x,decreasing = T)})
  get.count = function(x){
    ct = rep(0,6)
    t = table(x)
    ct[as.numeric(rownames(t))] = t
    return(ct)
  }
  C2 = apply(R2,1,get.count)
  C2 = data.frame(C2)
  colnames(C2) = as.character(1:6)
  
  #Add an id variable for the filled regions
  datm <- melt(cbind(C2, ind = rownames(C2)), id.vars = c('ind'))
  colnames(datm) = c("Treatments", "rank","Probability")
  datm$Treatments = as.character(datm$Treatments)
  l1 = c("citalopram","escitalopram", "fluoxetine","fluvoxamine", "paroxetine", "sertraline")
  datm$Treatments =mapvalues(datm$Treatment, from =1:6, to=l1)
  
  
  # calculate sucra
  result_6_fifty <- cal_sucra(datm, l1)
  
  ######################### Safety Only ##################################
  m = -mu2[drug_id_list]
  Varm =Varmatrix[22+drug_id_list,22+drug_id_list]
  y = rmvnorm(Nsim,m,Varm)
  R3 = apply(y,1,function(x){order(x,decreasing = T)})
  get.count = function(x){
    ct = rep(0,6)
    t = table(x)
    ct[as.numeric(rownames(t))] = t
    return(ct)
  }
  C3 = apply(R3,1,get.count)
  C3 = data.frame(C3)
  colnames(C3) = as.character(1:6)
  
  #Add an id variable for the filled regions
  datm <- melt(cbind(C3, ind = rownames(C3)), id.vars = c('ind'))
  colnames(datm) = c("Treatments", "rank","Probability")
  datm$Treatments = as.character(datm$Treatments)
  l1 = c("citalopram","escitalopram", "fluoxetine","fluvoxamine", "paroxetine", "sertraline")
  datm$Treatments =mapvalues(datm$Treatment, from =1:6, to=l1)
  
  
  # calculate sucra
  result_6_safety <- cal_sucra(datm, l1)
  
  ##########################################################################
  ########################## five drugs #####################################
  ##########################################################################
  #SUCRA plot
  drug_id_list2 = c(6, 7, 11, 12, 19)
  ########################### Efficacy only ################################
  m = mu1[drug_id_list2]
  Varm = Varmatrix[drug_id_list2,drug_id_list2]
  y = rmvnorm(Nsim,m,Varm)
  R1 = apply(y,1,function(x){order(x,decreasing = T)})
  get.count = function(x){
    ct = rep(0,5)
    t = table(x)
    ct[as.numeric(rownames(t))] = t
    return(ct)
  }
  C1 = apply(R1,1,get.count)
  C1 = data.frame(C1)
  colnames(C1) = as.character(1:5)
  
  
  #Add an id variable for the filled regions
  datm <- melt(cbind(C1, ind = rownames(C1)), id.vars = c('ind'))
  colnames(datm) = c("Treatments", "rank","Probability")
  datm$Treatments = as.character(datm$Treatments)
  l2 = c("desvenlafaxine", "duloxetine", "levomilnacipran", "milnacipran", "venlafaxine")
  datm$Treatments =mapvalues(datm$Treatment, from =1:5, to=l2)
  
  # calculate sucra
  result_5_eff <- cal_sucra(datm, l2)
  
  ########################### Efficacy 50% +Safety 50% ################################
  m = mu1[drug_id_list2]-mu2[drug_id_list2]
  Varm = Varmatrix[drug_id_list2,drug_id_list2]+Varmatrix[drug_id_list2+22,drug_id_list2+22]-
    Varmatrix[drug_id_list2,drug_id_list2+22]-t(Varmatrix[drug_id_list2,drug_id_list2+22])
  y = rmvnorm(Nsim,m,Varm)
  R2 = apply(y,1,function(x){order(x,decreasing = T)})
  get.count = function(x){
    ct = rep(0,5)
    t = table(x)
    ct[as.numeric(rownames(t))] = t
    return(ct)
  }
  C2 = apply(R2,1,get.count)
  C2 = data.frame(C2)
  colnames(C2) = as.character(1:5)
  
  #Add an id variable for the filled regions
  datm <- melt(cbind(C2, ind = rownames(C2)), id.vars = c('ind'))
  colnames(datm) = c("Treatments", "rank","Probability")
  datm$Treatments = as.character(datm$Treatments)
  l2 = c("desvenlafaxine", "duloxetine", "levomilnacipran", "milnacipran", "venlafaxine")
  datm$Treatments =mapvalues(datm$Treatment, from =1:5, to=l2)
  
  # calculate sucra
  result_5_fifty <- cal_sucra(datm, l2)
  
  ######################### Safety Only ##################################
  m = -mu2[drug_id_list2]
  Varm =Varmatrix[22+drug_id_list2,22+drug_id_list2]
  y = rmvnorm(Nsim,m,Varm)
  R3 = apply(y,1,function(x){order(x,decreasing = T)})
  get.count = function(x){
    ct = rep(0,5)
    t = table(x)
    ct[as.numeric(rownames(t))] = t
    return(ct)
  }
  C3 = apply(R3,1,get.count)
  C3 = data.frame(C3)
  colnames(C3) = as.character(1:5)
  
  #Add an id variable for the filled regions
  datm <- melt(cbind(C3, ind = rownames(C3)), id.vars = c('ind'))
  colnames(datm) = c("Treatments", "rank","Probability")
  datm$Treatments = as.character(datm$Treatments)
  l2 = c("desvenlafaxine", "duloxetine", "levomilnacipran", "milnacipran", "venlafaxine")
  datm$Treatments =mapvalues(datm$Treatment, from =1:5, to=l2)
  
  # calcualte sucra
  result_5_safety <- cal_sucra(datm, l2)
  
  # assign the result to result matrix
  result_matrix_eff[rep, c(1:6)] = result_6_eff
  result_matrix_eff[rep, c(7:11)] = result_5_eff
   
  result_matrix_fifty[rep, c(1:6)] = result_6_fifty
  result_matrix_fifty[rep, c(7:11)] = result_5_fifty
  
  result_matrix_safety[rep, c(1:6)] = result_6_safety
  result_matrix_safety[rep, c(7:11)] = result_5_safety
  
}

write.csv(result_matrix_eff, "100_eff.csv")
write.csv(result_matrix_fifty, "100_fifty.csv")
write.csv(result_matrix_safety, "100_safety.csv")

save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/SUCRA_Results_JT_1.RData")



