# this code is used to draw the barplot to show the ranking of different drugs
# Aim1: barplot for the following six drugs (selective serotonin reuptake inhibitor (SSRI))
#   "Paroxetine, Fluvoxamine, Escitalopram, Sertraline, Fluoxetine, Citalopram"
# Aim2: barplot for the following five drugs (serotonin–norepinephrine reuptake inhibitor (SNRI))
#   "Duloxetine, Venlafaxine, Milnacipran, Levomilnacipran, Desvenlafaxine" 
# last update: 05/27/2019

library(mvtnorm)
library(ggplot2)
library(reshape)
library(scales)
source("/Users/Jessie/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis")

#Run analysis
out = CLNMA.equal.tau.fullout(dataout)

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
Nsim = 50000
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

# ggplot
pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_prob_ranking_6drugs_efficacy_50000.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color='black') +
  scale_fill_manual("Drugs", values = c("white","#cfcfcd", "#a7a6a3", "#787878", "#575757","black")) +
  scale_y_continuous(labels = percent_format())+
  xlab("Ranks") + ylab("% probability to rank at each place") +
  ggtitle("Efficacy Only")+theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15)
dev.off()

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


datm <- melt(cbind(C2, ind = rownames(C2)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l1 = c("citalopram","escitalopram", "fluoxetine","fluvoxamine", "paroxetine", "sertraline")
datm$Treatments =mapvalues(datm$Treatment, from =1:6, to=l1)

pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_prob_ranking_6drugs_50+50_50000.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color='black') +
  xlab("Ranks") + ylab("% probability to rank at each place") +
  scale_fill_manual("Drugs", values = c("white","#cfcfcd", "#a7a6a3", "#787878", "#575757","black")) +
  scale_y_continuous(labels = percent_format())+ggtitle("50%Efficacy+50%Safety")+theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15)
dev.off()

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


datm <- melt(cbind(C3, ind = rownames(C3)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l1 = c("citalopram","escitalopram", "fluoxetine","fluvoxamine", "paroxetine", "sertraline")
datm$Treatments =mapvalues(datm$Treatment, from =1:6, to=l1)

pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_prob_ranking_6drugs_safety_50000.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color='black') +
  xlab("Ranks") + ylab("% probability to rank at each place") +
  scale_fill_manual("Drugs", values = c("white","#cfcfcd", "#a7a6a3", "#787878", "#575757","black")) +
  scale_y_continuous(labels = percent_format())+ggtitle("Safety Only")+theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
  theme_classic(base_size = 15)
dev.off()
##########################################################################
########################## done  #########################################
##########################################################################




##########################################################################
########################## five drugs #####################################
##########################################################################
#SUCRA plot
Nsim = 50000
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

# ggplot
pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_prob_ranking_5drugs_efficacy_50000.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color="black") +
  scale_fill_manual("Drugs", values = c("white","#cfcfcd", "#a7a6a3", "#787878", "black")) +
  scale_y_continuous(labels = percent_format())+
  xlab("Ranks") + ylab("% probability to rank at each place") +
  ggtitle("Efficacy Only")+theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15)
dev.off()


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


datm <- melt(cbind(C2, ind = rownames(C2)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l2 = c("desvenlafaxine", "duloxetine", "levomilnacipran", "milnacipran", "venlafaxine")
datm$Treatments =mapvalues(datm$Treatment, from =1:5, to=l2)

# ggplot
pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_prob_ranking_5drugs_50+50_50000.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color="black") +
  scale_fill_manual("Drugs", values = c("white","#cfcfcd", "#a7a6a3", "#787878", "black")) +
  scale_y_continuous(labels = percent_format())+
  xlab("Ranks") + ylab("% probability to rank at each place") +
  ggtitle("50%Efficacy+50%Safety")+theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15)
dev.off()

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


datm <- melt(cbind(C3, ind = rownames(C3)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l2 = c("desvenlafaxine", "duloxetine", "levomilnacipran", "milnacipran", "venlafaxine")
datm$Treatments =mapvalues(datm$Treatment, from =1:5, to=l2)

# ggplot
pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_prob_ranking_5drugs_safety_50000.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color="black") +
  scale_fill_manual("Drugs", values = c("white","#cfcfcd", "#a7a6a3", "#787878","black")) +
  scale_y_continuous(labels = percent_format())+
  xlab("Ranks") + ylab("% probability to rank at each place") +
  ggtitle("Safety Only")+theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15)
dev.off()
##########################################################################
########################## done  #########################################
##########################################################################


