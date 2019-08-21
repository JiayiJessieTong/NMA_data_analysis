# this code is used to draw the barplot to show the ranking of different drugs
# draw bar plot for Create a utility function of the outcome 
# change "Vaginal delivery not achieved within 24 hours" mu2
# to "serious maternal morbidity or death" (mu5)
# and 
# "Caesarean section". mu1
# last update: 08/21/2019

library(mvtnorm)
library(ggplot2)
library(reshape)
library(scales)
library(plyr)
library(randomcoloR)
# source("/Users/jiayito/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis")
# n = 11
# set.seed(1)
# palette <- distinctColorPalette(n)


#Run analysis
out = load("/Users/jiayito/Dropbox/Multi_NMA/2_dataset2/Labor_Induction/Alfirevic2015BMJ_data/Results.RData")

newlist = c(2,3,5,7,8,9,10,11)

Drugs_full = c("Vaginal PGE2 (tablet)",				
               "Vaginal PGE2 (gel)",			
               "Vaginal PGE2 pessary (slow release)",			
               "PGF2 gel",			
               "Intracervical PGE2",			
               "Vaginal PGE2 pessary (normal release)",				
               "Vaginal misoprostol (Dose less than 50 mcg)",				
               "Vaginal misoprostol (Dose 50 mcg or more)",				
               "Oral misoprostol tablet (Dose less than 50 mcg)",				
               "Oral misoprostol tablet (dose 50mcg or more)",				
               "Titrated (low dose) oral misoprostol solution",				
               "Sustained release misoprostol vaginal pessary",	
               "No treatment")

Drugs = Drugs_full[newlist]

Varmatrix = V[c(1:13, 51:61), c(1:13, 51:61)]
# Varmatrix = Varmatrix[-c(4,13,25),-c(4,13,25)]


#Confidence interval
mu1u = mu1+1.96*sqrt(diag(Varmatrix)[1:13])
mu1l = mu1-1.96*sqrt(diag(Varmatrix)[1:13])

# mu2u = mu2+1.96*sqrt(diag(Varmatrix)[14:25])
# mu2l = mu2-1.96*sqrt(diag(Varmatrix)[14:25])

mu5u = mu5+1.96*sqrt(diag(Varmatrix)[14:24])
mu5l = mu5-1.96*sqrt(diag(Varmatrix)[14:24])


##########################################################################
########################## eight drugs #####################################
##########################################################################
#SUCRA plot
Nsim = 50000
drug_id_list = c(2,3,5,7,8,9,10,11)
drug_id_list_2 = c(2,3,4,5,6,7,8,9)
# drug_id_list = c(1,2,3,5,6,7,8,9,10,11,12)
# drug_id_list_2 = c(1,2,3,4,5,6,7,8,9,10,11)
########################### 10% + 90% ################################
m = 0.1*mu1[drug_id_list] + 0.9*mu5[drug_id_list_2]
# Varm = 0.1*0.1*Varmatrix[1:11,1:11]+
#   0.9*0.9*Varmatrix[12:22,12:22]+
#   0.1*0.9*Varmatrix[1:11,12:22]+
#   0.1*0.9*t(Varmatrix[1:11,12:22])
Varm = 0.1*0.1*Varmatrix[newlist,newlist]+
  0.9*0.9*Varmatrix[newlist+13,newlist+13]+
  0.1*0.9*Varmatrix[newlist,newlist+13]+
  0.1*0.9*t(Varmatrix[newlist,newlist+13])
y = rmvnorm(Nsim,m,Varm)
R1 = apply(y,1,function(x){order(x,decreasing = T)})
get.count = function(x){
  ct = rep(0,8)
  t = table(x)
  ct[as.numeric(rownames(t))] = t
  return(ct)
}
C1 = apply(R1,1,get.count)            

C1 = data.frame(C1)
colnames(C1) = as.character(1:8)


#Add an id variable for the filled regions
datm <- melt(cbind(C1, ind = rownames(C1)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l1 = Drugs
datm$Treatments = mapvalues(datm$Treatment, from =1:8, to=l1)

levels(as.factor(datm$Treatments))
# ggplot
pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_2nd_barplot_10+90.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color='black') +
  # scale_fill_manual("Drugs", values = palette) +
  scale_y_continuous(labels = percent_format())+
  xlab("Ranks") + ylab("% probability to rank at each place") +
  ggtitle("90% Vaginal delivery not achieved within 24 hours\n+10% Caesarean section")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15) 
dev.off()

########################### 50% + 50% ################################
m = mu1[drug_id_list]-mu5[drug_id_list]
Varm = Varmatrix[newlist,newlist]+
  Varmatrix[newlist+13,newlist+13]+
  Varmatrix[newlist,newlist+13]+
  t(Varmatrix[newlist,newlist+13])
# Varm = Varmatrix[1:11,1:11]+Varmatrix[12:22,12:22]-
#   Varmatrix[1:11,12:22]-t(Varmatrix[1:11,12:22])
y = rmvnorm(Nsim,m,Varm)
R2 = apply(y,1,function(x){order(x,decreasing = T)})
get.count = function(x){
  ct = rep(0,8)
  t = table(x)
  ct[as.numeric(rownames(t))] = t
  return(ct)
}
C2 = apply(R2,1,get.count)

C2 = data.frame(C2)
colnames(C2) = as.character(1:8)


datm <- melt(cbind(C2, ind = rownames(C2)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l1 = Drugs
datm$Treatments =mapvalues(datm$Treatment, from =1:8, to=l1)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_2nd_barplot_50+50.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color='black') +
  xlab("Ranks") + ylab("% probability to rank at each place") +
  scale_fill_manual("Drugs", values = palette) +
  scale_y_continuous(labels = percent_format())+
  ggtitle("50% Vaginal delivery not achieved within 24 hours\n+50% Caesarean section")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  theme_classic(base_size = 15)
dev.off()

######################### 90% + 10% ##################################
m = 0.9*mu1[drug_id_list] + 0.1*mu5[drug_id_list_2]
# Varm = 0.9*0.9*Varmatrix[1:11,1:11]+
#   0.1*0.1*Varmatrix[12:22,12:22]+
#   0.1*0.9*Varmatrix[1:11,12:22]+
#   0.1*0.9*t(Varmatrix[1:11,12:22])
Varm = 0.9*0.9*Varmatrix[newlist,newlist]+
  0.1*0.1*Varmatrix[newlist+13,newlist+13]+
  0.1*0.9*Varmatrix[newlist,newlist+13]+
  0.1*0.9*t(Varmatrix[newlist,newlist+13])
y = rmvnorm(Nsim,m,Varm)
R3 = apply(y,1,function(x){order(x,decreasing = T)})
get.count = function(x){
  ct = rep(0,8)
  t = table(x)
  ct[as.numeric(rownames(t))] = t
  return(ct)
}
C3 = apply(R3,1,get.count)

C3 = data.frame(C3)
colnames(C3) = as.character(1:8)


datm <- melt(cbind(C3, ind = rownames(C3)), id.vars = c('ind'))
colnames(datm) = c("Treatments", "rank","Probability")
datm$Treatments = as.character(datm$Treatments)
l1 = Drugs
datm$Treatments =mapvalues(datm$Treatment, from =1:8, to=l1)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_2nd_barplot_90+10.pdf",height=7,width=11)
ggplot(datm,aes(x = rank, y = Probability,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity",color='black') +
  xlab("Ranks") + ylab("% probability to rank at each place") +
  scale_fill_manual("Drugs", values =  palette) +
  scale_y_continuous(labels = percent_format())+
  ggtitle("10% Vaginal delivery not achieved within 24 hours\n+90% Caesarean section")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
  theme_classic(base_size = 15)
dev.off()
##########################################################################
########################## done  #########################################
##########################################################################




