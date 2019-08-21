# This code is used to plot the results for NMA project with Rui, Lifeng, and Yong
# Late update by Jessie Tong, July 11

library(forestplot)

# load the results by Rui
load("/Users/jiayito/Dropbox/Multi_NMA/2_dataset2/Labor_Induction/Alfirevic2015BMJ_data/Results.RData")

# 13 labels for T1 and T2, omit last one ("placebo")
label = c("1.Vaginal PGE2 (tablet)",				
          "2.Vaginal PGE2 (gel)",			
          "3.Vaginal PGE2 pessary  (slow release)",			
          "4.PGF2 gel",			
          "5.Intracervical PGE2",			
          "6.Vaginal PGE2 pessary (normal release)",				
          "7.Vaginal misoprostol (Dose less than 50 mcg)",				
          "8.Vaginal misoprostol (Dose 50 mcg or more)",				
          "9.Oral misoprostol tablet (Dose less than 50 mcg)",				
          "10.Oral misoprostol tablet (dose 50mcg or more)",				
          "11.Titrated (low dose) oral misoprostol solution",				
          "12.Sustained release misoprostol vaginal pessary",	
          "13.No treatment")
##################### data1 ################################
l1 = length(mu1) # 13 drugs

# mu1 is the log odds ratio for data1
# exp(mu1) is the odds ratio for data1
exp_mu1 = exp(mu1)

# Varmatrix is the covariance matrix for data 1~5
# cov(mu1) = matrix(1:21, 1:21), entry on 22x22 is (tau_1)^2
cov_matrix_mu1 = V[1:l1, 1:l1]

# variance of efficacy
var_mu1 = diag(cov_matrix_mu1)

# Confidence interval for efficacy (log odds ratio)
mu1u = mu1+1.96*sqrt(var_mu1)
mu1l = mu1-1.96*sqrt(var_mu1)

# Confidence interval for efficacy (odds ratio)
upper_mu1 = exp(mu1u)
lower_mu1 = exp(mu1l)

# forestplot
CI = NULL
for (i in 1:length(mu1)){
  CI[i] = paste(round(exp_mu1[i],2),"(",round(lower_mu1[i],2),",",round(upper_mu1[i],2),")")
}

dataplot <- structure(list(
  mean  = c(NA,NA,exp_mu1), 
  lower = c(NA,NA,lower_mu1),
  upper = c(NA,NA,upper_mu1)),
  .Names = c("mean", "lower", "upper"), 
  row.names = c(NA, -15L), 
  class = "data.frame")

tabletext <- cbind(
  c('','',label),
  c("","OR (95% CI)",CI)
)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/table10.pdf",
    height = 8, width = 11)
forestplot(tabletext, 
           dataplot,new_page = TRUE,
           is.summary=c(TRUE,TRUE,rep(FALSE,13)),
           clip=c(-1,2.5), 
           zero =1,
           graphwidth = unit(3, "inches"),
           xlog=FALSE, 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title = "Caesarean section")
dev.off()


##################### data2 ################################
l2 = length(mu2) # 12 drugs

exp_mu2 = exp(mu2)

cov_matrix_mu2 = V[(l1+1):(l1+l2), (l1+1):(l1+l2)]

var_mu2 = diag(cov_matrix_mu2)

mu2u = mu2+1.96*sqrt(var_mu2)
mu2l = mu2-1.96*sqrt(var_mu2)

upper_mu2 = exp(mu2u)
lower_mu2 = exp(mu2l)

# forestplot
CI2 = NULL
for (i in 1:length(mu2)){
  CI2[i] = paste(round(exp_mu2[i],2),"(",round(lower_mu2[i],2),",",round(upper_mu2[i],2),")")
}

dataplot <- structure(list(
  mean  = c(NA,NA,exp_mu2), 
  lower = c(NA,NA,lower_mu2),
  upper = c(NA,NA,upper_mu2)),
  .Names = c("mean", "lower", "upper"), 
  row.names = c(NA, -14L), 
  class = "data.frame")

label2 = label[-4]

tabletext <- cbind(
  c('','',label2),
  c("","OR (95% CI)",CI2)
)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/table11.pdf",
    height = 8, width = 11)
forestplot(tabletext, 
           dataplot,new_page = TRUE,
           is.summary=c(TRUE,TRUE,rep(FALSE,l2)),
           clip=c(-1,2.5), 
           zero =1,
           graphwidth = unit(3, "inches"),
           xlog=FALSE, 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title = "Vaginal delivery not within 24 hours")
dev.off()


##################### data3 ################################
l3 = length(mu3) # 13 drugs

exp_mu3 = exp(mu3)

cov_matrix_mu3 = V[(l1+l2+1):(l1+l2+l3), (l1+l2+1):(l1+l2+l3)]

var_mu3 = diag(cov_matrix_mu3)

mu3u = mu3+1.96*sqrt(var_mu3)
mu3l = mu3-1.96*sqrt(var_mu3)

upper_mu3 = exp(mu3u)
lower_mu3 = exp(mu3l)

# forestplot
CI3 = NULL
for (i in 1:length(mu3)){
  CI3[i] = paste(round(exp_mu3[i],2),"(",round(lower_mu3[i],2),",",round(upper_mu3[i],2),")")
}

dataplot <- structure(list(
  mean  = c(NA,NA,exp_mu3), 
  lower = c(NA,NA,lower_mu3),
  upper = c(NA,NA,upper_mu3)),
  .Names = c("mean", "lower", "upper"), 
  row.names = c(NA, -14L), 
  class = "data.frame")

label3 = label[-4]

tabletext <- cbind(
  c('','',label3),
  c("","OR (95% CI)",CI3)
)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/table12.pdf",
    height = 8, width = 11)
forestplot(tabletext, 
           dataplot,new_page = TRUE,
           is.summary=c(TRUE,TRUE,rep(FALSE,l3)),
           clip=c(-1,2.5), 
           zero =1,
           graphwidth = unit(3, "inches"),
           xlog=FALSE, 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title = "Hyperstimulation (with fetal heart rate changes)")
dev.off()


##################### data4 ################################
l4 = length(mu4) # 13 drugs

exp_mu4 = exp(mu4)

cov_matrix_mu4 = V[(l1+l2+l3+1):(l1+l2+l3+l4), (l1+l2+l3+1):(l1+l2+l3+l4)]

var_mu4 = diag(cov_matrix_mu4)

mu4u = mu4+1.96*sqrt(var_mu4)
mu4l = mu4-1.96*sqrt(var_mu4)

upper_mu4 = exp(mu4u)
lower_mu4 = exp(mu4l)

# forestplot
CI4 = NULL
for (i in 1:length(mu4)){
  CI4[i] = paste(round(exp_mu4[i],2),"(",round(lower_mu4[i],2),",",round(upper_mu4[i],2),")")
}

dataplot <- structure(list(
  mean  = c(NA,NA,exp_mu4), 
  lower = c(NA,NA,lower_mu4),
  upper = c(NA,NA,upper_mu4)),
  .Names = c("mean", "lower", "upper"), 
  row.names = c(NA, -15L), 
  class = "data.frame")

label4 = label

tabletext <- cbind(
  c('','',label4),
  c("","OR (95% CI)",CI4)
)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/table13.pdf",
    height = 8, width = 11)
forestplot(tabletext, 
           dataplot,new_page = TRUE,
           is.summary=c(TRUE,TRUE,rep(FALSE,l4)),
           clip=c(-1,2.5), 
           zero =1,
           graphwidth = unit(3, "inches"),
           xlog=FALSE, 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title = "Serious neonatal morbidity or perinatal death")
dev.off()

##################### data5 ################################
l5 = length(mu5) # 11 drugs

exp_mu5 = exp(mu5)

cov_matrix_mu5 = V[(l1+l2+l3+l4+1):(l1+l2+l3+l4+l5), 
                   (l1+l2+l3+l4+1):(l1+l2+l3+l4+l5)]

var_mu5 = diag(cov_matrix_mu5)

mu5u = mu5+1.96*sqrt(var_mu5)
mu5l = mu5-1.96*sqrt(var_mu5)

upper_mu5 = exp(mu5u)
lower_mu5 = exp(mu5l)

# forestplot
CI5 = NULL
for (i in 1:length(mu5)){
  CI5[i] = paste(round(exp_mu5[i],2),"(",round(lower_mu5[i],2),",",round(upper_mu5[i],2),")")
}

dataplot <- structure(list(
  mean  = c(NA,NA,exp_mu5), 
  lower = c(NA,NA,lower_mu5),
  upper = c(NA,NA,upper_mu5)),
  .Names = c("mean", "lower", "upper"), 
  row.names = c(NA, -13L), 
  class = "data.frame")

label5 = label[-c(4,6)]

tabletext <- cbind(
  c('','',label5),
  c("","OR (95% CI)",CI5)
)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/table14.pdf",
    height = 8, width = 11)
forestplot(tabletext, 
           dataplot,new_page = TRUE,
           is.summary=c(TRUE,TRUE,rep(FALSE,l5)),
           clip=c(-1,2.5), 
           zero =1,
           graphwidth = unit(3, "inches"),
           xlog=FALSE, 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title = "Serious maternal morbidity or death")
dev.off()


