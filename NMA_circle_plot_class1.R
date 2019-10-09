# this code is used to draw the confidence interval for 21 drugs
# x-axis is the efficacy
# y-axis is the safety 
# last update: 05/27/2019
# load the results from Rui
# load("/Users/Jessie/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis/Results.RData")
load("/Users/jiayito/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis/Results.RData")

library(ggplot2)
library(gridExtra)
library(grid)
library(randomcoloR)

# efficacy
# mu1 is the log odds ratio for efficacy
# exp(mu1) is the odds ratio for efficacy
exp_mu1 = exp(mu1)

# Varmatrix is the covariance matrix for efficacy and safety
# cov(mu1) = matrix(1:21, 1:21), entry on 22x22 is (tau_1)^2
cov_matrix_mu1 = Varmatrix[1:21, 1:21]

# variance of efficacy
var_mu1 = diag(cov_matrix_mu1)

# Confidence interval for efficacy (log odds ratio)
mu1u = mu1+1.96*sqrt(var_mu1)
mu1l = mu1-1.96*sqrt(var_mu1)

# Confidence interval for efficacy (odds ratio)
upper_mu1 = exp(mu1u)
lower_mu1 = exp(mu1l)


# safety
# mu2 is the log odds ratio for safety
# exp(mu2) is the odds ratio for efficacy
exp_mu2 = exp(mu2)

# cov(mu2) = matrix(23:43, 23:43), entry on 44x44 is (tau_2)^2
cov_matrix_mu2 = Varmatrix[23:43, 23:43]

# variance of safety
var_mu2 = diag(cov_matrix_mu2)

#Confidence interval for safety (log odds ratio)
mu2u = mu2+1.96*sqrt(var_mu2)
mu2l = mu2-1.96*sqrt(var_mu2)

# Confidence interval for safety (odds ratio)
upper_mu2 = exp(mu2u)
lower_mu2 = exp(mu2l)

# drug names
Drugs = rownames(table(Data$treatment))[-16]

# class2 name list
class1_list = c(4,8,9,10,15,17)
Drugs = Drugs[class1_list]
exp_mu1 = exp_mu1[class1_list]
exp_mu2 = exp_mu2[class1_list]
lower_mu1 = lower_mu1[class1_list]
lower_mu2 = lower_mu2[class1_list]
upper_mu1 = upper_mu1[class1_list]
upper_mu2 = upper_mu2[class1_list]


# x-axis is efficacy; y-axis is safety
df = data.frame(exp_mu1,exp_mu2, Drugs)

n = 6
set.seed(4321)
palette <- distinctColorPalette(n)

# ggplot
pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/1001_class1_update_eff_safety_2dim.pdf",height=7,width=11)
ggplot(df, aes(x = exp_mu1, y = exp_mu2, fill = Drugs, color = Drugs)) + 
  geom_point(aes(x=1, y=1),color='red',size=5)+
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept=1, linetype="dashed", color = "red") +
  geom_text(label = "placebo", color = "black", aes(x=1.06, y=0.98)) +
  geom_segment(aes(x = lower_mu1, y = exp_mu2, xend = upper_mu1 , yend = exp_mu2)) +
  geom_segment(aes(x = exp_mu1, y = lower_mu2, xend = exp_mu1 , yend = upper_mu2)) +
  geom_point(aes(fill=Drugs,colour = Drugs),size=5) +
  geom_text(label = as.character(1:6), color = "black", size = 3.5) +
  scale_x_continuous(limits = c(1,2.3),name ="Efficacy (Response Rate)", breaks = c(1,1.25,1.5,1.75,2,2.25)) +
  scale_y_reverse(name ="Safety (Dropout Rate)") +
  scale_colour_manual(values = palette)+
  guides(color=guide_legend(ncol=1))+
  theme_classic(base_size = 15)

grid.text("1",
          x = unit(0.872, "npc"), y = unit(0.602, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("2",
          x = unit(0.872, "npc"), y = unit(0.602-0.035, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("3",
          x = unit(0.872, "npc"), y = unit(0.602-2*0.035, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("4",
          x = unit(0.872, "npc"), y = unit(0.602-3*0.035, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("5",
          x = unit(0.872, "npc"), y = unit(0.602-4*0.035, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("6",
          x = unit(0.872, "npc"), y = unit(0.602-5*0.035, "npc"),just = "left",  gp=gpar(fontsize=10))

dev.off() 

# save image
# save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/Results_JT.RData")

