# this code is used to draw the confidence interval labor induction data 
# x-axis is the efficacy
# y-axis is the safety 
# last update: 09/06/2019
# change all names of treatments to lower case letter
# Create a two dimensional visualization using "Caesarean section" and "Serious neonatal morbidity or perinatal death".
# load the results from Rui
# with new list of drugs from Lisa
## #2, #3, #5, #7, #8, #9, #10, #11
# load("/Users/Jessie/Dropbox/2_Rui_and_Jessie/Summer_2019/1_Multi_NMA/Data_analysis/Results.RData")
load("/Users/jiayito/Dropbox/Multi_NMA/2_dataset2/Labor_Induction/Alfirevic2015BMJ_data/Results.RData")

library(ggplot2)
library(gridExtra)
library(grid)
library(randomcoloR)

newlist = c(2,3,5,7,8,9,10,11)
##############################################################
# Caesarean section is data1 in "2th_data_analysis_forestplot.R"
##############################################################
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


##############################################################
# Serious neonatal morbidity or perinatal death is data4 in "2th_data_analysis_forestplot.R"
##############################################################
##################### data4 ################################
l2 = length(mu2)
l3 = length(mu3)
l4 = length(mu4) # 13 drugs

exp_mu4 = exp(mu4)

cov_matrix_mu4 = V[(l1+l2+l3+1):(l1+l2+l3+l4), (l1+l2+l3+1):(l1+l2+l3+l4)]

var_mu4 = diag(cov_matrix_mu4)

mu4u = mu4+1.96*sqrt(var_mu4)
mu4l = mu4-1.96*sqrt(var_mu4)

upper_mu4 = exp(mu4u)
lower_mu4 = exp(mu4l)

#################
# revised by jessie 
# 5 is too closed to 4, then change OR of CS a little bit larger
tmp_list = exp_mu1[newlist][order(exp_mu1[newlist], decreasing = FALSE)]
exp_mu1[which(exp_mu1 == tmp_list[4])] = 0.73

# 1 is too closed to 3 and 7, then change OR of serious neoxxx a little bit 
tmp_list_2 = exp_mu4[newlist][order(exp_mu4[newlist], decreasing = FALSE)]
exp_mu4[which(exp_mu4 == tmp_list_2[1])] = 0.84
exp_mu4[which(exp_mu4 == tmp_list_2[2])] = 0.85
exp_mu4[which(exp_mu4 == tmp_list_2[3])] = 0.857
exp_mu4[which(exp_mu4 == tmp_list_2[4])] = 0.867

###############################################
# drug names
# 13 labels for T1 and T2, omit last one ("placebo")
Drugs_full = c("vaginal PGE2 (tablet)",				
            "vaginal PGE2 (gel)",			
            "vaginal PGE2 pessary (slow release)",			
            "PGF2 gel",			
            "intracervical PGE2",			
            "vaginal PGE2 pessary (normal release)",				
            "vaginal misoprostol (dose less than 50 mcg)",				
            "vaginal misoprostol (dose 50 mcg or more)",				
            "oral misoprostol tablet (dose less than 50 mcg)",				
            "oral misoprostol tablet (dose 50mcg or more)",				
            "titrated (low dose) oral misoprostol solution",				
            "sustained release misoprostol vaginal pessary",	
            "no treatment")

Treatments = Drugs_full[newlist]
# x-axis is efficacy; y-axis is safety
df = data.frame(Treatments,exp_mu1[newlist],exp_mu4[newlist])

# color
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ggplot
pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/1008_reverse_update4_labor_circlr_plot.pdf",
    height=6,width=10)
upper_mu1[newlist][6] = 1.5
ggplot(df, aes(x = exp_mu1[newlist], y = exp_mu4[newlist])) + 
  geom_point(aes(x=1, y=1),size=5, color='red') +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept=1, linetype="dashed", color = "red") + 
  geom_text(label = "placebo", color = "black", aes(x=0.94, y=1.02)) +
  geom_segment(aes(x = lower_mu1[newlist], y = exp_mu4[newlist], xend = upper_mu1[newlist] , yend = exp_mu4[newlist], color = Treatments)) +
  geom_segment(aes(x = exp_mu1[newlist], y = lower_mu4[newlist], xend = exp_mu1[newlist] , yend = upper_mu4[newlist],color = Treatments)) +
  geom_point(aes(color = Treatments),size=5) +
  geom_text(label = as.character(1:8), size = 3.5, color = "black") +
  # scale_x_reverse(limits = c(0.5,1.5),name ="Caesarean section (OR)", breaks = c(0.5,0.75,1.0,1.25,1.5)) +
  scale_x_reverse(name ="Caesarean section (OR)") + 
  scale_y_reverse(name ="Serious neonatal morbidity or perinatal death (OR)") + 
  scale_color_manual(values = cbp1,breaks=c("vaginal PGE2 (gel)",			
                                "vaginal PGE2 pessary (slow release)",	
                                "intracervical PGE2",
                                "vaginal misoprostol (dose less than 50 mcg)",				
                                "vaginal misoprostol (dose 50 mcg or more)",				
                                "oral misoprostol tablet (dose less than 50 mcg)",				
                                "oral misoprostol tablet (dose 50mcg or more)",				
                                "titrated (low dose) oral misoprostol solution"		
                                )) +
  guides(color=guide_legend(ncol=1))+
  theme_classic(base_size = 15) 


grid.text("1",
          x = unit(0.608, "npc"), y = unit(0.658, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("2",
          x = unit(0.608, "npc"), y = unit(0.619, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("3",
          x = unit(0.608, "npc"), y = unit(0.619-0.040, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("4",
          x = unit(0.608, "npc"), y = unit(0.619-2*0.040, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("5",
          x = unit(0.608, "npc"), y = unit(0.619-3*0.040, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("6",
          x = unit(0.608, "npc"), y = unit(0.619-4*0.040, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("7",
          x = unit(0.608, "npc"), y = unit(0.619-5*0.040, "npc"),just = "left",  gp=gpar(fontsize=10))
grid.text("8",
          x = unit(0.608, "npc"), y = unit(0.619-6*0.040, "npc"),just = "left",  gp=gpar(fontsize=10))
# grid.text("9",
#           x = unit(0.707, "npc"), y = unit(0.652-7*0.028, "npc"),just = "left",  gp=gpar(fontsize=10))
# grid.text("10",
#           x = unit(0.704, "npc"), y = unit(0.652-8*0.028, "npc"),just = "left",  gp=gpar(fontsize=10))
# grid.text("11",
#           x = unit(0.704, "npc"), y = unit(0.652-9*0.028, "npc"),just = "left",  gp=gpar(fontsize=10))
# grid.text("12",
#           x = unit(0.704, "npc"), y = unit(0.652-10*0.028, "npc"),just = "left",  gp=gpar(fontsize=10))
# grid.text("13",
#           x = unit(0.704, "npc"), y = unit(0.651-11*0.028, "npc"),just = "left",  gp=gpar(fontsize=10))
dev.off() 

# save image
# save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/Results_JT.RData")

