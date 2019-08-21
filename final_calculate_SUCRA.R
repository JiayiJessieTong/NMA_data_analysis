# this code is to calculate final SUCRA
# set path
mypath = "/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/"
setwd(mypath)

# name list 
l =  c("citalopram","escitalopram", "fluoxetine","fluvoxamine", "paroxetine", "sertraline",
       "desvenlafaxine", "duloxetine", "levomilnacipran", "milnacipran", "venlafaxine")

# read .csv files
# eff
eff1 = read.csv('100_eff.csv')[,-1]
eff2 = read.csv('200_eff.csv')[,-1]
eff3 = read.csv('300_eff.csv')[,-1]
eff4 = read.csv('400_eff.csv')[,-1]
eff5 = read.csv('500_eff.csv')[,-1]

eff_full = rbind(eff1, eff2, eff3, eff4, eff5)
dim(eff_full) 
colnames(eff_full) = l
SUCRA_efficacy = apply(eff_full, 2, mean)

# 50+50
fifty1 = read.csv('100_fifty.csv')[,-1]
fifty2 = read.csv('200_fifty.csv')[,-1]
fifty3 = read.csv('300_fifty.csv')[,-1]
fifty4 = read.csv('400_fifty.csv')[,-1]
fifty5 = read.csv('500_fifty.csv')[,-1]

fifty_full = rbind(fifty1, fifty2, fifty3, fifty4, fifty5)
dim(fifty_full) 
colnames(fifty_full) = l
SUCRA_fifty = apply(fifty_full, 2, mean)

# safety
safety1 = read.csv('100_safety.csv')[,-1]
safety2 = read.csv('200_safety.csv')[,-1]
safety3 = read.csv('300_safety.csv')[,-1]
safety4 = read.csv('400_safety.csv')[,-1]
safety5 = read.csv('500_safety.csv')[,-1]

safety_full = rbind(safety1, safety2, safety3, safety4, safety5)
dim(safety_full) 
colnames(safety_full) = l
SUCRA_safety = apply(safety_full, 2, mean)

final_results = round(rbind(SUCRA_efficacy, SUCRA_fifty, SUCRA_safety), 4)
row.names(final_results) = c("efficaty", "50+50", "safety")

print(final_results)
write.csv(final_results, "SUCRA_results.csv")
