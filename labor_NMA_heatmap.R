# this code is used to draw the correlation heatmap
# written by Jessie
# last update: 08/06/2019

tmp =  diag(1/sqrt(diag(V)))
corr_matrix <- tmp %*% V %*% tmp

# corr_matrix = data.frame(corr_matrix)
corr_matrix_update = cbind(corr_matrix, c(c(1:13),c(1:3,5:13),c(1:3,5:13),c(1:13),c(1:3,5,7:13)))
list_order = order(corr_matrix_update[,62])
updated_corr_matrix_tmp = corr_matrix[list_order,][,-62]
updated_corr_matrix = updated_corr_matrix_tmp[,list_order]

table = table(corr_matrix_update[,62])

# melted_corr_matrix <- melt(corr_matrix)
melted_corr_matrix <- melt(updated_corr_matrix)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/update_labor_heatmap.pdf", width = 10, height = 9)
p <- ggplot(melted_corr_matrix, aes(factor(X1), factor(X2))) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradient(low = "white",high = "steelblue") + 
  theme_classic() +
  labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0,0), labels = c(1:61)) + 
  scale_y_discrete(expand = c(0,0), labels = c(1:61))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=20)) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 5.5) +
  geom_hline(yintercept = 5.5+5) +
  geom_hline(yintercept = 5.5+5+5) +
  geom_hline(yintercept = 5.5+5+5+2) +
  geom_hline(yintercept = 5.5+5+5+2+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5+5+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5+5+5+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5+5+5+5+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5+5+5+5+5+5) +
  geom_hline(yintercept = 5.5+5+5+2+5+4+5+5+5+5+5+5+5) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 5.5) +
  geom_vline(xintercept = 5.5+5) +
  geom_vline(xintercept = 5.5+5+5) +
  geom_vline(xintercept = 5.5+5+5+2) +
  geom_vline(xintercept = 5.5+5+5+2+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5+5+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5+5+5+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5+5+5+5+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5+5+5+5+5+5) +
  geom_vline(xintercept = 5.5+5+5+2+5+4+5+5+5+5+5+5+5)
p
dev.off()
