# this code is used to draw the correlation heatmap
# written by Jessie
# last update: 08/06/2019

tmp =  diag(1/sqrt(diag(V)))
corr_matrix <- tmp %*% V %*% tmp

melted_corr_matrix <- melt(corr_matrix)

pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_with_Rui/summer_2019_with_Rui/0_NMA/labor_heatmap.pdf", width = 10, height = 9)
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
  geom_hline(yintercept = 13.5) +
  geom_hline(yintercept = 13.5+12) +
  geom_hline(yintercept = 13.5+12+12) +
  geom_hline(yintercept = 13.5+12+12+13) +
  geom_hline(yintercept = 13.5+12+12+13+11) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 13+0.5) +
  geom_vline(xintercept = 13.5+12) +
  geom_vline(xintercept = 13.5+12+12) +
  geom_vline(xintercept = 13.5+12+12+13) +
  geom_vline(xintercept = 13.5+12+12+13+11)
p
dev.off()
