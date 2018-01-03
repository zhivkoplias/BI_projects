############ Plots for sum of Numreads for Kallisto, Salmon, RSEM

ggplot(DATAFR_kal, aes(Tool, Total_sum))+
  geom_bar(stat = 'identity', aes(fill = Tool))+ 
  guides(colour = guide_legend(override.aes = list(size=14)))+
  xlab('Tool')+
  ylab('Sum')+
  theme_bw()+
  ggtitle('Total sum of quantified reads for Heart_6a')+
  scale_fill_manual(values=c('black','red', 'blue'))+
  theme(axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        panel.background = element_blank(),
        plot.title = element_text(size = 20, colour = "black"))

DATAFR_kal <- data.frame(variable1=rep(c("RSEM", "Salmon", "kallisto")), value1=c(263564593,17255423,15460330))

colnames(DATAFR_kal) <- c('Tool', 'Total sum')

colnames(DATAFR_kal) <- c('Tool', 'Total_sum')


