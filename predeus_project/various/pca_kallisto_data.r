# Make a PCA biplot of kallisto quant output

dire = readline('Where are TSVs?   ')
setname = readline('What are your samples?   ')
kallisto_data = data.frame(c(1:216139))
for (exp in dir(dire)){
    tab = read.table(paste0(dire, exp), header=T)
    tpms = ifelse(tab[,5] > 0, log(tab[,5], 2), 0)
    kallisto_data = cbind(kallisto_data, tpms)
}
kallisto_data = kallisto_data[, c(2:ncol(kallisto_data))]
colnames(kallisto_data) = c(1:ncol(kallisto_data))

pcaout = prcomp(t(as.matrix(kallisto_data)))
samples = factor(sapply(strsplit(dir(dire), '_rep'), function(elt) elt[1]))

pcagraph = ggplot(as.data.frame(pcaout$x), aes(x=PC1, y=PC2, size=20, col=samples)) + 
  geom_point() + 
  scale_size_area(max_size = 15, guide=F) + 
  scale_color_hue(labels=levels(samples), name=setname) +
  guides(colour=guide_legend(override.aes=list(size=10))) +
  xlab(paste0('PC1 - ', summary(pcaout)$importance[2, 1]*100, ' % of variance')) +
  ylab(paste0('PC2 - ', summary(pcaout)$importance[2, 2]*100, ' % of variance')) +
  ggtitle(paste('PCA of', setname)) + theme_bw() +
  theme(plot.title=element_text(face='bold', size=rel(2)), legend.text=element_text(size=18), 
        legend.title=element_text(size=18), axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=18))

print(pcagraph)
