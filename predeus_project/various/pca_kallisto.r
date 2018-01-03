# Make a PCA biplot of kallisto quant output

dire = readline('Where are TSVs?   ')
setname = readline('What are your samples?   ')
kallisto_data = data.frame(c(1:59366))
for (exp in dir(dire)){
    tab = read.table(paste0(dire, exp), header=T)
    tpms = ifelse(tab[,3] > 0, log(tab[,3], 2), 0)
    kallisto_data = cbind(kallisto_data, tpms)
}
kallisto_data = kallisto_data[, c(2:ncol(kallisto_data))]
colnames(kallisto_data) = c(1:ncol(kallisto_data))

pcaout = prcomp(t(as.matrix(kallisto_data)))
samples = factor(sapply(strsplit(dir(dire), '_'), function(elt) elt[1]))

# Aesthetics should be set up manually!!
pcagraph = ggplot(as.data.frame(pcaout$x), aes(x=PC1, y=PC2, size=20, col=samples)) + 
  geom_point() + geom_text(aes(label=samples, hjust=0, vjust=1.7, size=16)) + 
  scale_size_area(max_size = 8, guide=F) + 
  scale_color_hue(labels=levels(samples), name=setname) +
  guides(colour=guide_legend(override.aes=list(size=15), ncol=2)) +
  xlab(paste0('PC1 - ', summary(pcaout)$importance[2, 1]*100, ' % of variance')) +
  ylab(paste0('PC2 - ', summary(pcaout)$importance[2, 2]*100, ' % of variance')) +
  ggtitle(paste('PCA of', setname)) + theme_bw() +
  theme(plot.title=element_text(face='bold', size=rel(4)), legend.text=element_text(size=30), 
        legend.title=element_text(size=30), axis.title.x=element_text(size=30), 
        axis.title.y=element_text(size=30))

cairo_pdf(paste0('~/pca_', gsub(" ", "_", setname), '.pdf'), width=30, height=18)
print(pcagraph)
dev.off()
