################# TPM PRC
View(USArrests)
H1299 = read.table('DRR016697/quant.sf', header = TRUE)
IMR90_PD32.1 = read.table('SRR1051229/quant.sf', header = TRUE)
IMR90_PD32.2 = read.table('SRR1051230/quant.sf', header = TRUE)
IMR90_PD88.1 = read.table('SRR1051231/quant.sf', header = TRUE)
IMR90_PD88.2 = read.table('SRR1051232/quant.sf', header = TRUE)
A431_replicate_a = read.table('SRR629557/quant.sf', header = TRUE)
A431_replicate_b = read.table('SRR629559/quant.sf', header = TRUE)
U251_MG_replicate_a = read.table('SRR629561/quant.sf', header = TRUE)
U251_MG_replicate_b = read.table('SRR629562/quant.sf', header = TRUE)
U2_OS_replicate_a = read.table('SRR629563/quant.sf', header = TRUE)
A549_replicate_a = read.table('SRR629564/quant.sf', header = TRUE)
A549_replicate_b = read.table('SRR629565/quant.sf', header = TRUE)
CACO2_replicate_a = read.table('SRR629566/quant.sf', header = TRUE)
CACO2_replicate_b = read.table('SRR629568/quant.sf', header = TRUE)
HEK293_replicate_a = read.table('SRR629569/quant.sf', header = TRUE)
HEK293_replicate_b = read.table('SRR629570/quant.sf', header = TRUE)
HeLa_replicate_a = read.table('SRR629571/quant.sf', header = TRUE)
HeLa_replicate_b = read.table('SRR629572/quant.sf', header = TRUE)
HepG2_replicate_a = read.table('SRR629573/quant.sf', header = TRUE)
HepG2_replicate_b = read.table('SRR629576/quant.sf', header = TRUE)
MCF7_replicate_a = read.table('SRR629577/quant.sf', header = TRUE)
MCF7_replicate_b = read.table('SRR629578/quant.sf', header = TRUE)
PC3_replicate_a = read.table('SRR629580/quant.sf', header = TRUE)
# log(IMR90_PD32.1$TPM, 2), log(IMR90_PD32.2$TPM, 2), log(IMR90_PD88.1$TPM, 2), log(IMR90_PD88.2$TPM, 2)
# frame = cbind(H1299$TPM, beigerep2$TPM, brownrep1$TPM, brownrep2$TPM, whiterep1$TPM, whiterep2$TPM)
log_frame = cbind(log(H1299$TPM, 2), log(A431_replicate_a$TPM, 2),
                  log(A431_replicate_b$TPM, 2), log(U251_MG_replicate_a$TPM, 2), log(U251_MG_replicate_b$TPM, 2), log(U2_OS_replicate_a$TPM, 2),
                  log(A549_replicate_a$TPM, 2), log(A549_replicate_b$TPM, 2), log(CACO2_replicate_a$TPM, 2), log(CACO2_replicate_b$TPM, 2),
                  log(HEK293_replicate_a$TPM, 2), log(HEK293_replicate_b$TPM, 2), log(HeLa_replicate_a$TPM, 2), log(HeLa_replicate_b$TPM, 2),
                  log(MCF7_replicate_a$TPM, 2), log(MCF7_replicate_b$TPM, 2), log(PC3_replicate_a$TPM, 2))

is.na(log_frame) <- do.call(cbind,lapply(log_frame, is.infinite))
clear_frame <- na.omit(log_frame)

mtab <- data.frame(clear_frame)
topca = as.matrix(mtab)
pcaout <- prcomp(t(topca))
summary(pcaout)
plot(pcaout, type = 'lines')
gg <- biplot(pcaout)
library(ggplot2)
ggplot(as.data.frame(pcaout$x), aes(x=PC1, y=PC2, label = row.names(pcaout$x))) + 
  geom_point(aes(col=factor(c('H1299', 'A431_replicate_a',
                              'A431_replicate_b', 'U251_MG_replicate_a', 'U251_MG_replicate_b', 'U2_OS_replicate_a',
                              'A549_replicate_a', 'A549_replicate_b', 'CACO2_replicate_a', 'CACO2_replicate_b',
                              'HEK293_replicate_a', 'HEK293_replicate_b', 'HeLa_replicate_a', 'HeLa_replicate_b',
                              'MCF7_replicate_a', 'MCF7-replicate_b', 'PC3_replicate_a')), size = 9)) + 
  scale_size_area(guide = F, max_size = 11) + 
  scale_colour_discrete(name = "Cell line") + 
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18)) + 
  guides(colour = guide_legend(override.aes = list(size=13)))+
  xlab('PC1 (15.3 explained)')+
  ylab('PC2 (12.5 explained)')+
  theme_bw()

png("myplot.png", width=4, height=4, units="in", res=300)

  

install.packages("ggbiplot")
library("ggbiplot")
install_github("ggbiplot", "vqv")
g <- ggbiplot(pcaout, pc.biplot = gg)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
              legend.position = 'top')
print(g)

install.packages('libxml-2.0.pc')
install_github("ggbiplot", "vqv")

install.packages("devtools")
setRepositories(libxml2-dev)

install.packages('libxml2-dev', dependencies=TRUE, repos='http://cran.rstudio.com/')

############################## The same with NumReads

beigerep1 = read.table('quantERR525589.sf', header = TRUE)
beigerep2 = read.table('quantERR525590.sf', header = TRUE)
brownrep1 = read.table('quantERR525591.sf', header = TRUE)
brownrep2 = read.table('quantERR525592.sf', header = TRUE)
whiterep1 = read.table('quantERR525593.sf', header = TRUE)
whiterep2 = read.table('quantERR525594.sf', header = TRUE)
frame = cbind(beigerep1$NumReads, beigerep2$NumReads, brownrep1$NumReads, brownrep2$NumReads, whiterep1$NumReads, whiterep2$NumReads)
log_frame = cbind(log(beigerep1$NumReads, 2), log(beigerep2$NumReads, 2), log(brownrep1$NumReads, 2), log(brownrep2$NumReads, 2), log(whiterep1$NumReads, 2), log(whiterep2$NumReads, 2))

is.na(log_frame) <- do.call(cbind,lapply(log_frame, is.infinite))
clear_frame <- na.omit(log_frame)

mtab <- data.frame(clear_frame)
topca = as.matrix(mtab)
pcaout = prcomp(t(topca))
library(ggplot2)
ggplot(as.data.frame(pcaout$x), aes(x=PC1, y=PC2, label = row.names(pcaout$x))) + 
  geom_point(aes(col=factor(c('bg', 'bg', 'brwn', 'brwn', 'white', 'white')), size = 9)) + 
  scale_size_area(guide = F, max_size = 15) + 
  scale_colour_discrete(name = "Cell type") + 
  theme(legend.text=element_text(size=15), legend.title=element_text(size=15)) + 
  guides(colour = guide_legend(override.aes = list(size=9)))

################# Now BAM-quant analysis

brownrep1 = read.table('Brown_1.sf', header = TRUE)
whiterep1 = read.table('White_1.sf', header = TRUE)
frame = cbind(brownrep1$NumReads, whiterep1$NumReads)
log_frame = cbind(log(brownrep1$NumReads, 2), log(whiterep1$NumReads, 2))

is.na(log_frame) <- do.call(cbind,lapply(log_frame, is.infinite))
clear_frame <- na.omit(log_frame)

mtab <- data.frame(clear_frame)
topca = as.matrix(mtab)
pcaout = prcomp(t(topca))
library(ggplot2)
ggplot(as.data.frame(pcaout$x), aes(x=PC1, y=PC2, label = row.names(pcaout$x))) + 
  geom_point(aes(col=factor(c('brwn', 'white')), size = 9)) + 
  scale_size_area(guide = F, max_size = 15) + 
  scale_colour_discrete(name = "Cell type") + 
  theme(legend.text=element_text(size=15), legend.title=element_text(size=15)) + 
  guides(colour = guide_legend(override.aes = list(size=9)))


################################# TPM > 100

beigerep1 = read.table('quantERR525589.sf', header = TRUE)
beigerep2 = read.table('quantERR525590.sf', header = TRUE)
brownrep1 = read.table('quantERR525591.sf', header = TRUE)
brownrep2 = read.table('quantERR525592.sf', header = TRUE)
whiterep1 = read.table('quantERR525593.sf', header = TRUE)
whiterep2 = read.table('quantERR525594.sf', header = TRUE)

TPM_1 <- beigerep1[order(-beigerep1$TPM),]
TPM_2 <- beigerep2[order(-beigerep2$TPM),]
TPM_3 <- brownrep1[order(-brownrep1$TPM),]
TPM_4 <- brownrep2[order(-brownrep2$TPM),]
TPM_5 <- whiterep1[order(-whiterep1$TPM),]
TPM_6 <- whiterep2[order(-whiterep2$TPM),]

TPM_1 <- TPM_1$TPM[1:100]
TPM_2 <- TPM_2$TPM[1:100]
TPM_3 <- TPM_3$TPM[1:100]
TPM_4 <- TPM_4$TPM[1:100]
TPM_5 <- TPM_5$TPM[1:100]
TPM_6 <- TPM_6$TPM[1:100]

log_frame = cbind(log(TPM_1, 2), log(TPM_2, 2), log(TPM_3, 2), log(TPM_4, 2), log(TPM_5, 2), log(TPM_6, 2))

is.na(log_frame) <- do.call(cbind,lapply(log_frame, is.infinite))
clear_frame <- na.omit(log_frame)

mtab <- data.frame(clear_frame)
topca = as.matrix(mtab)
pcaout = prcomp(t(topca))
library(ggplot2)
ggplot(as.data.frame(pcaout$x), aes(x=PC1, y=PC2, label = row.names(pcaout$x))) + 
  geom_point(aes(col=factor(c('bg', 'bg', 'brwn', 'brwn', 'white', 'white')), size = 9)) + 
  scale_size_area(guide = F, max_size = 15) + 
  scale_colour_discrete(name = "Cell type") + 
  theme(legend.text=element_text(size=15), legend.title=element_text(size=15)) + 
  guides(colour = guide_legend(override.aes = list(size=9)))