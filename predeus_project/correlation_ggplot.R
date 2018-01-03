RSEM <- read.table('heart_6a_tech1.isoforms.results', header = TRUE)
Salmon <- read.table('brain_3c.sf', header = TRUE)
Kallisto <- read.table('abundance.tsv', header = TRUE)
str(RSEM)

RSEM_sorted <- RSEM[with(RSEM, order(RSEM$transcript_id, RSEM$expected_count)), ]
Salmon_sorted <- Salmon[with(Salmon, order(Salmon$Name, Salmon$NumReads)), ]

coverage_RSEM <- sum(RSEM$effective_length)
sum(Salmon$NumReads)
coverage_Kallisto <- sum(Kallisto$est_counts)
coverage_Kallisto

cor(RSEM_sorted$expected_count, Salmon_sorted$NumReads)

test1 <- cbind(RSEM_sorted$expected_count, Salmon_sorted$NumReads)
colnames(test1) <- c('RSEM', 'Salmon')
test1 <- as.data.frame(test1)

cor(test1$RSEM, test1$Salmon)

ggplot(test1, aes(RSEM, Salmon))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

fit2 <- lm(RSEM ~ Salmon, test1)
library(car)
influencePlot(fit2)

fit2_new <- subset(test1, !(row.names(test1) %in%
                               c('5239', '23357')))

ggplot(fit2_new, aes(RSEM, Salmon))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()