
read_data <- function(){
  numb <- seq(1:216139)
  all_df <- as.data.frame(numb)
  for (i in dir()){
    tmp_df <- read.table(i, header = FALSE)
    colnames(tmp_df) <- c(i, i, i, i)
    all_df <- cbind(all_df, tmp_df[4])
  }
  print(paste("Merged", length(dir()), "files"))
  all_df$mean <- apply(all_df[,2:ncol(all_df)], 1, mean)
  
  return(all_df$mean)
}


getwd()
# run in folder where all folders with raw counts are
big_read_data <- function(){
  numb <- seq(1:216139)
  all_df <- as.data.frame(numb)
  all_folders <- list.files(getwd())
  counter = seq(length(all_folders))
  for (j in counter){
    local_df <- as.data.frame(numb)
    setwd(all_folders[j])
    setwd('./R')
    for (i in dir()){
      tmp_df <- read.table(i, header = F)
      local_df <- cbind(local_df, tmp_df[4])
    }
    name_t <- apply(local_df[,2:ncol(local_df)], 1, sum)
    all_df <- cbind(all_df, name_t)
    
    setwd('../../')
  }
  colnames(all_df) <- c('num', all_folders)
  print(paste("Merged", length(all_folders), "files"))
  return(all_df)
}

test_numreads <- big_read_data()

View(b)

library_sizes <- read.csv('norm_by_tissues.csv', header = TRUE)

mean_by_tissues <- apply(library_sizes, 2, mean, na.rm = T)
mean_by_tissues <- as.data.frame(test)
nod <- mean(test[,1])
coeff_for_norm <- mean_by_tissues / nod
coeff_for_norm <- t(coeff_for_norm)
coeff_for_norm <- as.vector(coeff_for_norm)

all_numreads_test1 <- test_numreads[,2:33]
sum(all_numreads_test1)

##### new normalisation

all_numreads_test1 <- test_numreads[,2:33]
coverage_by_tissues <- apply(all_numreads_test1, 2, sum, na.rm = T)
mean_coverage <- mean(coverage_by_tissues)
coverage_by_tissues <- coverage_by_tissues / mean_coverage
View(coverage_by_tissues)

# normed_data <- sweep(all_numreads_test, MARGIN=2, coverage_by_tissues,`*`)

mean_depth <- mean_coverage
normed_data <- as.data.frame(apply(all_numreads_test1, 2, function(x) x*(mean_depth/sum(x))))

### generating the table with raw counts and names
TCT <- read.table('adrenal_4a_tech1.sf', header = FALSE)
TCT <- TCT[1]
colnames(TCT) <- 'Name of transcript'

test_numreads_with_TCT <- cbind(TCT, test_numreads)

write.table(all_numreads_with_TCT, "Salmon_Numreads_PT_test.txt", sep="\t") 


######## возвращаемся к проекту

length_of_TCT <- read.table('adrenal_4a_tech1.sf', header = FALSE)
length_of_TCT <- length_of_TCT[2]

length_of_TCT_by_length_of_read <- 101 / length_of_TCT

length_of_TCT_by_length_of_read <- as.vector(length_of_TCT_by_length_of_read[,1])
str(length_of_TCT_by_length_of_read)

##### Кол-во нормализованных ридов умноженное на отношение длина рида/длина транскрипта

ak_data <- sweep(normed_data, MARGIN=1, length_of_TCT_by_length_of_read,`*`)

ak_data_ok <- subset(ak_data > 5)

ak_data_ok <- ak_data_ok + 0


salmon_binary <- ak_data_ok
covered = apply(salmon_binary, 1, function(x) sum(x) > 0)
number_of_covered = sum(covered)
rownames(salmon_binary) = length_of_TCT[,1]

contributors = c()
contributions = c()
for (i in 1:31){
  sums = apply(salmon_binary, 2, sum)
  the_best_tissue = colnames(salmon_binary)[which.max(sums)]
  print(the_best_tissue)
  contributors = c(contributors, the_best_tissue)
  contributions = c(contributions, max(sums))
  transcripts_to_keep = salmon_binary[, the_best_tissue] == 0
  salmon_binary = salmon_binary[transcripts_to_keep, ]
  salmon_binary = salmon_binary[, -which.max(sums)]
}

contributions = c(contributions, 1)
View(contributors)
contributors <- c(contributors, 'appendix')
output_s <- cbind(contributions, contributors)
View(output_s)

colnames(output_s) <- c('Num of Transcripts', 'Tissue')

write.table(output_s, "Basis_Salmon.txt", sep="\t") 

basis <- read.table('Salmon_Numreads_PT_test.txt', header = TRUE)

######### Pie chart

ggplot(df[1:10,], aes(x = factor(1), y = value[1:10], fill = tissue))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15))+
  xlab('T H I S  I S')+
  ylab('S P A R T A ! ')+
  scale_colour_discrete(name = "Tissue")+
  scale_fill_discrete(breaks=c("testis - 59996", "spleen - 13352", "brain - 7308", "bonemarrow - 5161", "duodenum - 4170", "kidney - 2801", "prostate - 2119", "skin - 1750", "appendix - 1413", "fallopiantube - 1240"))

df <- data.frame(variable, value)
value <- c(59996, 13352, 7308, 5161, 4170, 2801, 2119, 1750, 1413, 1240)
variable <- c("testis", "spleen", "brain", "bonemarrow", "duodenum", "kidney", "prostate", "skin", "appendix", "fallopiantube")
variable

tissue <- variable[1:10]

tissue <- tissue_numbers

tissue_numbers <- c("testis - 59996", "spleen - 13352", "brain - 7308", "bonemarrow - 5161", "duodenum - 4170", "kidney - 2801", "prostate - 2119", "skin - 1750", "appendix - 1413", "fallopiantube - 1240")

########## Vein diagram

install.packages(VennDiagram)

rm(list=U251_MG_replicate_b,envir=sys.frame(-1))

rm_obj <- function(obj){
  a <-deparse(substitute(obj))
  print (a)
  print(ls(envir=sys.frame(-1)))  
  rm(list=a,envir=sys.frame(-1))
}

rm_obj(A431_replicate_a)

library(VennDiagram)

RSEM <- read.table('testis_7a.isoforms.results', header = TRUE)
Salmon <- read.table('testis_7a.sf', header = TRUE)
Kallisto <- read.table('testis_7a.tsv', header = TRUE)
str(RSEM)

RSEM_sorted <- RSEM[with(RSEM, order(RSEM$expected_count, RSEM$transcript_id)), ]
Salmon_sorted <- Salmon[with(Salmon, order(Salmon$NumReads, Salmon$Name)), ]
Kallisto_sorted <- Kallisto[with(Kallisto, order(Kallisto$est_counts, Kallisto$target_id)), ]

coverage_RSEM <- sum(RSEM_sorted$expected_count)
coverage_Kallisto <- sum(Kallisto_sorted$est_counts)
coverage_Salmon <- sum(Salmon_sorted$NumReads)

coverage <- c(coverage_RSEM, coverage_Kallisto, coverage_Salmon)

rr <- calculate.overlap(coverage)

########## Ggplot_density
library(ggplot2)
colnames(Salmon_rows$V2) <- Salmon

Salmon_rows <- Salmon$NumReads
Salmon_rows <- cbind(Salmon_rows, 'Salmon')
Salmon_rows <- as.data.frame(Salmon_rows)
Salmon_rows$V2 <- as.factor(Salmon_rows$V2)
Salmon_rows$Salmon_rows <- as.numeric(Salmon_rows$Salmon_rows)
colnames(Salmon_rows) <- c('Num_of_reads', 'Tool')
str(Salmon_rows)

RSEM_rows <- RSEM$expected_count
RSEM_rows <- cbind(RSEM_rows, 'RSEM')
RSEM_rows <- as.data.frame(RSEM_rows)
RSEM_rows$V2 <- as.factor(RSEM_rows$V2)
RSEM_rows$RSEM_rows <- as.numeric(RSEM_rows$RSEM_rows)
colnames(RSEM_rows) <- c('Num_of_reads', 'Tool')
str(RSEM_rows)

Kallisto_rows <- Kallisto$est_counts
Kallisto_rows <- cbind(Kallisto_rows, 'Kallisto')
Kallisto_rows <- as.data.frame(Kallisto_rows)
Kallisto_rows$V2 <- as.factor(Kallisto_rows$V2)
Kallisto_rows$Kallisto_rows <- as.numeric(Kallisto_rows$Kallisto_rows)
colnames(Kallisto_rows) <- c('Num_of_reads', 'Tool')
str(Kallisto_rows)

final_rows <- rbind(Salmon_rows, RSEM_rows, Kallisto_rows)
final_rows <- subset(final_rows, final_rows$Num_of_reads>100)

d <- ggplot(final_rows, aes(Num_of_reads, fill = Tool)) +
  geom_density(alpha = 0.6)
d + scale_x_log10()

############# Boxplot Sum of Reads

sum(RSEM$expected_count)
sum(Salmon$NumReads)
sum(Kallisto$est_counts)

summm <- ggplot(DATAFR_kal, aes(Tool, Total_sum))+
  geom_bar(stat = 'identity', aes(fill = Tool))+ 
  guides(colour = guide_legend(override.aes = list(size=14)))+
  xlab('Tool')+
  ylab('Sum')+
  theme_bw()+
  ggtitle('Total number of good-covered transcripts (Testis 7a)')+
  scale_fill_manual(values=c('yellow','brown', 'blue'))+
  theme(axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        panel.background = element_blank(),
        plot.title = element_text(size = 20, colour = "black"))

summm

DATAFR_kal <- data.frame(variable1=rep(c("RSEM", "Salmon", "kallisto")), value1=c(30079,31917,29229))
# DATAFR_kal$Total_sum <- log(DATAFR_kal$Total_sum)

colnames(DATAFR_kal) <- c('Tool', 'Total_sum')

######## How much it would be in RSEM?

xdata <- cbind(RSEM$expected_count, Salmon$NumReads, Kallisto$est_counts)
colnames(xdata) <- c('RSEM', 'Salmon', 'Kallisto')

x_data <- as.data.frame(apply(xdata, 2, function(x) x*length_of_TCT_by_length_of_read))

x_data_ok <- subset(x_data > 10)

x_data_ok <- x_data_ok + 0

x_data_ok <- as.data.frame(x_data_ok)
sum(x_data_ok$RSEM)
sum(x_data_ok$Salmon)
sum(x_data_ok$Kallisto)


### TRASH

most_by_data_ok <- apply(ak_data_ok, 2, sum, na.rm = T)
View(most_by_data_ok)
most_by_data_ok <- as.data.frame(most_by_data_ok)
hist(most_by_data_ok)

big_read_data <- function(dataframe){
  counting_1 <- apply(ak_data_ok, 2, sum, na.rm = T)
  counting_1 <- as.data.frame(counting_1)
  sorted_counting_1 <- counting_1[with(counting_1, order(counting_1)), ]
  sorted_counting_1 <- as.data.frame(sorted_counting_1)
  for (j in counter){
    local_df <- as.data.frame(numb)
    setwd(all_folders[j])
    setwd('./R')
    for (i in dir()){
      tmp_df <- read.table(i, header = F)
      local_df <- cbind(local_df, tmp_df[4])
    }
    name_t <- apply(local_df[,2:ncol(local_df)], 1, mean)
    all_df <- cbind(all_df, name_t)
    
    setwd('../../')
  }
  colnames(all_df) <- c('num', all_folders)
  print(paste("Merged", length(all_folders), "files"))
  return(all_df)
}




