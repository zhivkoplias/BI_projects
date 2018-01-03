# Optimal basis
# Reading data
dire = readline('Where are TSVs?   ')
setname = readline('What are your samples?   ')
kallisto_data = data.frame(c(1:216139))
for (exp in dir(dire)){
   tab = read.table(paste0(dire, exp), header=T)
   tpms = tab[,4]
   kallisto_data = cbind(kallisto_data, tpms)
 }
kallisto_data = kallisto_data[, c(2:ncol(kallisto_data))]
# samples = factor(sapply(strsplit(dir(dire), '_'), function(elt) elt[1]))
# colnames(kallisto_data) = c(as.character(samples))
# kallisto_data$LENGTH = tab[,2]
#Output here is the data frame of all the bullshit

#Preparing binary matrix (covered/not)
kallisto_data <- all_numreads_with_TCT


kallisto_summed = data.frame(kallisto_data$`Name of transcript`)
for (i in 1:length(levels(samples))){
  this_tissue_cols = which(samples == levels(samples)[i])
  this_tissue_sum = apply(kallisto_data[, this_tissue_cols], 1, sum)
  kallisto_summed = cbind(kallisto_summed, this_tissue_sum)
}
colnames(kallisto_summed) = c('LENGTH', levels(samples))
mean_depth = mean(apply(kallisto_summed[, 2:ncol(kallisto_summed)], 2, sum))
kallisto_norm = as.data.frame(apply(kallisto_summed[, 2:ncol(kallisto_summed)], 2, function(x) x*(mean_depth/sum(x))))

# Threshold in the next line specifies the point of coverage cutoff
kallisto_binary = as.data.frame(apply(kallisto_norm, 2, function(elt) ifelse((elt*101)/kallisto_summed$LENGTH >= 100, 1, 0)))

salmon_binary <- ak_data_ok
covered = apply(salmon_binary, 1, function(x) sum(x) > 0)
number_of_covered = sum(covered)
rownames(salmon_binary) = length_of_TCT[,1]

# Finding algo
salmon_binary_filtered = salmon_binary[covered, ]
contributors = c()
contributions = c()
for (i in 1:31){
  sums = apply(salmon_binary_filtered, 2, sum)
  the_best_tissue = colnames(salmon_binary_filtered)[which.max(sums)]
  print(the_best_tissue)
  contributors = c(contributors, the_best_tissue)
  contributions = c(contributions, max(sums))
  transcripts_to_keep = salmon_binary_filtered[, the_best_tissue] == 0
  salmon_binary_filtered = salmon_binary_filtered[transcripts_to_keep, ]
  salmon_binary_filtered = salmon_binary_filtered[, -which.max(sums)]
}

contributions = c(contributions, sum(salmon_binary_filtered))
View(contributors)
output_s <- cbind(contributions, contributors)
output_s <- as.data.frame(output_s)
colnames(output_s) <- c('Num of Transcripts', 'Tissue')

contributions <- c(contributions, 1)
contributors <- c(contributors, 'appendix')

plot(t(output_s))

write.table(output_s, "Basis_Salmon.txt", sep="\t") 



contributions
contributors = c(contributors, levels(ak_data_ok)[-which(levels(ak_data_ok) %in% contributors)])
names(contributions) = contributors
basis = as.data.frame(contributions)
colnames(basis) = "number"

#Using protein-coding ones
#This part needs grepped protein-coding transcripts names
prc = read.table('/media/DISK1/reference/Gencode/human_23/protein_coding_ids')
covered_prc = as.character(prc[prc$V1 %in% rownames(kallisto_binary)[covered], ])
kallisto_binary_prc = kallisto_binary[covered_prc, ]
contributors_prc = c()
contributions_prc = c()
for (i in 1:31){
  sums_prc = apply(kallisto_binary_prc, 2, sum)
  the_best_tissue_prc = colnames(kallisto_binary_prc)[which.max(sums_prc)]
  print(the_best_tissue_prc)
  contributors_prc = c(contributors_prc, the_best_tissue_prc)
  contributions_prc = c(contributions_prc, max(sums_prc))
  transcripts_to_keep_prc = kallisto_binary_prc[, the_best_tissue_prc] == 0
  kallisto_binary_prc = kallisto_binary_prc[transcripts_to_keep_prc, ]
  kallisto_binary_prc = kallisto_binary_prc[, -which.max(sums_prc)]
}
contributions_prc = c(contributions_prc, sum(kallisto_binary_prc))
contributors_prc = c(contributors_prc, levels(samples)[-which(levels(samples) %in% contributors_prc)])
names(contributions_prc) = contributors_prc
basis_prc = as.data.frame(contributions_prc)
colnames(basis_prc) = "number"

#Totalling
total_basis = data.frame(all_names = rownames(basis), all_counts = basis$number)

#Plotting the final data as the cumgraph
total_basis$all_cumulate = cumsum(total_basis$all_counts)
total_basis$prc_cumulate = cumsum(total_basis$prc_counts)
trf_basis = as.data.frame(rbind(c(1, 1), total_basis[, 5:6]))
colnames(trf_basis) = c('all', 'prc')
library(ggplot2)
library(reshape2)
toplot = melt(trf_basis)
# From here needs running manually, for some bullshit happens
toplot$coord = rep(1:33, 2)
ggplot(toplot, aes(x = coord, y = value, group = variable, col = variable)) + geom_line(size = 1.5) + theme_bw() + 
  ggtitle('Cumulative number of represented transcripts') + xlab('tissues') + scale_y_continuous(limits = c(1, 120000)) + 
  scale_color_manual(values = c('red', 'green'), guide = F) + geom_abline(aes(intercept = number_of_covered, slope = 0), linetype = 2, size = 1, col = 'red') +
  geom_abline(aes(intercept = length(covered_prc), slope = 0), linetype = 2, size = 1, col = 'green')