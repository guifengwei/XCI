
library(edgeR)

data <- read.table('GenesCounts.txt', header = FALSE, sep = "\t", row.names = 1)

data_cpms = round(cpm(data), 4)

write.table(data_cpms, 'CPM_data.txt',  sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
