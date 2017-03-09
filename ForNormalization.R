
library(edgeR)

data <- read.table('CPM_data_filtered.txt', sep = '\t', header = FALSE, row.names = 1)

group <- factor(c('NoDox', '72hA', '72hB'))

d <- DGEList(counts=data, group=group)

data_norm = round( cpm(data, normalized.lib.size = T), 4 )

write.table(data_norm, 'CPM_data_filtered_libsizeNormed.txt', sep = '\t', row.names = T, col.names = F, quote = F)
