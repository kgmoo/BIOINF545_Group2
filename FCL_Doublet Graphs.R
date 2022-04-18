#Load data (skip this if data is already in global environment)
# epithelial counts
epi.counts = readMM('gene_sorted-Epi.matrix.mtx')
rownames(epi.counts) = readLines('Epi.genes.tsv')
colnames(epi.counts) = readLines('Epi.barcodes2.tsv')

# stromal counts
fib.counts = readMM('gene_sorted-Fib.matrix.mtx')
rownames(fib.counts) = readLines('Fib.genes.tsv')
colnames(fib.counts) = readLines('Fib.barcodes2.tsv')

# immune counts
imm.counts = readMM('gene_sorted-Imm.matrix.mtx')
rownames(imm.counts) = readLines('Imm.genes.tsv')
colnames(imm.counts) = readLines('Imm.barcodes2.tsv')

# Load metadata for discovery and validation cohorts
meta = read.table('all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
##########################################


#we can look at library complexity (the number of genes detected per cell) to identify potential shortcomings in the data before doing any manipulation
#first, make objects so we can pull the counts to a graph
#use Matrix::--- to force function from Matrix library (ensures correct handling of sparse matrix)
epi_counts_per_cell <- Matrix::colSums(epi.counts)
epi_counts_per_gene <- Matrix::rowSums(epi.counts)
epi_genes_per_cell <- Matrix::colSums(epi.counts>0)
hist(log10(epi_counts_per_cell+1))
hist(log10(epi_genes_per_cell+1))
plot(epi_counts_per_cell, epi_genes_per_cell, log='xy')
#low end outliers indicate failed libraries while high end outliers indicate doublets

plot(sort(epi_genes_per_cell), xlab='Cell', log='y', main='Genes per Cell (ordered)', ylab='Number of Genes')
imm_counts_per_cell <- Matrix::colSums(imm.counts)
imm_counts_per_gene <- Matrix::rowSums(imm.counts)
imm_genes_per_cell <- Matrix::colSums(imm.counts>0)
hist(log10(imm_counts_per_cell+1))
hist(log10(imm_genes_per_cell_imm+1))
plot(imm_counts_per_cell, imm_genes_per_cell, log='xy')
plot(sort(imm_genes_per_cell), xlab = 'cell', log='y', main='Genes per Cell (ordered)', ylab='Number of Genes')

fib_counts_per_cell <- Matrix::colSums(fib.counts)
fib_counts_per_gene <- Matrix::rowSums(fib.counts)
fib_genes_per_cell <- Matrix::colSums(fib.counts>0)
hist(log10(fib_counts_per_cell+1))
hist(log10(fib_genes_per_cell+1))
plot(fib_counts_per_cell, fib_genes_per_cell, log='xy')
plot(sort(fib_genes_per_cell), xlab = 'cell', log='y', main='Genes per Cell (ordered)', ylab='Number of Genes')