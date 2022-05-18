
# Read read counts created by featureCounts
counts = read.delim("report/featureCounts/counts_against_genome.txt", skip=1, header=TRUE, row.names=1)
names(counts)

# Slice only read counts
counts = counts[-c(1:5)]

# reorder columns
counts = counts[,c(5,6,7,8,1,2,3,4)]


# import edger library
library(edgeR)

# Grouping by control and metal-exposed
group = factor(c(rep("control", 4), rep("metal", 4)))



# Convert read count table to dge object
dge = DGEList(counts=counts, group=group)
str(dge)
levels(dge$samples$group)

# Filter out low expressed records
keep = filterByExpr(dge)
str(keep)
class(keep)
dge = dge[keep, , keep.lib.sizes=FALSE]
str(dge)

# Normalize library sizes
dge = calcNormFactors(dge)
str(dge)

# Estimate dispresion
dge = estimateDisp(dge)
str(dge)

# Testing for DE genes using Fisher's exact test
et = exactTest(dge)
str(et)
topTags(et)
summary(decideTests(et))
jpeg("report/edger/DE_highlighted.jpg")
plotMD(et)
dev.off()

# Draw heatmap 
scale = dge$samples$lib.size * dge$samples$norm.factors
normed = round(t(t(counts)/scale) * mean(scale))
data = topTags(et, n=nrow(dge$counts))$table
total = merge(data, normed, by='row.names')
col = c('control_gill', 'control_intestine', 'control_kidney', 'control_liver','metal_gill', 'metal_intestine', 'metal_kidney', 'metal_liver')
count_data = subset(total, total$FDR <= 0.05)
values = count_data[,6:ncol(count_data)] #data.frame
colnames(values) = col
values = as.matrix(values) # matrix
values = jitter(values, factor=1, amount=0.0001)
zscores = NULL
for (i in 1:nrow(values)) {
    row = values[i,]
    zrow = (row - mean(row)) / sd(row)
    zscores = rbind(zscores, zrow)
}
row_names = count_data[,1]
row.names(zscores) = row_names
library(gplots)
jpeg("report/edger/heatmap.jpg")
heatmap.2(zscores, col=greenred, density.info='none', Colv=NULL, dendrogram="row", trace="none", margins=c(12,18), lhei=c(1,5))
dev.off()

head(total)
total$baseMean = rowMeans(total[,6:13])
total$baseMeanControl = rowMeans(total[,6:9])
total$baseMeanMetal= rowMeans(total[,10:13])
total$foldChange = 2 ^ total[,2]
names(total)[1] = 'name'
names(total)[6:13] = col
total = total[, c(1,14,15,16,17,2,3,4,5,6,7,8,9,10,11,12,13)]
total = total[with(total, order(FDR, -foldChange)),]
write.csv(total, file='report/edger/DE_MetalvsControl.csv', row.names=FALSE, quote=FALSE)
# ==== Gill ====
# Read read counts created by featureCounts
counts = read.delim("report/featureCounts/counts_against_genome.txt", skip=1, header=TRUE, row.names=1)
names(counts)

# Slice only read counts
counts = counts[c(6,10)]
counts[1:10,]

# import edger library
library(edgeR)

Grouping by control and metal-exposed
group = c("metal", "control")
group


# Convert read count table to dge object
dge = DGEList(counts=counts, group=group)
str(dge)

# Filter out low expressed records
keep = filterByExpr(dge)
str(keep)
class(keep)
dge = dge[keep, , keep.lib.sizes=FALSE]
str(dge)

# Normalize library sizes
dge = calcNormFactors(dge)
str(dge)

# Estimate dispresion
dge = estimateDisp(dge)
str(dge)

# Testing for DE genes using Fisher's exact test
et = exactTest(dge)
str(et)
topTags(et)

