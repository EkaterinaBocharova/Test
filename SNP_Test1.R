setwd('/Users/katya/SNP')
library(vcfR)

vcf <- read.vcfR(file = 'test_data1.vcf')

library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")

# Трансформация в genlight (подходит для SNPs)
gl <- vcfR2genlight(vcf)
gl

# Расчет эвклидовых генетических дистанций
distgenEUCL <- dist(gl, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
distgenEUCL

# Визуализация Pheatmap
library("pheatmap")
sampleDistMatrix = as.matrix(distgenEUCL)
pheatmap(sampleDistMatrix, clustering_distance_rows=distgenEUCL, clustering_distance_cols=distgenEUCL, fontsize_row = 5, fontsize_col = 5)

# Визуализация NJ
tre <- nj(dist(as.matrix(gl)))
plot(tre, cex=0.4)
title("NJ tree")

# Визуализация PCA
pca1 <- glPca(gl)
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=1)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topleft", inset=.05, ratio=.3)
title("PCA")



