setwd('/Users/katya/SNP')
library(vcfR)

vcf <- read.vcfR(file = 'test_data1.vcf')

library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")

gl <- vcfR2genlight(vcf)
gl

gl <- as.snpclone(gl)
gl

# Расчет эвклидовых генетических дистанций
distgenEUCL <- dist(gl, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
distgenEUCL
min(distgenEUCL)
max(distgenEUCL)

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



