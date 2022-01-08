###This is a code to reproduce statistical analysis of proteomics comparison of HCA and HIT cells

##Openinig_the_data
#We reccomed to set the woring directory to make easy to reproduce the code
#setwd("your directory")
dat <- data.frame(read.csv("data.csv", header = TRUE))

##preparing the data
head(dat)
str(dat)
rownames(dat) <- dat[,3]
dat <- dat[,12:17]
head(dat)

fact <- data.frame(colnames(dat), c("HIT_2", "HIT_4","HIT_6","HCA_2","HCA_4","HCA_6"), c("HIT", "HIT","HIT","HCA","HCA","HCA"))
rownames(fact) <- fact[,2]
colnames(dat) <- fact[,2]
colnames(fact) <- c("old_names", "new_names", "type")
fact$type <- as.factor(fact$type)

head(fact)
head(dat)


##Qualittive analysis
##creating_the_venn_diagram
dat_HIT <- dat[which(rowMeans(!is.na(dat[,1:3])) >= 2/3), ]
dat_HCA <- dat[which(rowMeans(!is.na(dat[,4:6])) >= 2/3), ]


if(!require(gplots)){
  install.packages("gplots")
  library(gplots)
}
#list of group-specific names
v.table1 <- venn(list(rownames(dat_HIT), rownames(dat_HCA)))
print(v.table1)
str(v.table1)
#Then we extract group-specific accessions to excell manually



## Quantitative analysis
#Removing rows with a lot of missing values
dat_full <- dat[which(rowMeans(!is.na(dat)) >= 4/6), ]
mean(complete.cases(dat_full))
colSums(is.na(dat_full))

#variance stabilization normalization of protein expression data
#raw data
if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$type]
boxplot(dat_full, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)
#dev.off()

#imputation_of_missed_values_by_knn
if(!require(impute)){
  install.packages("impute")
  library(impute)
}
dat_knn1 <- impute.knn(t(dat_full), k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)

#log2 transformation
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)

#Quantile normalization
if(!require(limma)){
  install.packages("limma")
  library(limma)
}
dat_norm <- normalizeQuantiles(dat_log)
head(dat_norm)
#tiff('Norm_dat.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
boxplot(dat_norm, col = cols, main = "Normalized data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)
#dev.off()
mean(complete.cases(dat_norm))
head(dat_norm)
tail(dat_norm)

#MA-plot
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(dat_knn[, c(1:3)], dat_knn[, c(4:6)], main = "Raw data")
maplot(dat_log[, c(1:3)], dat_log[, c(4:6)], main = "Log-expression data")
maplot(dat_norm[, c(1:3)], dat_norm[, c(4:6)], main = "Normalized data")

##Cluster analysis
t_dat_norm <- t(dat_norm)
# nMDS
if(!require(vegan)){
  install.packages("vegan")
  library(vegan)
}
dat_ord <- metaMDS(t_dat_norm,
                   distance = "euclidean",
                   autotransform = FALSE)
dat_ord$stress

ord_ggplo <- data.frame(fact, scores(dat_ord, display = "sites"))
head(ord_ggplo, 2)

ord_ggplo_sp <- data.frame(scores(dat_ord, display = "species"))
ord_ggplo_sp$Species <- rownames(ord_ggplo_sp)
head(ord_ggplo_sp, 2)

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
#tiff('nMDS.tiff', units="in", width=12, height=8, res=600, compression = 'lzw')
ggplot() +
  geom_point(data = ord_ggplo, 
             aes(x = NMDS1, y = NMDS2, colour = type, 
                 shape = type), size = 4)
#dev.off()


#PCA
if(!require(mixOmics)){
  install.packages("mixOmics")
  library(mixOmics)
}


dat_pca <- pca(t(dat_norm), ncomp = 4, center = TRUE)

#tiff('PCA_cell_type_2D.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 3), ind.names = F, 
          group = fact$type, legend = TRUE, ellipse = TRUE,
          title = 'PCA')
#dev.off()

# sPLS-DA
ordination.optimum.splsda <- splsda(t_dat_norm, fact$type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PLSDA.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = T, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)


##Differential expression analysis by limma
fact$type
#HCA is a basic level!
X <- model.matrix(~ fact$type)
X

fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit <- eBayes(fit)

# Dif_expr_table
topTable(efit, coef = 2)
full_list <- topTable(efit, coef = 2, number = length(dat_knn))
#write.csv(full_list,'Dif_expr_HCA_vs_HIT.csv')
head(full_list)


#Vulcano
library(EnhancedVolcano)

head(full_list)
rownames(full_list) <- substr(rownames(full_list),1,6)
head(full_list)

#tiff('Vulcano_HCA_vs_HIT.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(full_list,
                lab = rownames(full_list),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 2,
                xlim = c(-10, 10),
                ylim = c(0, 5),
                title ="HCA_vs_HIT",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()



