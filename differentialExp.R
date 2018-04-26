source("https://bioconductor.org/biocLite.R")
biocLite("affy") ## Pre-processing Affy chips
biocLite("limma") ## Differential expression
library(affy)
library(limma)

limmaUsersGuide()

## Read in .CEL files

df <- ReadAffy()

## Diagnostic plots

pdf("boxplot_arrays_beforeNRM.pdf")
boxplot(df) # Boxplot of intensity values before normalization
dev.off()

png("spatial_intensity.png")
image(df) # 2D spatial colour images of log intensities
dev.off()

png("mva_pairs_beforeNRM.png")
mva.pairs(exprs(df)) # MVA plots
dev.off()

png("exp_intensitiy_hist_beforeNRM.png")
hist(df, type = "l", lty = 1, lwd = 3) # Density plots of log intensities
dev.off()

## Normalize data - rma = model based background correction, quantile normalisation,
## median polish expression measures

eset <- rma(df)

## Plots after normalisation

png("mva_pairs_afterNRM.png")
mva.pairs(exprs(eset)) 
dev.off()

pdf("boxplot_arrays_afterNRM.pdf")
boxplot(exprs(eset)) 
dev.off()

## differential expression
# Create a design set - 2 pairs of replicate array so we should esimate 2 parameters in the linear
# model

design <- model.matrix(~0+factor(c(1, 1, 2, 2))) # describe model to be fit

colnames(design) <- c("Pop1", "Pop2")

fit <- lmFit(eset, design)  # fit each probeset to model

cont.matrix <- makeContrasts(Pop1vsPop2 = Pop1-Pop2, levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)

efit <- eBayes(fit2)        # empirical Bayes adjustment (will throw warning message - you can ignore)
topTable(efit)      # table of differentially expressed probesets


## Plot the distribution of p-values over the experiemnt. If all p-values have nearly the same
## frequency, there is no change in expression in the experiment. A high frequency of p-values
## at zero is indicating a number of differentially expressed genes

pdf("histo_pval.pdf")
hist(efit$p.value)
dev.off()

## Make a heatmap 

# Take top 100 genes associated with phenotype

require(Heatplus)

geneID <- rownames(topTable(efit, num=100, adjust="bonferroni", p.value=0.05))

eset2 <- eset[geneID, ]                 

heat <- regHeatmap(exprs(eset2))
pdf("heat_map_top100.pdf")
plot(heat)
dev.off()
