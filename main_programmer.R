library('affy')
library("affyPLM")
library("AnnotationDbi")
library("hgu133plus2.db")
library("sva")
library("ggplot2")

#3.
cel_filenames <- list.files(pattern='\\.gz$') 
#all FILE NAMES that end with .gz will get stored into variable cel_filenames as a list

Data <- ReadAffy(filenames=cel_filenames)
# this will read all CEL files and store the results in Data variable
normalized_data <- rma(Data)
#this will  normalize all of the CEL files together and store it in variable normalized_data

#4.
PLMSet <- fitPLM(Data,  normalize=TRUE, background=TRUE)

RLE_stats <- RLE(PLMSet,type="stats")
RLE_medians = RLE_stats["median", ]

NUSE_stats <- NUSE(PLMSet, type = "stats")
NUSE_medians <- NUSE_stats["median", ]

hist(RLE_medians,col="lightgreen",nclass=8,main="Histogram of RLE medians",xlab="median RLE")
hist(NUSE_medians,col="lightgreen",nclass=8,main="Histogram of NUSE medians",xlab="median NUSE")

exprs_data = exprs(normalized_data)
write.csv(exprs_data, "rma_norm_exprs_data.csv")

#5.
metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
mod <- model.matrix(~normalizationcombatmod, data = metadata)
combat_exprs_data <- ComBat(dat = exprs_data, batch = metadata$normalizationcombatbatch, mod = mod)
write.csv(combat_exprs_data, 'combat_exprs_data.csv')

trans_combat_exprs_data <- t(combat_exprs_data)
scaled_trans_combat_exprs_data <- scale(trans_combat_exprs_data, center = TRUE, scale = TRUE)
scaled_combat_exprs_data <- t(scaled_trans_combat_exprs_data)

#6.
pca <- prcomp(scaled_combat_exprs_data, center = FALSE, scale = FALSE)
summary(pca)
pca$rotation
summary(pca)$importance[,1:2]

#7.
pc1_label  <- paste0("PC1 ",summary(pca)$importance[2,1]*100, "%")
pc2_label  <- paste0("PC2 ",summary(pca)$importance[2,2]*100, "%")

pca_rotation_df <- as.data.frame(pca$rotation)
ggplot(data = pca_rotation_df, mapping = aes(x = PC1, y = PC2)) +
  geom_point()+
  theme_bw() +
  labs(title = 'PCA plot', x= pc1_label, y=pc2_label)