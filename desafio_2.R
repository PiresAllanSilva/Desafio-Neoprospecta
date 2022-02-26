#install.packages("tidyverse", dependencies = TRUE)
library(tidyverse)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
#install.packages("ggfortify", dependencies = TRUE)
library(ggfortify)

#setwd("Path/to/work/directory")

otu_table <- read.csv("tables/otu_table_tax_amostras.tsv", sep = "\t")
metadata <- read.csv("metadata.csv", sep = "\t") # A tabela foi modificada manualmente para contar os dados na mesma ordem da OTU
cols <- (colnames(otu_table[2:20]))
df <- data.frame(DWP=character(), Abundance=integer(), ID=character())
names(df)<-c("DWP","Abundance", "ID")
df2 <- data.frame(1:50)
c <- 1

for(i in cols){
  temp <- otu_table %>%
    arrange(desc(otu_table[i]))
  mabu <- dplyr::slice(temp,1:50)
  new_columns <- data.frame(cbind(mabu[1],mabu[i]))
  if (c == 1){
    df2 <- data.frame(cbind(mabu[1],mabu[i]))
    c <- 2
  }
  if (c == 2){
    df2 <- full_join(df2, new_columns, by = "OTU")
  }
  count <- colSums(unname(mabu[i]))
  days <- (str_split(i, "D", n=2))
  days <- as.array(days[[1]])  
  days <- str_split(days[2], "_")
  days <- days[[1]]  
  de <- data.frame(days[1], count, i)
  names(de)<-c("DWP","Abundance", "ID")
  df <- rbind(df, de)
}

df2[is.na(df2)] <- 0
# Tabela com as união dos dados de abundância das 50 espécies mais abundantes por dia
write.csv(df2,"tables/", row.names = FALSE)

#Análises diferenciais
deseq <- otu_table[,-1]
rownames(deseq) <- otu_table[,1]
coldata <- sample[, c(1, 311)]
dds <- DESeqDataSetFromMatrix(countData = deseq, colData = coldata, design = ~time)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="time_Late_vs_Early")
res <- lfcShrink(dds, coef="time_Late_vs_Early")
resNorm <- lfcShrink(dds, coef=2, type="normal")
res05 <- results(dds, alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE) #Number of bacterial species with a padj <0.05
# Exporta a tabela com os dados gerados pelo DESeq2
write.csv(res05,"tables/", row.names = FALSE)


# Gráfico de PCoA das amostras. O arquivo otu.csv foi preparado por meio do 
# script transpose.py para facilitar a análise. Foi utilizada a biblioteca
# pandas no script

sample <- read.csv("otu.csv", sep = ",", header = FALSE)
pca_sample <- prcomp(sample[2:310], center = TRUE)
sample <- cbind(sample, metadata[3])

pdf("PCoA.pdf")
PCA<-data.frame(pca_sample$x,time=sample$time)
ggplot(PCA,aes(x=PC1,y=PC2,col=time, label=sample$V1))+
  geom_point(size=4,alpha=0.8)+ 
  scale_color_manual(values = c("#000090","#000000"))+
  geom_text(vjust = 2, nudge_y = 100)+
  theme_classic()
dev.off
pdf("PCoA_nonames.pdf")
PCA<-data.frame(pca_sample$x,time=sample$time)
ggplot(PCA,aes(x=PC1,y=PC2,col=time))+
  geom_point(size=4,alpha=0.8)+ 
  scale_color_manual(values = c("#000090","#000000"))+
  theme_classic()
dev.off


# Gráfico das contagens absolutas das 50 bactérias mais abundantes
# por dia depois do desmame
df$DWP = as.numeric(df$DWP)
df <- df[order(df$DWP),]
df$DWP = as.character(df$DWP)
pdf("Total_abundance.pdf")
ggplot(df, aes(y = Abundance, x = fct_inorder(DWP), fill = fct_inorder(DWP))) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 10) +
  xlab("Days after weaning") + 
  ylab("Abundance")+
  guides(fill=guide_legend(title="DWP"))
dev.off()


