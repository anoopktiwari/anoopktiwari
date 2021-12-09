
# link for fpkm calculation
# https://github.com/AAlhendi1707/countToFPKM/blob/master/R/countToFPKM-internal.R
# package counttoFPKM
# https://cran.r-project.org/web/packages/countToFPKM/countToFPKM.pdf
# https://rdrr.io/bioc/DESeq2/man/fpkm.html

(# Copy first 1000 lines to a separate test file
awk 'NR<1001' counts-mod-ak.txt > counts-test.txt) this is for shell 

# Copy the counts-test.txt to your computer using Dropbox or other means and read it into the R
read.table(file = "~/Dropbox/tmp-for-transfer-to-server/counts-test.txt", sep = "\t", header = T) -> counts

#Create a list of geneIds
genes <- counts$Geneid 
#install biomaRt package in R and run


#to check package is available 
library("package name")


#Setup mart so that biomaRt uses Feb 2021 version of Ensembl
mart103.hs <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "https://feb2021.archive.ensembl.org")

#Search Ensembl 103 and create a separate list 
hgnc.list <- getBM(filters = "ensembl_gene_id", 
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      values = genes, mart = mart103.hs)


#install tidyverse for next step
full_join(hgnc.list, counts, by=c("ensembl_gene_id"="Geneid")) -> counts.with.hgnc
View(counts.with.hgnc)


path=(C:/Users/Anoop Kumar Tiwari/OneDrive/Documents/R/win-library/4.0/biomaRt/html/)


#calculation of RPKM using edgeR 
counts.DGElist <- DGEList(counts=counts.with.hgnc[,8:87], genes=counts.with.hgnc[,1:7])
counts.DGElist.normF <- calcNormFactors(counts.DGElist)
counts.rpkm <- rpkm(counts.DGElist.normF)
class(counts.rpkm)
View(counts.rpkm)


#to merge 2 data frames of hgnc data and rpkm data
total <- merge.data.frame(hgnc.list, counts.rpkm)


#to save a df
save(file.test, file = "data.Rda")


#to write a df from environment


#to drop a column from recent file 
df <- subset(counts.test, select = c(Chr, start))


#to remove columns from file 
within(source-file, rm(Chr, End)) -> output-name


#to find out similarity in files in selected columns
intersect(1stfile$column-name, 2nd-file$colummn-name)


#to find union in 2 columns of different files
setdiff(1stfile$column-name, 2nd-file$colummn-name)


#to select or slice some rows from the df
source-file %>% slice(3:17) -> out-pu-file
https://rdrr.io/bioc/DESeq2/man/fpkm.html



#to run join function with full join
full_join(df-x, df-y, by=c("x-column-name"="y-column-name")) -> output-name
#link for help https://dplyr.tidyverse.org/reference/join.html


#to plot a graph with ggplot (x=should-have-mathematical-values, y=any-column-with-character-also)
ggplot(data = pot, aes(x=Length, y=hgnc_symbol)) + geom_tile(size=10, fill="red")
#geom_tile function in 
ggplot(data = (filter(sialyl_20.scaled, hgnc_symbol == "ST6GALNAC1")), aes(x = samples, y = hgnc_symbol)) + geom_tile(aes(fill=rpkm)) + scale_fill_gradient2(high = "red", mid = "white", low = "blue") + labs(title = "rpkm_samples", x = "Samples", y = "Genes")
#with filled mid v alue
ggplot(data = sialyl_scaled, aes(x = samples, y = hgnc_symbol)) + geom_tile(aes(fill=rpkm)) + labs(title = "rpkm_samples", x = "rpkm_values", y = "hgnc_symbol") + scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0)




#to merge df by cbind function
cbind(file_name, files_name) -> output_name


#to set rownames as 1st column
df <- tibble::rownames_to_column(df, "VALUE")
#to set colnames in my df
colnames(trial2_transpose) <- rownames(trial$hgnc_symbol)
#set colnames by specify column no.
colnames(rawdata_oc)[1] <- "ens_id"
or
#to set 1st row as colname
names(trial2_transpose) <- lapply(trial2_transpose[1, ], as.character)
trial2_transpose <- trial2_transpose[-1,]

#transpose
trial2_transpose <- as.data.frame(t(as.matrix(trial2)))

#to save a file
write.table(trial, file="transposed.txt", row.names=T, sep="\t", quote=F)


#########################################################


library(tidyverse)
Anoop.Scaling.data <- read.table(file = "~/Downloads/transposed.txt", sep = "\t", header = T, row.names = 1)
View(Anoop.Scaling.data)
Anoop.Scaling.data <- read.table(file = "~/Downloads/transposed.txt", sep = "\t", header = T, row.names = 1)
?any
any(is.na(Anoop.Scaling.data))
is.na(Anoop.Scaling.data) %>% View()
dim(Anoop.Scaling.data)
class(Anoop.Scaling.data$OSCC.GB_00490122)
Anoop.Scaling.data %>% t() %>%  scale( center = TRUE, scale = TRUE) %>% t() -> Anoop.Scaling.data.transposed
View(Anoop.Scaling.data.transposed)

# hclust calculation 
Anoop.Scaling.data.transposed[,2:81] %>% t() -> Anoop.Scaling.data.transposed.hclust

Anoop.hclust.order <- hclust(dist(Anoop.Scaling.data.transposed.hclust, 
                                  method = "euclidean"),
                             method = "ward.D")$order
#getting the rownames in an order to be plotted
rownames(Anoop.Scaling.data.transposed.hclust[Anoop.hclust.order,]) -> Anoop.hclust.rownames


# Data Pivotlonger 
pivot_longer(data = Anoop.Scaling.data.transposed,
             !Genes,
             names_to = "samples",
             values_to = "Expression_level") -> Anoop.Scaling.data.transposed.longer


## heatmap using ggplot2
ggplot(Anoop.Scaling.data.transposed.longer, 
       aes(x = factor(samples,
                      levels = Anoop.hclust.rownames),
           y= Genes)) + 
  geom_tile(aes(fill=Expression_level), color = "black") +
  scale_fill_gradientn(colours=c("#5B6F9B", "#BECFDA", "#FBF6C7", 
                                 "#DDB780", "#a0433c"),
                       guide="colorbar") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")


#to clustring
help link:- https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html


#to calculate z-score
help-link:- https://www.r-bloggers.com/2020/02/how-to-compute-the-z-score-with-r/



#####################################################



library(tidyverse)
library(reshape2)
anoop.sialylT.rpkm <- read.table(file = "~/Downloads/sialyl_rpkm.txt", sep = "\t", header = T, row.names = 1)
anoop.sialylT.rpkm <- t(anoop.sialylT.rpkm)
anoop.sialysT.rpkm.cor <- round(cor(anoop.sialylT.rpkm),2)


#### clustering the data ####
#function
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

# reordering 
anoop.sialysT.rpkm.cor <- reorder_cormat(anoop.sialysT.rpkm.cor)



#Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
  
anoop.sialysT.rpkm.cor.upper_tri <- get_upper_tri(anoop.sialysT.rpkm.cor)

#Reshaping data using melt
anoop.sialysT.rpkm.cor.upper_tri.melt <- melt(anoop.sialysT.rpkm.cor.upper_tri, na.rm = TRUE)


#Reshaping the data using tidyverse pivot_longer inspired by this webpage.. still not  working
# https://stackoverflow.com/questions/47475897/correlation-matrix-tidyr-gather-v-reshape2-melt
#anoop.sialysT.rpkm.cor.upper_tri %>%
 # as.data.frame() %>%
 # mutate(Gene1 = factor(row.names(.),
 #                      levels=row.names(.))) %>%
 # pivot_longer(names_to = "Gene2",
 #              values_to = "correlation",
 #              cols = 1:20,
 #              values_drop_na = TRUE) -> anoop.sialysT.rpkm.cor.upper_tri.longer
  
#Rounding of the correlation values - Not needed anymore
#anoop.sialysT.rpkm.cor.upper_tri.melt$value <- round(anoop.sialysT.rpkm.cor.upper_tri.melt$value,2)
  
  
ggheatmap <- ggplot(anoop.sialysT.rpkm.cor.upper_tri.melt, aes(Var1,
                                                      Var2,
                                                      fill = value)) +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue",
                         high = "red",
                         mid = "white", 
                         midpoint = 0,
                         limit = c(-1,1),
                         space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1, 
                                     size = 12,
                                     hjust = 1)) +
    coord_fixed()

ggheatmap + 
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0, 0),
    legend.position = c(0.9, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))




  
