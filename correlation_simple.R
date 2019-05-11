A1 <- read.csv("COAD_Level_3__RSEM_genes_normalized_combined.csv", header = TRUE)

rownames(Liver) = make.names(Liver$ï..X, unique=TRUE) ## if duplicate use this one
B1=A1[,c(2):ncol(A1)] ## Remove the fist row

data <- read.csv("norm_counts_COAD_mRNA.csv", header=TRUE)
row.names(Liver)=Liver$ï..ID_REF ## Set column 1 as the row names
Liver$ï..ID_REF <- NULL ## Remove column by a col name

B2<- t(B2)
B2 <- as.data.frame(B2)
B3 <- cor(B2, method="pearson") ## "kendall", "spearman", "pearson"

B3 <- t(B3)
write.csv(B3, "correlation.csv")

B1 <- log2(A1)
B2 <- do.call(data.frame,lapply(B1, function(x) replace(x, is.infinite(x),0)))
write.csv(B3, "LIHC_hiseq_RSEM_genes_log2_normalized.csv" )
