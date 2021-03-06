#source("http://bioconductor.org/biocLite.R")   
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages("gplots")

foldChange=1
padj=0.05

setwd("")                    #设置工作目录
library("edgeR")
rt=read.table("genesymbolmatrix.txt",sep="\t",header=T,check.names=F)  #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

#group=c("normal","tumor","tumor","normal","tumor")
group=c(rep("normal",4),rep("tumor",178)) #按照癌症和正常样品数目修改
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）

heatmapData <- newData[rownames(diffSig),]

#volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

#heatmap
hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",width=60,height=90)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()


DESeq
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
#biocLite("limma")
#install.packages("gplots")

foldChange=1
padj=0.05


setwd("") #设置工作目录（需修改）
library("DESeq")
library("limma")

rt=read.table("genesymbolmatrix.txt",sep="\t",header=T,check.names=F)  #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,] 
data=round(data,0)

group=c(rep("normal",4),rep("tumor",178))   #按照癌症和正常样品数目修改
design = factor(group)
newTab = newCountDataSet( data, design )
newTab = estimateSizeFactors(newTab)
newData=counts(newTab, normalized=TRUE )

#have replicates
newTab = estimateDispersions( newTab, fitType = "local")
diff = nbinomTest( newTab, "normal", "tumor")
diff = diff[is.na(diff$padj)==FALSE,]
diff = diff[order(diff$pval),]
write.table( diff, file="DESeqOut.xls",sep="\t",quote=F,row.names=F)
diffSig = diff[(diff$padj < padj & (diff$log2FoldChang>foldChange | diff$log2FoldChange<(-foldChange))),]
write.table( diffSig, file="diffSig.xls",sep="\t",quote=F,row.names=F)
diffUp = diff[(diff$padj < padj & (diff$log2FoldChang>foldChange)),]
write.table( diffUp, file="up.xls",sep="\t",quote=F,row.names=F)
diffDown = diff[(diff$padj < padj & (diff$log2FoldChange<(-foldChange))),]
write.table( diffDown, file="down.xls",sep="\t",quote=F,row.names=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）
diffExp=rbind(id=colnames(newData),newData[diffSig$id,])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmiRNAExp.txt）

#heatmap
hmExp=log10(newData[diffSig$id,]+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",width=60,height=30)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()

#volcano
pdf(file="vol.pdf")
allDiff=diff[is.na(diff$padj)==FALSE,]
xMax=max(-log10(allDiff$padj))+1
yMax=10
plot(-log10(allDiff$padj), allDiff$log2FoldChange, xlab="-log10(padj)",ylab="log2FoldChange",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange>foldChange,]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange<(-foldChange),]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()


## 提取生存分析中的基因表达量
#install.packages("hash")

gene='BRCA1'   #基因名（需修改）
clinicalFile="time.txt"
expFile="normalizeExp.txt"

library(hash)
setwd("")     #工作目录（需修改）
rt=read.table(clinicalFile,header=T,check.names=F,sep="\t")
h = hash(keys = rt$id, values = paste(rt$futime,rt$fustat,sep="\t"))

exp=read.table(expFile,header=T,check.names=F,row.names=1,sep="\t")
geneExp=t(exp[gene,])
write.table("sample\tfutime\tfustat\texpression",file="survivalInput.txt",sep="\t"
  ,quote=F,row.names=F,col.names=F)

for(i in rownames(geneExp)){
  j=unlist(strsplit(i,"\\-"))
  if(grepl("^0",j[4])){
    name4=paste(j[1],j[2],j[3],j[4],sep="-")
    name3=paste(j[1],j[2],j[3],sep="-")
    if(has.key(name3,h)){
      write.table(paste(name4,h[[name3]],geneExp[i,],sep="\t"),file="survivalInput.txt",sep="\t",
                  quote=F,append=TRUE,row.names=F,col.names=F)
                  }
    }
}


setwd("C:\\Users\\DELL\\Desktop\\TCGA\\time") 

library(survival)
rt=read.table("survivalInput.txt",header=T,sep="\t")
#rt$futime=rt$futime/365       #如果以月为单位，除以30；以年为单位，除以365
a=rt[,"expression"]<median(rt[,"expression"])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,5)
fit <- survfit(Surv(futime, fustat) ~ a, data = rt)

summary(fit)    #查看五年生存率
pdf(file="survival.pdf")
plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
     main=paste("surival curve (p=", pValue ,")",sep=""))
legend("topright", c("BRCA1 high expression", "BRCA1 low expression"), lty = 2:3, col=c("red","blue"))
dev.off()