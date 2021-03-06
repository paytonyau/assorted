## install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("affy", "affyPLM","RColorBrewer","impute","limma","pheatmap"))
install.packages("ggplot2")

library(affyPLM)

Data<-ReadAffy() ## input data
Pset<-fitPLM (Data) ##对数据集进行回归计算
image(Data[,1]) ## 质量控制：查看灰度图
image(Pset,type="weights",which=1,main="Weights") ##根据计算结果，画权重图
image(Pset,type="resids",which=1,main="Residuals") ##根据计算结果，画残差图
image(Pset,type="sign.resids",which=1,main="Residuals.sign")##根据计算结果，画残差符号图


## 质量控制:相对对数表达（RLE）
## 一个探针组在某个样品的表达值除以该探针组在所有样品中表达之的中位数后取对数
## 反映平行实验的一致性
#调用R包
library(affyPLM)
library(RColorBrewer)

Pset<-fitPLM (Data) #对数据集进行回归计算
colors<-brewer.pal(12,"Set3") #载入颜色
Mbox(Pset,col=colors,main="RLE",las=3) #绘制RLE箱线图

## 质量控制：相对标准差（NUSE）
## 一个探针组在某个样品的PM值的标准差除以该探针组在各样品中的PM值标准差的中位数后取对数。
##反映平行实验的一致性比RLE更为敏感。

library(affyPLM)
library(RColorBrewer)

Pset<-fitPLM (Data) #对数据集进行回归计算
colors<-brewer.pal(12,"Set3") #载入颜色
boxplot(Pset,col=colors,main="NUSE",las=3) #绘制NUSE箱线图

## 质量控制：RNA降解图
## 原理：RNA降解从5’端开始，因为芯片结果5’端荧光强度要远低于3’端

#调用R包
library(affy)

data.deg<-AffyRNAdeg(Data) #获取降解数据
plotAffyRNAdeg(data.deg,col=colors) #绘制RNA降解图
legend("topleft",sampleNames(Data),col=colors,lwd=1,inset=0.05,cex=0.2) #在左上部位添加图注


## RMA法预处理normal样本
setwd("")
library(affyPLM)
library(affy)

Data<-ReadAffy()
sampleNames(Data)
N=length(Data)

eset.rma<-rma(Data) ##用RMA预处理数据
##获取表达数据并输出到表格
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="normal.expres.txt",sep='\t',quote=F,row.names=F)

## RMA法预处理tumor样本
setwd("")
library(affyPLM)
library(affy)

Data<-ReadAffy()
sampleNames(Data)
N=length(Data)

## 用RMA预处理数据
eset.rma<-rma(Data)
## 获取表达数据并输出到表格
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="tumor.expres.txt",sep='\t',quote=F,row.names=F)

## 合并N和T的数据
#setwd(" ")
normal_exprs<-read.table("normal.expres.txt",header=T,sep="\t")
tumor_exprs<-read.table("tumor.expres.txt",header=T,sep="\t")
#讲T和N合并
probe_exprs<-merge(normal_exprs,tumor_exprs,by="probeid")
write.table(probe_exprs,file="cancer.probeid.exprs.txt",sep='\t',quote=F,row.names=F)

## Probe ID conversion Gene symbol（对平台文件进行整理）
setwd("")

probe_exp<-read.table("cancer.probeid.exprs.txt",header=T,sep="\t",row.names=1) #读取基因表达文件
probeid_geneid<-read.table("GPL570-55999.txt",header=T,sep="\t") #读取探针文件
probe_name<-rownames(probe_exp)
loc<-match(probeid_geneid[,1],probe_name) #probe进行匹配
probe_exp<-probe_exp[loc,] #确定能匹配上的probe表达值
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3])) #每个probeid应对的geneid
index<-which(!is.na(raw_geneid)) #找出有geneid的probeid并建立索引

geneid<-raw_geneid[index] #提取有geneid的probe

#找到每个geneid的表达值
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)

#多个探针对应1个基因的情况，取平均值
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))

#geneid作为行名
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="Gastric.cancer.geneid.exprs.txt",sep='\t',quote=F,row.names=F)

#将gene id 转换为gene symbol
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="Gastric.cancer.genesyb.exprs.txt",sep='\t',quote=F,row.names=F)

## 补充缺失值 (需要对genesyb这个文件需要处理)
## 最近邻居法（KNN，k-Nearest Neighbor）法：
## 此方法是寻找和有缺失值的基因的表达谱相似的其他基因，
## 通过这些基因的表达值（依照表达谱相似性加权）来填充缺失值
#调用函数impute

library(impute)

#读取表达值
gene_exp_matrix<-read.table("Gastric.cancer.genesyb.exprs.txt",header=T,sep="\t",row.names=1)
gene_exp_matrix<-as.matrix(gene_exp_matrix)

#KNN法计算缺失值
imputed_gene_exp<-impute.knn(gene_exp_matrix,k=10,rowmax=0.5,colmax=0.8,maxp=3000,rng.seed=362436069)

#读出经过缺失值处理的数据
GeneExp<-imputed_gene_exp$data

#写入表格
genesymbol<-rownames(GeneExp)
GeneExp<-cbind(genesymbol,GeneExp)
write.table(GeneExp,file="Gastric.cancer.gene.exprs.txt",sep='\t',quote=F,row.names=F)


library(limma)
rt<-read.table("Gastric.cancer.gene.exprs.txt",header=T,sep="\t",row.names="genesymbol")

#differential
class<-c(rep("normal",10),rep("tumor",10))
design<-model.matrix(~factor(class))
colnames(design)<-c("normal","tumor")
fit<-lmFit(rt,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffLab<-allDiff[with(allDiff, ((logFC>1 |logFC<(-1)) & adj.P.Val<0.05)),]
write.table(diffLab,file="diffEXp.xls",sep="\t",quote=F)

## co-expression
diffExpLevel<-rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLevel.xls",sep="\t",quote=F)

## Heatmap
hmExp=log10(diffExpLevel+0.00001)
library('gplots')

hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",height=120,width=90)
par(oma=c(3,3,3,5))
heatmap.2(hmMat,col='greenred',trace="none",cexCol=1)
dev.off()

## Volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$adj.P.Val),allDiff$logFC,xlab="adj.P.Val",ylab="logFC",main="Volcano",xlim=c(0,xMax),ylim=c(-yMax,yMax),pch=20,cex=0.4)
diffSub=subset(allDiff,allDiff$adj.P.Val<0.05 & abs(allDiff$logFC)>1)
points(-log10(diffSub$adj.P.Val),diffSub$logFC,pch=20,col="red",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()


BiocManager::install(c("clusterProfiler","pathview"))
setwd("C:\\Users\\DELL\\Desktop\\mo")

library("clusterProfiler")
rt=read.table("X.txt",sep="\t",head=T,check.names=F)
geneFC=rt$logFC
gene<-rt$ENTREZ_GENE_ID
names(geneFC)=gene

#kegg
kk<-enrichKEGG(gene=gene,organism="human",qvalueCutoff=0.05,readable=TRUE)
write.table(summary(kk),file="KEGG.xls",sep="\t",quote=F,row.names=F)
pdf(file="KEGG.barplot.pdf")
barplot(kk,drop=TRUE,showCategory=12)
pdf(file="KEGG.cnetplot.pdf")
cnetplot(kk,categorySize="geneNum",foldChange=geneFC)   

library("pathview")
keggxls=read.table("KEGG.xls",sep="\t",header=T)

for(i in keggxls$ID){
pv.out<-pathview(gene.data=geneFC,pathway.id=i,species="hsa",out.suffix="pathview")}

