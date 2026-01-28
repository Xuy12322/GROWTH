# Co-occurrence network
#清除环境变量数据
rm(list=ls())
#安装并加载包 
package.list=c("psych","reshape2","Hmisc","ggplot2")
for (package in package.list) {
  if (!require(package,character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


#读取otu数据,这里去除线粒体或叶绿体(细菌去除线粒体,真菌去除叶绿体)
OTU <- read.table(file.choose(), header=T, row.names=1)
#读取otu表注释信息
tax <- read.table(file.choose(), sep="\t", header=T)
names(tax)[1] <- "Id"

#例如只保留在 5 个及以上样本中出现的ASV
OTU1 <- OTU
OTU1[OTU1>0] <- 1
OTU <- OTU[which(rowSums(OTU1) >= 3), ]
OTU = t(OTU)

#设置分析阈值,默认保留丰度>0.05%的OTU；
abundance=0.05
#相对丰度表
OTU <- OTU[,colSums(OTU)/sum(OTU)>=(abundance/100)]


#重要:网络分析的关联阈值(自行可以设置)
r.cutoff=0.6
p.cutoff=0.01
#检验,这么采用spearman,可适应非正态数据,一般微生物数据不符合正态的，要么分析前进行检验，采用spearman不会有问题。
occor=corr.test(OTU, use="pairwise", method="spearman", adjust="fdr", alpha=0.05)
#提取相关矩阵的r、p值；
r_matrix=occor$r
p_matrix=occor$p
#在转换为长格式前，将相关矩阵的对角线元素设置为0
diag(r_matrix) <- NA
#确定物种间存在相互关系的阈值，将相关性R矩阵不符合的数据转换为0；
r_matrix[p_matrix>p.cutoff|abs(r_matrix)<r.cutoff]=NA
#转换数据为长格式形式，方便下游分析；
p_value=melt(p_matrix)
r_value=melt(r_matrix)
#将r、p两个表合并；
r_value=cbind(r_value, p_value$value)
#删除含有r_value=0的行；
r_value=subset(r_value, r_value[,3]!=0)
#删除含有r_value=NA的行；
r_value=na.omit(r_value)
#对r表格增补绝对值、正负型等信息
abs=abs(r_value$value)
linktype=r_value$value
linktype[linktype>0]=1
linktype[linktype<0]=-1
r_value=cbind(r_value, abs, linktype)
#重命名r、p表头
names(r_value) <- c("Source","Target","r_value","p_value", "abs_value", "linktype")
names(p_value) <- c("Source","Target","p_value")
#获取节点数据
node_OTU <- as.data.frame(as.data.frame(r_value[,1])[!duplicated(as.data.frame(r_value[,1])), ])
names(node_OTU)="Id"
#OTU ID合并生成节点索引表,用于检索注释信息;
list <- rbind(node_OTU)
#筛选节点对应的注释信息;
list=subset(tax,Id %in% list$Id)
#复制一列当节点Lable
list$Label <- list$Id

write.csv(list,file=".csv", row.names=FALSE)
write.csv(r_value, file=".csv", row.names=FALSE)
