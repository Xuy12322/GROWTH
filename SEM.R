#SEM
rm(list=ls())
#未安装lavaan包可以释放下面的安装命令
#install.packages("lavaan")
library(lavaan)
library(haven)
library(Hmisc)
library(openxlsx)
#读取数据
#读取数据
setwd("")
mydata11 <- read.csv(".csv", header = TRUE)#linear_mixed model data
head(mydata11)
mydata11 <- mydata11[c("Maizebiomass","AMFdiversity","Bacterialdiversity","BacteriaNetworkcomplexity","Myceliumdensity")]
#对数据进行标准化
mydata<-data.frame(scale(mydata11))
head(mydata)
model3 <- '
Bacterialdiversity~ AMFdiversity
BacteriaNetworkcomplexity~  AMFdiversity+Bacterialdiversity
Myceliumdensity~ AMFdiversity+Bacterialdiversity+BacteriaNetworkcomplexity
Maizebiomass~  Myceliumdensity+Bacterialdiversity+BacteriaNetworkcomplexity
'
#拟合模型
fit1 <- sem(model3,data = mydata)

#模型摘要结果
summary(fit1,standardized = TRUE,fit.measures = TRUE,rsquare = T)#fit.measures = TRUE
fitMeasures(fit1,c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea","AIC"))
standardizedSolution(fit1)
