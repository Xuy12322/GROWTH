#zipi

rm(list=ls())#好习惯，确保有干净的 R 环境
library(ggClusterNet)#加载包

library(igraph)#后读入edge和node文件
edge <- read.csv(".csv",header = T, sep = ",", check.names = FALSE)
node <- read.csv(".csv",header = T, sep = ",", check.names = FALSE)

#从edge、node的数据框中创建igraph对象
igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
igraph <- simplify(igraph, remove.multiple = TRUE, remove.loops = TRUE)
#ZiPiplot
res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")

p <- res[[1]]# res[[1]]是图
p
