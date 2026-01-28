#Microbial co-occurrence subnetwork
# 清除当前工作环境
rm(list = ls())
gc()
# 定义需要加载的包列表
packages <- c(
  "Hmisc",     
  "vegan",     
  "igraph",    
  "phyloseq",  
  "WGCNA",     
  "microeco",  
  "ggthemes",  
  "ggplot2",   
  "dplyr",     
  "tidyr",     
  "reshape2",  
  "agricolae", 
  "car",       
  "multcompView" 
)
# 批量加载包
lapply(packages, library, character.only = TRUE)
# 验证包是否成功加载
loaded_packages <- packages[packages %in% (.packages())]

# 读取ASV（微生物序列变体）丰度表
otu_table <- read.csv(".csv", header = T, row.names = 1) 

# 读取样本分组信息
sample_info <- read.csv(".csv", header = T, row.names = 1) 

# 读取分类学注释表
otu_tax <- read.csv(".csv", header = T, row.names = 1)
####################################################################
#1.建立微生物网络相关性矩阵
# 数据预处理
df <- otu_table  # 创建数据副本
colSums(df)  # 检查每个样本的测序深度

# 计算相对丰度（标准化）
# 将每列除以该列的总和，得到相对丰度
res1 = t(t(df)/colSums(df))  
colSums(res1)  # 验证相对丰度每列和为1

# 导出相对丰度表
write.table(res1, file =".csv", sep =",", quote =FALSE)  

# 读取相对丰度表
genus <- read.csv('b.csv', row.name = 1)  

# 过滤低丰度微生物（去除丰度低于0.1%的序列）
genus <- genus[which(rowSums(genus) > 0.001), ]    

# 使用Spearman相关性进行微生物间关联分析
genus_corr <- rcorr(t(genus), type = 'spearman')  

# 处理相关性矩阵
# 阈值设置：|r| > 0.8为显著相关
r <- genus_corr$r
r[abs(r) < 0.6] <- 0  # 小于0.8的相关性置零
r[r >= 0.6] <- 1      # 正相关
r[r <= -0.6] <- -1    # 负相关

# p值校正：使用BH方法控制假阳性
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    
p[p>=0.01] <- 0   
p[p<0.01] <- 1    

# 生成网络邻接矩阵
# 相关性矩阵与p值矩阵的乘积
z <- r * p
diag(z) <- 0  # 对角线置零（节点与自身无关联）

# 导出相关性矩阵
write.table(data.frame(z, check.names = FALSE), 'ba2corr.matrix.csv',   
            col.names = NA, sep = ',', quote = FALSE)  

#############################################################
#2.微生物子网络性质计算

# 构建子网络拓扑属性分析函数
process_subnetwork_topology <- function(Genus, g1) {
  # 查找共同的节点名称
  common_nodes <- intersect(V(g1)$name, rownames(Genus))
  
  # 筛选genus数据
  otu <- genus[common_nodes, ]
  
  # 转换为二值矩阵（大于平均丰度的设为1）
  otu_binary <- apply(otu, 2, function(x) {
    # 转换为二值：大于该列平均丰度的设为1
    x > mean(x)
  })
  
  rownames(otu_binary) <- rownames(otu)
  
  sub_graph_stats <- list()
  
  for (i in colnames(otu_binary)) {
    sample_i <- otu_binary[, i]
    select_node <- names(sample_i)[which(sample_i)]
    
    # 仅在节点数大于5时构建子网络
    if (length(select_node) > 5) {
      # 使用共同节点构建子网络
      sub_g <- induced_subgraph(g1, select_node)
      
      sub_graph_stats[[i]] <- data.frame(
        sample_name = i,
        nodes_num = length(V(sub_g)),
        edges_num = length(E(sub_g)),
        clustering_coefficient = transitivity(sub_g, type = "global"),
        degree = mean(degree(sub_g)),
        network_density = edge_density(sub_g, loops = FALSE),
        average_path_length = mean_distance(sub_g, directed = FALSE),
        betweenness_centralization = centr_betw(sub_g)$centralization,
        connectance = edge_density(sub_g, loops = FALSE),
        modularity = tryCatch(
          modularity(sub_g, membership(cluster_walktrap(sub_g))),
          error = function(e) NA
        )
      )
    }
  }
  
  # 合并所有子网络统计
  if (length(sub_graph_stats) > 0) {
    sub_graph_stat <- do.call(rbind, sub_graph_stats)
    return(sub_graph_stat)
  } else {
    warning("No valid subnetworks found")
    return(NULL)
  }
}

# 主分析流程
tryCatch({
  
  # 读取相关性矩阵
  adjacency_unweight <- read.csv('ba2corr.matrix.csv', row.names = 1, check.names = FALSE)
  
  # 处理负值
  adjacency_unweight[adjacency_unweight < 0] <- 1
  
  # 构建无权网络
  g1 <- graph_from_adjacency_matrix(
    as.matrix(adjacency_unweight),   
    mode = 'undirected',   
    weighted = NULL,   
    diag = FALSE
  )
  
  # 读取相对丰度表
  genus <- read.csv('.csv', row.name = 1)
  
  # 过滤低丰度微生物（可根据需要调整阈值）
  genus <- genus[which(rowSums(genus) > 0.001), ]
  
  # 执行子网络拓扑属性分析
  sub_graph_stat <- process_subnetwork_topology(genus, g1)
  
  # 保存详细子网络拓扑属性
  if (!is.null(sub_graph_stat)) {
    # 保存详细子网络拓扑属性
    write.csv(sub_graph_stat, '子网络拓扑属性.csv', quote = FALSE, row.names = FALSE)
    # 输出处理摘要
    cat("分析完成！\n")
    cat("详细子网络拓扑属性已保存到 '子网络拓扑属性.csv'\n")
  }
})

