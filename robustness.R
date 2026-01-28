#robustness
rm(list = ls())
library(microeco)
library(meconetcomp)
library(magrittr)
library(dplyr)
setwd("")

# 导入数据
sample_info_16S <- read.csv(".csv", row.names = 1)
otu_table_16S <- read.csv(".csv", row.names = 1)

otu2 <- apply(otu_table_16S, 2, function(x) x / sum(x))
sample_presence <- apply(otu2, 1, function(x) sum(x > 0))
otu3 <- otu2[sample_presence >= 3, ]
otu_table_16S <- as.data.frame(otu3)

taxonomy_table_16S <- read.csv(".csv", row.names = 1)
taxonomy_table_16S <- taxonomy_table_16S[rownames(otu_table_16S), ]

# 构建 microtable
soil_amp <- microtable$new(
  sample_table = sample_info_16S,
  otu_table = otu_table_16S,
  tax_table = taxonomy_table_16S
)

soil_amp_network <- list()

## ---------------------
## A 组网络
## ---------------------
tmpA <- clone(soil_amp)
tmpA$sample_table %<>% subset(group == "A")
tmpA$tidy_dataset()

net_A <- trans_network$new(
  dataset = tmpA,
  cor_method = "spearman",
  filter_thres = 0.000001
)
net_A$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)

soil_amp_network$A <- net_A

## ---------------------
## B 组网络
## ---------------------
tmpB <- clone(soil_amp)
tmpB$sample_table %<>% subset(group == "B")
tmpB$tidy_dataset()

net_B <- trans_network$new(
  dataset = tmpB,
  cor_method = "spearman",
  filter_thres = 0.000001
)
net_B$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)

soil_amp_network$B <- net_B

## ---------------------
## 鲁棒性分析（只返回 A 与 B）
## ---------------------
rob <- robustness$new(
  soil_amp_network,
  remove_ratio = seq(0, 0.99, 0.1),
  measure = c("Eff", "Eigen", "Pcr"),
  run = 10
)

View(rob$res_table)
write.csv(rob$res_table, ".csv")
View(rob$res_summary)
write.csv(rob$res_summary, ".csv")
library(ggplot2)
g1 <- rob$plot(linewidth = 1) + scale_color_brewer(palette = "Set3")
g1
vul_table <- vulnerability(soil_amp_network)# 计算易损性
View(vul_table)
write.csv(vul_table, ".csv")
result_table <- vul_table %>%
  group_by(Network) %>%
  slice_max(vulnerability, n = 1) %>%
  ungroup()
rite.csv(result_table, ".csv")
t1 <- cohesionclass$new(soil_amp_network)
View(t1$res_list$sample)
write.csv(t1$res_list$sample,".csv")
View(t1$res_list$feature)
write.csv(t1$res_list$feature,".csv")
net<- t1$res_list$sample
library(dplyr)
net <- net%>%mutate(stability=abs(c_neg)/c_pos)
write.csv(net,".csv")
dir.create("D0_results", showWarnings = FALSE, recursive = TRUE)
