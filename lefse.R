#lefse
rm(list=ls())

library(tidyverse)
library(microeco)
library(magrittr)
library(patchwork)
library(tidygraph)
library(ggfun)
library(ggplot2)

setwd("")

otu <- read.csv("", row.names = 1)
group <- read.csv("", row.names = 1)
tax <- read.csv("", row.names = 1)

if (!all(rownames(group) == colnames(otu))) {
  stop("The row names of sample_table are not sample names!")
}
colnames(otu) <- rownames(group)

dataset <- microtable$new(
  sample_table = group,
  otu_table = otu,
  tax_table = tax
)

set.seed(123)

t1 <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  filter_thres = 0.00001,
  alpha = 0.05,
  p_adjust_method = "none",
  group = "group",
  taxa_level = "genus"
)

## ===== 差异丰度图（可选）=====
t1$plot_diff_abund(
  use_number = 1:30,
  color_values = c("#84c8B7", "#F57C6E"),
  simplify_names = TRUE,
  keep_prefix = FALSE,
  order_x_mean = TRUE,
  coord_flip = TRUE,
  add_sig = TRUE,
  xtext_angle = 45,
  xtext_size = 13,
  ytitle_size = 17
)

## ===== LEfSe LDA bar（关键图）=====
t1$plot_diff_bar(
  color_values = c("#84c8B7", "#F57C6E"),
  color_group_map = TRUE,
  use_number = 1:30,     # 显示 Top 30
  threshold = 2,        #LDA > 2
  keep_full_name = FALSE,
  keep_prefix = FALSE,
  group_aggre = TRUE,
  group_two_sep = TRUE,
  coord_flip = TRUE,
  add_sig = TRUE,
  add_sig_increase = 0.1,
  add_sig_text_size = 5,
  xtext_angle = 45,
  xtext_size = 10,
  axis_text_y = 12,
  heatmap_cell = "P.unadj",
  heatmap_sig = "Significance",
  heatmap_x = "Factors",
  heatmap_y = "Taxa",
  heatmap_lab_fill = "P value"
)
