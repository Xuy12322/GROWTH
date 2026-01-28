#OLS regression model
library(ggplot2)
library(ggpubr)
library(dplyr)

#============================
# 读取数据
#============================
data <- read.csv(".csv")

#============================
# 定义绘图函数
#============================
plot_corr <- function(data, xvar, yvar = "Maizebiomass") {
  # 构造公式
  formula <- as.formula(paste(yvar, "~", xvar))
  
  # 拟合线性模型
  model <- lm(formula, data = data)
  r2 <- summary(model)$r.squared
  p  <- summary(model)$coefficients[2, 4]
  
  # 格式化标注
  label_text <- paste0(
    "R² = ", round(r2, 3),
    "\nP = ", ifelse(p < 0.001, "< 0.001", round(p, 3))
  )
  
  # 绘图
  ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    annotate(
      "text",
      x = min(data[[xvar]], na.rm = TRUE),
      y = max(data[[yvar]], na.rm = TRUE),
      hjust = 0, vjust = 1,
      label = label_text,
      size = 5,
      color = "black",
      fontface = "bold"
    ) +
    labs(
      x = xvar,
      y = yvar,
      title = paste("Relationship between", yvar, "and", xvar)
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

#============================
# 批量绘图
#============================
xvars <- c("", "", "", "")

plots <- lapply(xvars, function(x) {
  plot_corr(data, xvar = x)
})

#============================
# 合并成一张图
#============================
ggarrange(
  plotlist = plots,
  ncol = 2,
  nrow = 2,
  labels = c("(a)", "(b)", "(c)", "(d)")
)
