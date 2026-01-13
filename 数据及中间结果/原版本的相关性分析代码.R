# ==========================================================
# Phase 3 (升级版): N=57 深度相关性分析 & 带星号棒棒糖图
# ==========================================================
library(dplyr)
library(ggplot2)
library(readr)
library(stringr) # 用于处理长字符串

# 1. 检查数据
if (!exists("analysis_data") || !exists("candidates")) {
  stop("数据缺失！请先运行 Phase 1 和 Phase 2 的代码。")
}

cat(">>> 开始计算相关性并生成星号...\n")

# 定义一个函数：将 P值 转换为 星号
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("") # 不显著的不标
}

# 2. 批量计算相关性
cor_data <- data.frame(
  Display_Name = character(),
  Rho = numeric(),
  P_Value = numeric(),
  Stars = character(),     # 新增：星号列
  Type = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(candidates)) {
  comp_name <- candidates$Computational_Name[i]
  disp_name <- candidates$Display_Name[i]
  
  if (comp_name %in% colnames(analysis_data)) {
    # 计算 Spearman
    res <- cor.test(analysis_data[[comp_name]], analysis_data$Sourness_Score, 
                    method = "spearman", exact = FALSE)
    
    rho <- res$estimate
    p   <- res$p.value
    
    # 新增：根据 Rho 的正负设定文字标签的位置偏移量 (为了画图好看)
    # 正相关星号标在右边，负相关标在左边
    
    cor_data <- rbind(cor_data, data.frame(
      Display_Name = disp_name,
      Rho = rho,
      P_Value = p,
      Stars = get_stars(p), # 转换星号
      Type = ifelse(rho > 0, "Positive (Promoting)", "Negative (Suppressing)")
    ))
  }
}

# 3. 整理排名 (只保留 P < 0.05)
final_ranking <- cor_data %>%
  filter(P_Value < 0.05) %>%
  arrange(desc(abs(Rho)))

# 导出 CSV (现在包含星号列了！)
write_csv(final_ranking, "3_N57_Correlation_Ranking_WithStars.csv")
cat(">>> 表格已更新！请查看 '3_N57_Correlation_Ranking_WithStars.csv' (包含星号)\n")

# ----------------------------------------------------------
# ==========================================================
## 4. 基于专家经验的代表性物质绘图-棒棒糖图
# ==========================================================

# 1. 定义专家筛选名单 (The "Golden List")
# 这些是我们精挑细选的 16 个代表性物质
expert_list <- c(
  # --- 1. 酸类 (Acids) ---
  "乙酸", "丙酸", "乳酸",         # 促进酸感的挥发酸/半挥发酸
  "柠檬酸", "草酸", "3-甲基戊酸", # 负相关的酸 (缓冲或杂气)
  
  # --- 2. 糖类 (Sugars) ---
  "总糖", "还原糖", "蔗糖",       # 甜味物质 (平衡酸感)
  
  # --- 3. 碱/氮类 (Alkalis) ---
  "尼古丁", "总氮",               # 碱性物质 (中和酸感)
  "糖碱比",                       # 关键平衡指标
  
  # --- 4. 关键前体/多酚 (Precursors) ---
  "绿原酸", "芸香苷",             # 多酚 (相关性Top)
  "DDMP",                         # 焦甜香/烘烤产物
  "脯氨酸"                        # 氨基酸代表
)

# 假设您的相关性结果表叫 cor_data (包含 Display_Name, Rho, Stars, Type, P_Value)
# 如果您是直接接着上一段代码跑，这个 cor_data 应该还在内存里

# 2. 筛选数据
plot_data <- cor_data %>%
  filter(Display_Name %in% expert_list) %>%
  arrange(desc(abs(Rho))) # 依然按相关性强度排序，画出来好看

# 3. 绘图 (样式优化)
# 增加文字换行，防止名字太长
plot_data$Display_Name <- str_wrap(plot_data$Display_Name, width = 20)

p <- ggplot(plot_data, aes(x = reorder(Display_Name, Rho), y = Rho)) +
  # 画线
  geom_segment(aes(xend = reorder(Display_Name, Rho), y = 0, yend = Rho, color = Type), 
               linewidth = 1.2) +
  # 画点
  geom_point(aes(color = Type, size = abs(Rho)), alpha = 1) +
  
  # 标星号 (位置自动偏移)
  geom_text(aes(label = Stars, y = Rho), 
            nudge_y = ifelse(plot_data$Rho > 0, 0.06, -0.06), # 稍微调大一点偏移量
            vjust = 0.75, size = 5, fontface = "bold") +
  
  # 颜色和标签
  scale_color_manual(values = c("Positive (Promoting)" = "#E41A1C", 
                                "Negative (Suppressing)" = "#377EB8")) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  labs(title = "Representative Sourness-Correlated Compounds",
       subtitle = "Selected Markers: Acids, Sugars, Alkaloids & Polyphenols",
       x = "", y = "Spearman Correlation Coefficient (Rho)",
       color = "Correlation Type", size = "Strength (|Rho|)") +
  
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom",
    plot.margin = unit(c(1,1,1,1), "cm")
  )

# 保存
ggsave("3_Selected_Representatives_Lollipop.png", p, width = 8, height = 8)
print(p)