# ==========================================================
#  数据体检自测脚本 (No Black Box)
# ==========================================================
library(readr)
library(dplyr)
library(tidyr)
library(rstatix) # 用于正态性检验

# 1. 读取数据
df <- read_csv("论文1，5.0/数据及中间结果/data-18.csv") 
num_cols <- df %>% select(where(is.numeric)) %>% colnames()


# ----------------------------------------------------------
# A. 正态性检验 (Shapiro-Wilk 方法)
# 原理：P > 0.05 代表正态；P < 0.05 代表不正态
# ----------------------------------------------------------
normality_report <- data.frame(Variable = num_cols, P_value = NA, Is_Normal = NA)

for (i in 1:length(num_cols)) {
  col <- num_cols[i]
  vals <- df[[col]]
  
  # Shapiro检验要求样本量 3-5000
  if (length(na.omit(vals)) >= 3) {
    res <- shapiro_test(vals)
    normality_report$P_value[i] <- res$p.value
    normality_report$Is_Normal[i] <- ifelse(res$p.value > 0.05, "Yes", "No")
  }
}

# 保存正态性报告
View(normality_report)



# ----------------------------------------------------------
# B. 鲁棒离群值检测 (IQR / Boxplot Method)
# 原理：不看平均值，看中位数。
# 标准：超过 (Q3 + 3*IQR) 或 低于 (Q1 - 3*IQR) 为极端异常
# ----------------------------------------------------------
robust_report <- data.frame()

for (col in num_cols) {
  vals <- df[[col]]
  
  # 计算四分位数 (25% 和 75%)
  Q1 <- quantile(vals, 0.25, na.rm = TRUE)
  Q3 <- quantile(vals, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  median_val <- median(vals, na.rm = TRUE)
  
  # 定义极端异常的界限 (Extreme Outlier Limits)
  # 这里的 3 倍 IQR 相当于抓 "极度离谱" 的值
  # 如果想抓得严一点，可以把 3 改成 1.5
  upper_limit <- Q3 + 3 * IQR_val
  lower_limit <- Q1 - 3 * IQR_val
  
  # 找到异常值的位置
  bad_idx <- which(vals > upper_limit | vals < lower_limit)
  
  if (length(bad_idx) > 0) {
    temp <- data.frame(
      Variable = col,
      SampleID = df$'样品编号'[bad_idx],
      Value = vals[bad_idx],
      Median = median_val,
      # 计算偏离程度 (倍数)：(值 - 中位数) / IQR
      # 这个值越大，说明越离谱
      Deviation_Score = round((vals[bad_idx] - median_val) / IQR_val, 2)
    )
    # 过滤掉 Deviation_Score 无穷大或太小的情况
    temp <- temp %>% filter(is.finite(Deviation_Score) & abs(Deviation_Score) > 0)
    
    robust_report <- rbind(robust_report, temp)
  }
}

# 按离谱程度排序 (从大到小)
robust_report <- robust_report %>% arrange(desc(Deviation_Score))

# 查看
View(robust_report)





