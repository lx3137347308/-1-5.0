这是一份为您定制的 **《加热卷烟酸感物质基础分析指南 v5.0 (旗舰精简版)》**。

这份指南严格遵照您的最新策略：**摒弃复杂的机器学习黑箱，回归经典统计学**。同时，我保留了您在 v3.3 版本中辛苦调试出的**“中文支持”、“特殊字符清洗映射”**以及**“发表级绘图美学”**代码。

您可以直接将以下内容保存为 `Analysis_Guide_v5.0.md` 使用。

---
  
  # 加热卷烟酸感物质基础分析指南 v5.0 (旗舰精简版)
  
  **版本核心策略：**
  
  1. **分组策略**：将连续评分转化为 4 个离散等级（微酸/稍酸/酸/很酸）。
2. **筛选逻辑**：漏斗式筛选 —— 差异显著 ()  相关显著 ()  交集确证。
3. **核心算法**：Kruskal-Wallis 非参数检验 + Spearman 相关性。
4. **功能保留**：保留了列名清洗/还原映射机制，以及中文字体支持。

---
  
  ## 〇、项目准备：环境与配置
  
  此部分代码整合了 v3.3 的核心配置，确保中文不乱码、特殊符号不报错。

```r
# ==========================================================
# Phase 0: 程序包加载与环境配置
# ==========================================================
# 核心数据处理
library(readr)       # 读取数据
library(dplyr)       # 数据清洗管道
library(tibble)      # 数据框操作
library(tidyr)       # 数据整洁化

# 统计分析
library(rstatix)     # 现代统计检验 (KW test等)
library(agricolae)   # 用于多重比较标记 (a, b, ab...)

# 可视化工具
library(ggplot2)     # 绘图核心
library(ggpubr)      # 专为发表设计的绘图包 (ggboxplot, ggbarplot)
library(pheatmap)    # 热图
library(RColorBrewer)# 配色方案
library(ggrepel)     # 避免标签重叠
library(showtext)    # 中文支持

# --- 字体配置 (保留 v3.3 设置) ---
# 自动搜索系统中的宋体 (Windows通常为 simsun.ttc, Mac可能不同)
# 如果报错，请根据系统实际字体路径修改
font_add("songti", "simsun.ttc") 
showtext_auto() 

# 设置全局绘图主题，移除网格线，统一字体
theme_set(theme_bw(base_family = "songti") + 
            theme(panel.grid = element_blank()))

```

---
  
  ## 第一步：数据准备与清洗 (保留映射机制)
  
  此处完美保留了您“先改名计算，再改回名画图”的逻辑，解决 `1-` 或 `()` 等符号导致的报错问题。

```r
# ==========================================================
# Phase 1: 数据加载、清洗与映射表构建
# ==========================================================

# 1.1 读取原始数据 (假设文件名为 data.csv)
# 请确保第一列是 SampleID，第二列是 Sourness_Score
raw_data <- read_csv("data.csv", locale = locale(encoding = "UTF-8"))

# 1.2 构建“计算名-显示名”映射表 (解决特殊字符报错问题)
# 获取原始列名 (显示用)
original_names <- colnames(raw_data)

# 生成清洗后的列名 (计算用：替换空格、括号、减号为下划线，数字开头加X)
clean_names <- make.names(original_names, unique = TRUE)

# 创建并保存映射字典
mapping_df <- tibble(
  Display_Name = original_names,
  Computational_Name = clean_names
)
# 将数据框的列名替换为“计算名”
colnames(raw_data) <- clean_names

# 1.3 数据分组 (Discretization)
# 根据 Sourness_Score 进行硬性分组
# 标准: 0-1.5(微酸), 1.5-3.0(稍酸), 3.0-4.5(酸), 4.5-5.0(很酸)
analysis_data <- raw_data %>%
  mutate(Group = cut(Sourness_Score, 
                     breaks = c(-Inf, 1.5, 3.0, 4.5, Inf),
                     labels = c("Trace", "Mild", "Moderate", "Strong"),
                     include.lowest = TRUE))

# 定义物质变量的范围 (假设从第3列开始是化学物质)
# 请根据实际情况修改列索引，这里假设 SampleID, Sourness_Score, Group 之后全是物质
chemical_vars <- clean_names[which(!clean_names %in% c("SampleID", "Sourness_Score", "Group"))]

cat("数据准备完成。已创建分组变量 'Group'。\n")
write_csv(mapping_df, "name_mapping_v5.csv") # 备份映射表

```

---
  
  ## 第二步：描述性统计 (新增需求)
  
  按照学术规范，输出各组的均值和标准差。

```r
# ==========================================================
# Phase 2: 描述性统计 (Descriptive Statistics)
# ==========================================================

desc_stats <- analysis_data %>%
  group_by(Group) %>%
  summarise(across(all_of(chemical_vars), 
                   list(Mean = ~mean(., na.rm = TRUE), 
                        SD = ~sd(., na.rm = TRUE)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(cols = -Group, names_to = c("Computational_Name", "Stat"), names_sep = "__") %>%
  pivot_wider(names_from = c(Group, Stat), values_from = value) %>%
  left_join(mapping_df, by = "Computational_Name") %>%
  select(Display_Name, everything(), -Computational_Name)

# 保存描述性统计表
write_csv(desc_stats, "Table1_Descriptive_Statistics.csv")
cat("描述性统计表 (Table 1) 已生成。\n")

```

---
  
  ## 第三步：漏斗式筛选 (差异 -> 相关 -> 交集)
  
  ### 3.1 差异性筛选 (Kruskal-Wallis)
  
  ```r
# ==========================================================
# Phase 3.1: 差异性筛选 (Kruskal-Wallis Test)
# ==========================================================

# 初始化结果列表
kw_results <- data.frame()

for (var in chemical_vars) {
  # 构建公式
  f <- as.formula(paste(var, "~ Group"))
  # 运行非参数检验
  test <- kruskal.test(f, data = analysis_data)
  
  kw_results <- rbind(kw_results, data.frame(
    Computational_Name = var,
    KW_p_value = test$p.value
  ))
}

# 筛选 P < 0.05 的物质
diff_vars <- kw_results %>%
  filter(KW_p_value < 0.05) %>%
  pull(Computational_Name)

cat("差异性筛选完成，共找到", length(diff_vars), "种差异显著物质。\n")

```

### 3.2 相关性筛选 (Spearman) & 核心物质锁定

```r
# ==========================================================
# Phase 3.2 & 3.3: 相关性筛选与核心物质确证
# ==========================================================

# 计算 Spearman 相关性 (针对 Sourness_Score 连续变量)
cor_results <- analysis_data %>%
  select(Sourness_Score, all_of(chemical_vars)) %>%
  cor_test(Sourness_Score, method = "spearman") %>%
  select(var, cor, p) %>%
  rename(Computational_Name = var, Rho = cor, Cor_p_value = p)

# 筛选 P < 0.05 的物质
cor_vars_sig <- cor_results %>%
  filter(Cor_p_value < 0.05) %>%
  pull(Computational_Name)

# === Phase 3.3: 取交集 (Intersection) ===
core_substances <- intersect(diff_vars, cor_vars_sig)

# 生成最终核心物质结果表 (包含原始名称、相关系数、P值)
final_core_table <- cor_results %>%
  filter(Computational_Name %in% core_substances) %>%
  left_join(mapping_df, by = "Computational_Name") %>%
  left_join(kw_results, by = "Computational_Name") %>%
  arrange(desc(abs(Rho))) %>%
  select(Display_Name, Rho, Cor_p_value, KW_p_value, Computational_Name)

write_csv(final_core_table, "Table2_Core_Biomarkers.csv")
cat("核心物质筛选完成。共确证", length(core_substances), "种关键物质。\n")
print(final_core_table$Display_Name)

```

---
  
  ## 第四步：发表级可视化 (v5.0 重点)
  
  ### 4.1 箱线图 (带显著性标记与抖动点)
  
  ```r
# ==========================================================
# Phase 4.1: 核心物质差异箱线图 (Boxplot)
# ==========================================================
# 选择前 6-9 个最相关的核心物质进行绘图，避免图太多
top_features <- final_core_table$Computational_Name[1:min(9, length(core_substances))]

# 准备绘图数据 (长格式)
plot_data_box <- analysis_data %>%
  select(Group, all_of(top_features)) %>%
  pivot_longer(cols = -Group, names_to = "Computational_Name", values_to = "Value") %>%
  left_join(mapping_df, by = "Computational_Name")

# 绘制箱线图
p_box <- ggplot(plot_data_box, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) + # 抖动点展示真实分布
  facet_wrap(~Display_Name, scales = "free_y") +       # 分面展示
  scale_color_npg() +                                  # Nature 配色
  stat_compare_means(label = "p.signif",               # 自动添加显著性标记 (*, **)
                     method = "kruskal.test", 
                     label.x.npc = "center") +
  labs(x = "酸感分组", y = "标准化含量 (Z-score)", title = "核心酸感物质的分布差异") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0))

print(p_box)
ggsave("Figure1_Boxplots.png", p_box, width = 12, height = 10, dpi = 300)

```

### 4.2 相关性条形图 (保留 v3.3 风格)

```r
# ==========================================================
# Phase 4.2: 相关性条形图 (Correlation Bar Plot)
# ==========================================================

# 准备数据
plot_data_cor <- final_core_table %>%
  mutate(Type = ifelse(Rho > 0, "正相关 (Promoter)", "负相关 (Inhibitor)")) %>%
  mutate(Signif_Label = case_when(
    Cor_p_value < 0.001 ~ "***",
    Cor_p_value < 0.01 ~ "**",
    Cor_p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

p_cor <- ggplot(plot_data_cor, aes(x = reorder(Display_Name, Rho), y = Rho, fill = Type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  geom_text(aes(label = Signif_Label, 
                y = Rho + 0.05 * sign(Rho)), 
            vjust = 0.7) +
  scale_fill_manual(values = c("正相关 (Promoter)" = "#E64B35", # NPG Red
                               "负相关 (Inhibitor)" = "#4DBBD5")) + # NPG Blue
  labs(x = "核心物质", y = "Spearman 相关系数 (Rho)", title = "物质与酸感的相关性强度") +
  theme(legend.position = "top")

print(p_cor)
ggsave("Figure2_Correlation_Bar.png", p_cor, width = 8, height = 6, dpi = 300)

```

### 4.3 聚类热图 (带分组注释)

这里使用了您之前代码中的聚类参数 (`correlation`, `ward.D2`)，并加入了**分组注释条**，这是验证聚类效果的关键。

```r
# ==========================================================
# Phase 4.3: 验证热图 (Clustered Heatmap)
# ==========================================================

# 1. 准备矩阵数据
heatmap_matrix <- analysis_data %>%
  select(all_of(core_substances)) %>%
  as.data.frame()
rownames(heatmap_matrix) <- analysis_data$SampleID # 确保行名是样品ID

# 2. 准备列名显示 (映射回中文)
current_names <- colnames(heatmap_matrix)
display_names <- mapping_df$Display_Name[match(current_names, mapping_df$Computational_Name)]
colnames(heatmap_matrix) <- display_names

# 3. 准备注释条 (Annotation Row) - 关键步骤！
# 这将在热图左侧显示该样品原本属于哪个酸感组
annotation_row <- data.frame(Acid_Group = analysis_data$Group)
rownames(annotation_row) <- rownames(heatmap_matrix)

# 定义注释条颜色
ann_colors <- list(
  Acid_Group = c(Trace = "#91D1C2B2", Mild = "#4DBBD5B2", 
                 Moderate = "#E64B35B2", Strong = "#3C5488B2")
)

# 4. 绘制热图
pheatmap(
  t(heatmap_matrix),           # 转置：让物质在行，样品在列 (更符合生物信息学习惯)
  scale = "row",               # 按行归一化
  annotation_col = annotation_row, # 添加分组注释
  annotation_colors = ann_colors,
  clustering_distance_rows = "correlation", # 沿用 v3.3 参数
  clustering_distance_cols = "euclidean",   
  clustering_method = "ward.D2",            # 沿用 v3.3 参数
  show_colnames = FALSE,       # 样品名太多通常不显示
  fontsize_row = 8,
  main = "核心物质的层次聚类回验",
  filename = "Figure3_Heatmap_Validation.png",
  width = 10, height = 8
)

cat("分析全部完成。请查看生成的 CSV 表格和 PNG 图片。\n")

```