# 加载所需 R 包
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ggpubr)

# 设置路径
input_file <- "C:/Users/24071/Desktop/sepsis_platelet/210and252DEGS/common_upregulated_no_dot_genes_ascending_order.xlsx"
output_excel <- "C:/Users/24071/Desktop/sepsis_platelet/210and252DEGS/GO_result_by_type.xlsx"
go_plot_output <- "C:/Users/24071/Desktop/sepsis_platelet/210and252DEGS/GO_three_types_combined_plot.pdf"
kegg_excel <- "C:/Users/24071/Desktop/sepsis_platelet/210and252DEGS/KEGG_result.xlsx"
kegg_plot_pdf <- "C:/Users/24071/Desktop/sepsis_platelet/210and252DEGS/KEGG_dotplot_barplot.pdf"

# 读取基因名
gene_data <- read_excel(input_file)
gene_list <- na.omit(gene_data$gene)

# SYMBOL 转 ENTREZID
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# ---------------- GO分析 ----------------
ego_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", 
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
ego_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", 
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
ego_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", 
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

# 写入Excel，不同sheet
write.xlsx(list(
  Biological_Process = as.data.frame(ego_bp),
  Cellular_Component = as.data.frame(ego_cc),
  Molecular_Function = as.data.frame(ego_mf)
), file = output_excel)

# ---------------- GO可视化（合并三类） ----------------
# 添加类别标签
ego_bp@result$Category <- "BP"
ego_cc@result$Category <- "CC"
ego_mf@result$Category <- "MF"

# 合并前三十的富集结果
go_all <- bind_rows(
  ego_bp@result %>% top_n(-10, p.adjust),
  ego_cc@result %>% top_n(-10, p.adjust),
  ego_mf@result %>% top_n(-10, p.adjust)
)

# 画气泡图
go_plot <- ggplot(go_all, aes(x = Category, y = reorder(Description, p.adjust), 
                              size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "Top GO Terms (BP, CC, MF)", x = "GO Category", y = "Term") +
  theme_minimal(base_size = 12)

ggsave(go_plot_output, plot = go_plot, width = 10, height = 8)

# ---------------- KEGG分析 ----------------
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", 
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)

# 保存KEGG分析结果
write.xlsx(as.data.frame(ekegg), file = kegg_excel)

# ---------------- KEGG可视化 ----------------
# 气泡图
dot_kegg <- dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Dotplot") +
  theme(axis.text.y = element_text(size = 10))

# 柱状图
bar_kegg <- barplot(ekegg, showCategory = 20, title = "KEGG Pathway Barplot") +
  theme(axis.text.y = element_text(size = 10))

# 合并图像保存
pdf(kegg_plot_pdf, width = 12, height = 8)
print(dot_kegg)
print(bar_kegg)
dev.off()