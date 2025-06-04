# 安装和加载所需的包
install.packages("R6")
install.packages("jsonlite")
install.packages("vscDebugger")
if (!requireNamespace("LEAP", quietly = TRUE)) {
  install.packages("LEAP")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

library(LEAP)
library(igraph)

# 输入数据准备
# 替换为您的基因表达矩阵 CSV 文件路径
expression_file <- "C:/Users/Administrator/Desktop/单细胞2/ZINB_GRAN-main/Datasets/500_ChIP-seq_hESC/500_ChIP-seq_hESC-ExpressionData原.csv"
pseudotime_file <-"C:/Users/Administrator/Desktop/单细胞2/ZINB_GRAN-main/Datasets/500_ChIP-seq_hESC/PseudoTime.csv"  # 替换为伪时间序列文件路径

# 从 CSV 文件中读取基因表达矩阵
# 假设行为基因，列为细胞
expression_matrix <- as.matrix(read.csv(expression_file, row.names = 1))

# 从 CSV 文件中读取伪时间序列
# 假设伪时间文件只有一列，并与表达矩阵列顺序一致
pseudotime_order <- as.numeric(read.csv(pseudotime_file, header = FALSE)[, 1])

# 确保伪时间序列与表达矩阵的列数匹配
if (length(pseudotime_order) != ncol(expression_matrix)) {
  stop("伪时间序列长度与基因表达矩阵列数不一致！")
}

# 构建基因共表达网络
# 使用 LEAP 包计算基因间的时间滞后相关性矩阵
leap_results <- LEAP(expression_matrix, pseudotime_order)

# 提取相关性矩阵
correlation_matrix <- leap_results$correlations
print("相关性矩阵的前几行：")
print(head(correlation_matrix))

# 筛选显著调控关系（相关性阈值）
threshold <- 0.5  # 设置阈值
significant_edges <- which(correlation_matrix > threshold, arr.ind = TRUE)

# 输出显著调控关系
print("显著调控关系（基因对）：")
print(significant_edges)

# 可视化基因调控网络
# 使用 igraph 包创建网络图
network <- graph_from_adjacency_matrix(correlation_matrix > threshold, mode = "undirected")
plot(network, vertex.label = rownames(expression_matrix), 
     main = "Gene Co-expression Network")

# 保存结果
# 保存相关性矩阵
write.csv(correlation_matrix, "correlation_matrix.csv")

# 保存显著边（基因对）
write.table(significant_edges, "significant_edges.txt", row.names = FALSE, col.names = FALSE)
