
import pandas as pd

# 定义文件路径
file_path = "c:/Users/Administrator/Desktop/PBCM/GSE118389_norm_data.txt"

# 读取文件，找到基因表达矩阵的起始行（忽略以 "!" 开头的元数据）
with open(file_path, "r") as f:
    lines = f.readlines()

# 找到数据起始行（第一个不以 "!" 开头的行）
data_start_idx = next(i for i, line in enumerate(lines) if not line.startswith("!"))

# 读取基因表达数据，并设置第一列为索引（基因名称）
df = pd.read_csv(file_path, sep="\t", skiprows=data_start_idx, index_col=0)

# 保存转换后的基因表达矩阵
df.to_csv("C:/Users/Administrator/Desktop/PBCM/GSE118389_expression_matrix.csv")

print("基因表达矩阵已转换并保存！")
'''
import pandas as pd

# 文件路径
file_path = "C:/Users/Administrator/Desktop/PBCM/GSE118390-GPL9052_series_matrix.txt"

# 读取文件内容
with open(file_path, "r") as f:
    lines = f.readlines()

# 查找基因表达数据的起始行
valid_data_start = None
for i, line in enumerate(lines):
    if not line.startswith("!") and "ID_REF" in line:  # 识别基因表达数据的表头
        valid_data_start = i
        break

# 如果找到基因表达数据
if valid_data_start:
    # 读取基因表达数据
    df_expression = pd.read_csv(file_path, sep="\t", skiprows=valid_data_start)

    # 确保第一列是基因 ID，其余列是样本表达值
    df_expression = df_expression.set_index(df_expression.columns[0])

    # 只保留数值型数据（去掉任何非数值行）
    df_expression = df_expression.apply(pd.to_numeric, errors='coerce').dropna()

    # 保存为 CSV 文件
    output_path = "C:/Users/Administrator/Desktop/PBCM/GSE118390_expression_matrix.csv"
    df_expression.to_csv(output_path)

    print(f"基因表达矩阵已成功提取，并保存到 {output_path}")

else:
    print("未找到基因表达矩阵，请检查文件格式！")
'''