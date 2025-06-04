##preprocessing single-cell expression data, filtering low-quality genes and cells, transforming data to the LTMG model, meanwhile sorting the cells by the pseudo-time 
import time
import argparse
import numpy as np
import pandas as pd
import pickle
import os.path
import scipy.sparse as sp
import scipy.io
'''
parser = argparse.ArgumentParser(description='Arguments for Tools')
parser.add_argument('--expressionFile', type=str, default='Use_expression.csv',
                    help='expression File in csv')
parser.add_argument('--processedFile', type=str, default='processed_expression.csv',
                    help='processed expression File in csv')
parser.add_argument('--cellLabelFile', type = str, default = None,
                    help = 'Cell labels file in csv')
parser.add_argument('--LTMGDir', type=str, default=None,
                    help='directory of LTMG model produced.')
parser.add_argument('--orderedCellsFile', type=str, default=None,
                    help='Directory of odered cell expression matrix by pseudotime.')
parser.add_argument('--filetype', type=str, default='CSV',
                    help='select input filetype, 10X or CSV')
parser.add_argument('--delim', type=str, default='comma',
                    help='File delim type, comma or space: default(comma)')
parser.add_argument('--transform', type=str, default='log',
                    help='Whether transform')
parser.add_argument('--cellRatio', type=float, default=0.99,
                    help='cell ratio')
parser.add_argument('--geneRatio', type=float, default=0.99,
                    help='gene ratio')
parser.add_argument('--geneCriteria', type=str, default='variance',
                    help='gene Criteria')
parser.add_argument('--transpose', action='store_true', default=False,
                    help='whether transpose or not')
parser.add_argument('--isUMI', action='store_true', default=False,
                    help='whether UMI or not')


args = parser.parse_args()

expressionFile = "C:/Users/Administrator/Desktop/单细胞2/ZINB_GRAN-main/Datasets500_ChIP-seq_mHSC-L\500_ChIP-seq_mHSC-L-ExpressionData.csv"  # 指定本地的输入文件路径
processedFile = '/mnt/data/processed_expression.csv'               # 指定预处理后文件的保存路径
cellLabelFile = None  # 如果有标签文件可以指定路径；否则保持为 None
LTMGDir = None  # 如果需要存储LTMG模型的输出目录，可以指定路径；否则保持为 None
orderedCellsFile = None  # 如果需要存储按伪时间排序后的细胞表达数据，可以指定路径；否则保持为 None
filetype = 'CSV'
delim = 'comma'
transform = 'log'
cellRatio = 0.99
geneRatio = 0.99
geneCriteria = 'variance'
transpose = False
isUMI = False

def preprocessingCSV(expressionFile, processedFile, delim = 'comma', transform = 'log', cellRatio=0.9, geneRatio=0.9, geneCriteria = 'variance', transpose = False):
    
    #preprocessing CSV files:
    #transform='log' or None
    
    if not os.path.exists(expressionFile):
        print('Dataset ' + expressionFile + ' not exists!')

    print("Input scRNA data in CSV format, start to reading...")

    if delim == 'space':
        df = pd.read_csv(expressionFile, index_col = 0, delim_whitespace = True)
    elif delim == 'comma':
        df = pd.read_csv(expressionFile, index_col = 0)

    print('Expression File loading complete, start to filter low-quality genes and cells')

    if transpose ==True:
        df = df.T

    df1 = df[df.astype('bool').mean(axis=1) >= (1-geneRatio)]
    print('After preprocessing, {} genes remaining'.format(df1.shape[0]))
    criteriaGene = df1.astype('bool').mean(axis=0) >= (1-cellRatio)
    df2 = df1[df1.columns[criteriaGene]]
    print('After preprocessing, {} cells have {} nonzero'.format(
        df2.shape[1], geneRatio))

    if transform == 'log':
        df2 = df2.transform(lambda x: np.log(x + 1))
    print(df2.shape)
    df.to_csv(processedFile)

def preprocessingH5(expressionFile, sc_gene_list, processedFile, cellRatio=0.9, geneRatio=0.9, geneCriteria = 'variance', transpose = False):
    if not os.path.exists(expressionFile):
        print('Dataset ' + expressionFile + ' not exists!')
    
    print("Input scRNA data in H5 format, start to reading...")
    f = pd.HDFStore(expressionFile)
    df = f['RPKMs']
    gene_name_dict = get_gene_list(sc_gene_list)
    genes = df.columns
    gene_symbols = []
    for gene in genes:
        gene_symbols.append(gene_name_dict[str(gene)])

    df.columns = gene_symbols

    print('Expression File loading complete, start to filter low-quality genes and cells')

    if transpose ==True:
        df = df.T

    df1 = df[df.astype('bool').mean(axis=1) >= (1-geneRatio)]
    print('After preprocessing, {} genes remaining'.format(df1.shape[0]))
    criteriaGene = df1.astype('bool').mean(axis=0) >= (1-cellRatio)
    df2 = df1[df1.columns[criteriaGene]]
    print('After preprocessing, {} cells have {} nonzero'.format(
        df2.shape[1], geneRatio))

    if transform == 'log':
        df2 = df2.transform(lambda x: np.log(x + 1))
    print(df2.shape)
    df.to_csv(processedFile)

def get_gene_list(file_name):
    import re
    h={}
    s = open(file_name,'r') #gene symbol ID list of sc RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(2)]=search_result.group(1).lower() # h [gene ID] = gene symbol
    s.close()
    return h

def transfOdiM2ExprsM(sparseMatrix, processedFile):
    df = pd.read_csv(sparseMatrix, header=None,
                     skiprows=1, delim_whitespace=True)
    counts = len(df) - 1
    proc_matrix = pd.read_csv(processedFile, header = 0, index_col = 0)
    for row in df.itertuples():
        # For the first row, it contains the number of genes and cells. Init the whole matrix
        if row[2] == counts:
            matrix = np.zeros((row[0], row[1]))
        else:
            matrix[row[2]-1][row[1]-1] = proc_matrix[row[2]-1][row[1]-1]

def computeCorr(expressionFile, corrMethod = 'pearson', threshold = 0.4):
    exprs = pd.read_csv(expressionFile, header = 0, index_col = 0, sep = "\t")
    exprsT = pd.DataFrame(exprs.T, columns = exprs.index, index = exprs.columns)
    corr = exprsT.corr(method= corrMethod)
    corr.values[corr.values < threshold] = 0.0
    return corr 

if __name__ == "__main__":
    start_time = time.time()

    # preprocessing
    print('Step1: Start filter and generating CSV')
    #if args.filetype == '10X':
        #expressionFilename = args.LTMGDir+args.datasetName+'/'+args.expressionFile
        # data = preprocessing10X(args.datasetDir, args.datasetName, args.LTMGDir+args.datasetName+'/'+args.expressionFile, args.transform, args.cellRatio, args.geneRatio, args.geneCriteria, args.geneSelectnum)
        #preprocessing10X(args.expressionFile, processedFile, args.transform,args.cellRatio, args.geneRatio, args.geneCriteria, args.geneSelectnum)
    if args.filetype == 'CSV':
        preprocessingCSV(args.expressionFile, args.processedFile, args.delim, args.transform,
                         args.cellRatio, args.geneRatio, args.geneCriteria, args.transpose)

    # start LTMG transforming
    from util.LTMG_Monocle_R import *
    print('Step2: Start infer LTMG and Pseudotime from Expression matrix')
    # run LTMG in R
    runLTMG(args.expressionFile, args.processedFile, args.cellLabelFile, args.LTMGDir, args.orderedCellsFile, args.isUMI)

    print("Preprocessing Done. Total Running Time: %s seconds" %
          (time.time() - start_time))
'''

'''
import os
import numpy as np
import pandas as pd
import time

# 直接在代码中定义参数，而不是使用 argparse
expressionFile = "C:/Users/Administrator/Desktop/PBCM/GSE118389_expression_matrix.csv"
processedFile = "C:/Users/Administrator/Desktop/PBCM/GSE118389_expression_matrix1.csv"
cellLabelFile = None  # 如果有标签文件可以指定路径，否则保持为 None
LTMGDir = None  # LTMG模型输出目录，如果不需要生成则保持为 None
orderedCellsFile = None  # 存储按伪时间排序的文件路径，保持为 None 时不生成
filetype = 'CSV'
delim = 'comma'
transform = 'log'
cellRatio = 0.9
geneRatio = 0.9
geneCriteria = 'variance'
transpose = False
isUMI = False

def preprocessingCSV(expressionFile, processedFile, delim='comma', transform='log', cellRatio=0.9, geneRatio=0.9, geneCriteria='variance', transpose=False):
    
    #预处理CSV文件：过滤低质量基因和细胞，选择性进行log变换
    
    if not os.path.exists(expressionFile):
        print('Dataset ' + expressionFile + ' does not exist!')
        return

    print("Input scRNA data in CSV format, start to reading...")

    # 根据指定的分隔符读取数据
    if delim == 'space':
        df = pd.read_csv(expressionFile, index_col=0, delim_whitespace=True)
    elif delim == 'comma':
        df = pd.read_csv(expressionFile, index_col=0)

    print('Expression File loading complete, start to filter low-quality genes and cells')

    # 如果需要，转置数据矩阵
    if transpose:
        df = df.T

    # 过滤低质量基因和细胞
    df1 = df[df.astype(bool).mean(axis=1) >= (1 - geneRatio)]
    print('After preprocessing, {} genes remaining'.format(df1.shape[0]))
    criteriaGene = df1.astype(bool).mean(axis=0) >= (1 - cellRatio)
    df2 = df1[df1.columns[criteriaGene]]
    print(f'After preprocessing, {df2.shape[1]} cells remain after filtering with {geneRatio} nonzero gene ratio')

    # 是否对数据进行log变换
    if transform == 'log':
        df2 = df2.transform(lambda x: np.log(x + 1))
    print(df2.shape)

    # 保存预处理后的数据
    df2.to_csv(processedFile)
    print(f'Processed data saved to {processedFile}')

    # 计算每个基因的非零表达比例
    nonzero_gene_ratio = (df2.iloc[:, 1:] > 0).mean(axis=1)

# 计算有多少基因的非零表达比例 < 10%（意味着这些基因可能是低表达基因）
    low_expr_genes = (nonzero_gene_ratio < 0.1).sum()
    total_genes = df2.shape[0]

    print(f"低表达基因数: {low_expr_genes} / {total_genes} ({(low_expr_genes/total_genes)*100:.2f}%)")

# 仅在文件类型为 CSV 时调用 CSV 预处理函数
if filetype == 'CSV':
    start_time = time.time()
    print('Step 1: Start filtering and generating processed CSV')
    preprocessingCSV(expressionFile, processedFile, delim, transform, cellRatio, geneRatio, geneCriteria, transpose)
    print("Preprocessing Done. Total Running Time: %s seconds" % (time.time() - start_time))

# 这里的LTMG步骤涉及 R 脚本的调用，根据需要可以添加适当的导入和调用
# 示例:
# from util.LTMG_Monocle_R import *
# if LTMGDir:
#     print('Step 2: Start inferring LTMG and Pseudotime from Expression matrix')
#     runLTMG(expressionFile, processedFile, cellLabelFile, LTMGDir, orderedCellsFile, isUMI)

'''
'''
import pandas as pd

# 1. 读取 TF 数据文件（第一列是 TF）
#file_path ="C:/Users/Administrator/Desktop/PBCM/GARNet_results/CD8 T-TF-Target.csv"
#tf_data = pd.read_csv(file_path)
#CD14
tf_list=['GATA2', 'ATF3', 'MYB', 'JUND', 'PAX5', 'BACH2', 'BATF', 'GATA3', 'MYC','STAT1', 'CEBPB', 'SPI1', 'RPS4Y2', 'FCER1A', 'BCL11A']
#CD8 T
#tf_list=['BACH2', 'NFE2', 'SPI1', 'MAX', 'GATA2', 'MEF2C', 'GATA3', 'ATF3']
# 2. 提取 TF 名称并去重（转换为大写）
#tf_column = tf_data.iloc[:, 0] 
#tf_list = set(tf_column.astype(str).str.strip().str.upper())  

# 3. 读取转录组表达数据
expression_matrix = pd.read_csv("C:/Users/Administrator/Desktop/PBCM/CD14Monocyte_gene_expression.csv", index_col=0)

# 4. 确保表达矩阵的索引统一为大写并去除空格
expression_matrix.index = expression_matrix.index.astype(str).str.strip().str.upper()

# 5. 计算基因表达变异性（标准差）
gene_variability = expression_matrix.std(axis=1)

# 6. 选取非 TF 基因
non_tf_genes = expression_matrix.drop(tf_list, errors='ignore')  # 移除 TF 基因
non_tf_variability = non_tf_genes.std(axis=1)

# 7. 选取前 500 个高变异非 TF 基因
selected_non_tfs = non_tf_variability.sort_values(ascending=False).head(500)

# 8. 合并最终用于 GRN 推断的数据
selected_genes = list(tf_list) + list(selected_non_tfs.index)

# 9. 提取最终的表达矩阵
selected_expression_matrix = expression_matrix.loc[selected_genes]

# 10. 确保基因名（索引）都是大写
selected_expression_matrix.index = selected_expression_matrix.index.str.upper()

# 11. 保存结果（确保基因名称为大写）
output_path = "C:/Users/Administrator/Desktop/PBCM/CD8 T_gene_expression筛选.csv"
selected_expression_matrix.to_csv(output_path)

# 输出信息
print(f"选取了 {len(tf_list)} 个 TF 基因和 {len(selected_non_tfs)} 个高变异非 TF 基因，总共 {len(selected_genes)} 个基因用于 GRN 推断。")
print(f"数据已保存至: {output_path}")
'''



import pandas as pd

# 定义 TF 基因列表，并去除空格
#CD14
tf_list1 = ['GATA2', 'ATF3', 'MYB', 'JUND', 'PAX5', 'BACH2', 'BATF', 'GATA3', 'MYC', 'STAT1', 'CEBPB', 'SPI1', 'RPS4Y2', 'FCER1A', 'BCL11A']
#B
tf_list=['EGR1', 'EGR2', 'TBX21', 'MEF2C', 'GATA2', 'BACH2', 'STAT1', 'MYC', 'NFE2','PAX5', 'RUNX3', 'FCER1A']
tf_list = [tf.strip().upper() for tf in tf_list]  # 去除空格并转换为大写

# 读取转录组表达数据
expression_matrix = pd.read_csv("C:/Users/Administrator/Desktop/PBCM/B - 副本.csv", index_col=0)

# 确保基因名索引为大写、去除空格，防止匹配错误
expression_matrix.index = expression_matrix.index.astype(str).str.strip().str.upper()

# **检查 TF 基因是否存在**
existing_tf_genes = list(set(tf_list) & set(expression_matrix.index))
missing_tf_genes = list(set(tf_list) - set(expression_matrix.index))

print(f"找到 {len(existing_tf_genes)} 个 TF 基因: {existing_tf_genes}")
if missing_tf_genes:
    print(f"⚠️ 缺失 {len(missing_tf_genes)} 个 TF 基因: {missing_tf_genes}")

# 计算基因表达变异性（标准差）
gene_variability = expression_matrix.std(axis=1)

# **筛选出非 TF 基因**
non_tf_genes = expression_matrix.drop(existing_tf_genes, errors='ignore')  # 确保不会因缺失基因报错
non_tf_variability = non_tf_genes.std(axis=1)

# **选取前 500 个高变异非 TF 基因**
selected_non_tfs = non_tf_variability.sort_values(ascending=False).head(500)

# **合并最终用于 GRN 推断的数据**
selected_genes = existing_tf_genes + list(selected_non_tfs.index)

# **确保基因名索引匹配**
existing_selected_genes = list(set(selected_genes) & set(expression_matrix.index))
missing_selected_genes = list(set(selected_genes) - set(expression_matrix.index))

if missing_selected_genes:
    print(f"⚠️ 部分选定的基因在数据集中丢失: {missing_selected_genes}")

# **提取最终的表达矩阵**
selected_expression_matrix = expression_matrix.loc[existing_selected_genes]

# **保存结果**
output_path = "C:/Users/Administrator/Desktop/PBCM/B 筛选.csv"
selected_expression_matrix.to_csv(output_path)

# **输出信息**
print(f"✅ 选取了 {len(existing_tf_genes)} 个 TF 基因和 {len(selected_non_tfs)} 个高变异非 TF 基因，总共 {len(existing_selected_genes)} 个基因用于 GRN 推断。")
print(f"📂 数据已保存至: {output_path}")
