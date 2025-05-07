import pandas as pd
import sys

# 读取 CSV 文件
df = pd.read_csv(sys.argv[1])

# 转置 DataFrame
df_transposed = df.T

# 保存为新文件
df_transposed.to_csv(sys.argv[2], header=False)
