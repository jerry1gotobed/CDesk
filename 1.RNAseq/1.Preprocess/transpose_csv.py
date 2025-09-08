import pandas as pd
import sys

df = pd.read_csv(sys.argv[1])
df_transposed = df.T
df_transposed.to_csv(sys.argv[2], header=False)
