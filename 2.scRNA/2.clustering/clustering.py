import scanpy as sc
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import warnings
import sys
warnings.filterwarnings("ignore", category=UserWarning, module='scanpy')
warnings.filterwarnings("ignore", category=UserWarning, module='umap')
#warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


# 初始化解析器
parser = argparse.ArgumentParser(description="Seurat-like Processing Pipeline")

# 添加参数
parser.add_argument('-i', '--input', type=str, required=True, help="Input directory")
parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")
parser.add_argument('-n', '--name', type=str, required=True, help="Prefix name")


# 参数可以有Default值
parser.add_argument('--min_cells', type=int, default=3, help="Default: 3")
parser.add_argument('--min_features', type=int, default=200, help="Default: 200")
parser.add_argument('--nFeature_RNA_min', type=int, default=200, help="Default: 200")
parser.add_argument('--nFeature_RNA_max', type=int, default=2500, help="Default: 2500")
parser.add_argument('--mt_percent', type=float, default=5, help="Default: 5")
parser.add_argument('--variable_features', type=int, default=2000, help="Default: 2000")
parser.add_argument('--dim_prefer', type=int, default=30, help="Default: 30")
parser.add_argument('--res', type=float, default=1, help="Default: 1")

# 解析参数
args = parser.parse_args()

# 输出获取到的参数
print(f"输入路径: {args.input}")
print(f"输出路径: {args.output}")
print(f"最小细胞数: {args.min_cells}")
print(f"最小特征数: {args.min_features}")
print(f"nFeature_RNA_min: {args.nFeature_RNA_min}")
print(f"nFeature_RNA_max: {args.nFeature_RNA_max}")
print(f"线粒体百分比: {args.mt_percent}")
print(f"变量特征数: {args.variable_features}")
print(f"首选维度数: {args.dim_prefer}")
print(f"分辨率: {args.res}")
print(f"文件前缀: {args.name}")

# Command line arguments (replace with your actual input and output)
input_seurat = args.input  # 输入路径
output_directory = args.output  # 输出路径
sample_name = args.name

# Parameters with default values
min_cells = args.min_cells  # Default: 3
min_features = args.min_features  # Default: 200
nFeature_RNA_min = args.nFeature_RNA_min # Default: 200
nFeature_RNA_max = args.nFeature_RNA_max  # Default: 2500
mt_percent = args.mt_percent  # Default: 5
variable_features = args.variable_features  # Default: 2000
dim_prefer = args.dim_prefer # Default: 30
res = args.res

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)  # 创建目录，包括多级目录
    print(f"Output directory created: {output_directory}")

# Read input data
print("Reading input data...")
if input_seurat.endswith('.h5'):
    adata = sc.read_10x_h5(input_seurat)
elif input_seurat.endswith(('.txt.gz', '.csv.gz', '.tsv.gz')):
    adata = sc.read_text(input_seurat,var_names="gene_symbols")
    #adata.var_names = adata.var_names.astype(str)  # Ensure gene names are strings
elif os.path.isdir(input_seurat):
    adata = sc.read_10x_mtx(input_seurat,var_names="gene_symbols")
elif input_seurat.endswith('.h5ad'):
    adata = sc.read_h5ad(input_seurat,var_names="gene_symbols")
elif input_seurat.endswith('.rds'):
    print(f"Input is an RDS file, only run in seurat")
    sys.exit(1)
else:
    raise ValueError(f"Input {input_seurat} does not meet our requirements")

print(f"Input Seurat file: {input_seurat}")
print(f"Output directory: {output_directory}")

# Create Seurat object (Scanpy's AnnData object)
print(f"Creating AnnData object...\nMin cells: {min_cells}\nMin features: {min_features}")
sc.pp.filter_cells(adata,min_genes=min_features) #200
sc.pp.filter_genes(adata,min_cells=min_cells) #3

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt'],percent_top=None,log1p=False,inplace=True
)

adata = adata[adata.obs.n_genes_by_counts >= nFeature_RNA_min, :]  # Remove cells with fewer than min_features
adata = adata[adata.obs.n_genes_by_counts < nFeature_RNA_max, :]  # Filter cells based on max features
adata = adata[adata.obs.pct_counts_mt < mt_percent, :]  # Filter cells based on mitochondrial percentage

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes= variable_features)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

sc.tl.pca(adata, svd_solver="arpack")

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=dim_prefer)

sc.tl.umap(adata)
sc.tl.tsne(adata)

sc.tl.leiden(
    adata,
    resolution=res,
    random_state=0,
    #flavor="igraph",
    n_iterations=2,
    directed=False,
)

with PdfPages(f'{output_directory}/{sample_name}_Scanpy_clustering.pdf') as pdf:
    # UMAP 图像保存到 PDF
    sc.pl.umap(adata, color="leiden", show=False)
    pdf.savefig()  # 保存当前图像
    plt.close()

    # t-SNE 图像保存到 PDF
    sc.pl.tsne(adata, color="leiden", show=False)
    pdf.savefig()  # 保存当前图像
    plt.close()

adata.write(f"{output_directory}/{sample_name}_scanpy_clustering.h5ad")
print('Finished! You can check the result now')
