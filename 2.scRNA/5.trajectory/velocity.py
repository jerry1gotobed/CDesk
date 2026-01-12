import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy as sc
import igraph
import loompy as lmp
import anndata
import scvelo as scv
import warnings
import os
import shutil
import subprocess
import sys
import json

warnings.filterwarnings('ignore')
# Set parameters for plots, including size, color, etc.
#pl.rcParams["figure.figsize"] = (10,10)
Colorss=["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF"]
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run scRNA velocity")
    
    parser.add_argument("--cellranger_output", type=str, help="Input cellranger output directory")
    parser.add_argument("--gtf",default='',type=str, help="Specify the gtf file")

    parser.add_argument("--rds_file", type=str, help="Input rds file")

    #parser.add_argument("--loom_file", type=str, help="Input velocity loom file")
    #parser.add_argument("--cellID_file", type=str, help="Input cell ID csv file")
    #parser.add_argument("--cellEmbedding_file", type=str, help="Input cell embedding csv file")
    parser.add_argument("-t", type=int,default=10,help="Number of threads")
    parser.add_argument("--genes", type=str,default='',help="Interested genes file")
    parser.add_argument("-o", type=str, help="Output directory")
    parser.add_argument("--meta", type=str, help="Meta column used")
    parser.add_argument("--width", type=int,default=10 ,help="Plot width, default: 10")
    parser.add_argument("--height", type=int, default=10,help="Plot height, default: 10")
    return parser.parse_args()

def main():
    args = parse_arguments()

    config = os.environ.get('CDesk_config')
    with open(config, "r") as f:
        config = json.load(f)

    output_dir = args.o
    threads = int(args.t)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir,'input'), exist_ok=True)
    input_data = os.path.join(output_dir,'input')
    meta = args.meta
    plot_width = int(args.width)
    plot_height = int(args.height)
    pl.rcParams["figure.figsize"] = (plot_width,plot_height)
    
    # Velocity analysis

    # Check if there is loom result
    cellranger_result = args.cellranger_output
    if os.path.exists(os.path.join(cellranger_result,'velocyto')):
        temp = os.path.join(cellranger_result,'velocyto')
        files = os.listdir(temp)
        loom_file = [f for f in files if f.endswith('.loom')][0]
        loom_file = os.path.join(temp,loom_file)
        shutil.copy2(loom_file, input_data)
    else:
        if args.gtf == '':
            print('Please specify the gtf file')
            sys.exit(1)
        gtf = args.gtf
        cmd = ['velocyto','run10x',cellranger_result,gtf,'-@',str(threads),'--verbose','-v']
        log_path = os.path.join(output_dir,'velocyto.log')
        print('Run velocyto ...')
        with open(log_path, 'w') as logfile:
            subprocess.run(
                cmd,
                stdout=logfile,          
                stderr=subprocess.STDOUT,
                check=True
            )
        temp = os.path.join(cellranger_result,'velocyto')
        files = os.listdir(temp)
        if len([f for f in files if f.endswith('.loom')]) == 0:
            print('No loom file found in the velocyto directory of cellranger result directory')
            sys.exit(1)
        loom_file = [f for f in files if f.endswith('.loom')][0]
        loom_file = os.path.join(temp,loom_file)
        shutil.copy2(loom_file, input_data)

    # Grab the seurat meta
    seurat_file = args.rds_file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cmd = ['Rscript',os.path.join(current_dir,'velocity.R'),seurat_file,input_data,meta]
    subprocess.run(cmd, check=True)

    # scVelo analysis
    # Read the loom file
    ldata = sc.read(loom_file,cache=False)
    ldata.var_names_make_unique()
    # Mofify the cell barcodes to keep the same with seurat file
    barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
    barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
    ldata.obs.index = barcodes
    ldata.var_names_make_unique()
    # Read the barcode, UMAP and Seurat Grouping result
    cellID_obs = pd.read_csv(os.path.join(input_data,'cellID_obs.csv'))
    umap_cord = pd.read_csv(os.path.join(input_data,'cell_embeddings.csv'))
    cell_clusters = pd.read_csv(os.path.join(input_data,'clusters.csv'))
    # Choose the cells we want according to the barcode
    filtered_ldata = ldata[cellID_obs['x']].copy()
    # Add UMAP information
    ldata_index = pd.DataFrame(filtered_ldata.obs.index)
    ldata_index = ldata_index.rename(columns={ldata_index.columns[0]: 'CellID'})
    umap_cord = umap_cord.rename(columns={'Unnamed: 0': 'CellID'})

    umap_ordered = ldata_index.merge(umap_cord, on="CellID")
    umap_ordered = umap_ordered.iloc[:, 1:]
    filtered_ldata.obsm['X_umap'] = umap_ordered.values
    # Add grouping information
    cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'CellID'})
    cell_clusters = ldata_index.merge(cell_clusters, on="CellID")
    cell_clusters = cell_clusters.iloc[:, 1:]
    filtered_ldata.obs['seurat_clusters'] = cell_clusters[meta].values

    adata = filtered_ldata
    adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

    # Calculate the proportion of spliced and unspliced
    scv.pl.proportions(adata, groupby='seurat_clusters',save=f'{output_dir}/splice_proportions.pdf')

    # Velocity analysis
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)

    # Velocity analysis result visualization
    scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters',scale=0.25,save=f'{output_dir}/velocity.pdf')
    scv.pl.velocity_embedding_stream(adata, basis='umap',color='seurat_clusters',save=f'{output_dir}/velocity_with_embeddings.pdf')
    
    # Identify important genes
    scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
    ranked_genes = adata.uns['rank_velocity_genes']['names']
    df = pd.DataFrame(ranked_genes)
    df.to_csv(f'{output_dir}/velocity_key_genes.csv', index=False)

    # Interested genes
    if args.genes != '':
        with open(args.genes, 'r') as file:
            genes_list = [line.strip() for line in file if line.strip()]
        try:
            scv.pl.velocity(adata, pd.Series(genes_list),ncols=2,save=f'{output_dir}/interested_genes.pdf')
        except Exception as e:
            print(f"An error occurred: {e}")

    # Pseudotime calculate and visualization
    scv.tl.velocity_pseudotime(adata)  # 0 - 1
    scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',save=f'{output_dir}/velocity_pseudotime.pdf') 

    # PAGA analysos
    scv.tl.paga(adata, groups='seurat_clusters')
    scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
                min_edge_width=2, node_size_scale=1.5,save=f'{output_dir}/PAGA_velocity.pdf')



if __name__ == "__main__":

    main()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    figures_dir = os.path.join(current_dir, 'figures')

    if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
        # Remove figures folder
        shutil.rmtree(figures_dir)
    print('Done, you can see the results now.')

# python velocity.py --cellranger_output /mnt/linzejie/CDesk/test/result/2.scRNA/1.preprocess/cellranger/pbmc --gtf /mnt/linzejie/scRNA/count/refdata-gex-GRCh38-2020-A/genes/genes.gtf --rds_file /mnt/linzejie/CDesk/test/data/2.scRNA/7.trajectory/RNAvelocity/velocity.rds --genes SP100,ACPP,IFNG,RGS3 -o /mnt/linzejie/CDesk/test/result/2.scRNA/7.trajectory/RNAvelocity
