#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
此脚本使用pySCENIC进行单细胞RNA-seq数据的基因调控网络分析，可通过命令行直接调用。
支持SCENIC工作流程，包括共表达网络推断、顺式作用元件富集分析和稳定状态评分。

作者: yutiancheng
日期: 2025
"""

import os
import sys
import argparse
import subprocess
import tempfile
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Union, Optional
from multiprocessing import cpu_count
import requests
import json

def filter_genes(expr_matrix: pd.DataFrame, min_number_cells: int = 10) -> pd.DataFrame:
    """
    过滤基因表达矩阵，去除在少于指定数量细胞中表达的基因
    
    参数:
        expr_matrix (pd.DataFrame): 基因表达矩阵，行为细胞，列为基因
        min_number_cells (int): 基因需要在至少多少个细胞中表达才会被保留
    
    返回:
        pd.DataFrame: 过滤后的基因表达矩阵
    """
    # 计算每个基因在多少个细胞中表达
    n_cells = (expr_matrix > 0).sum(axis=0)
    
    # 选择在足够多细胞中表达的基因
    genes_to_keep = n_cells[n_cells >= min_number_cells].index
    
    # 返回过滤后的矩阵
    return expr_matrix[genes_to_keep]

def download_motif_annotations(output_dir: str, organism: str = "hg19") -> str:
    """
    下载motif注释文件
    
    参数:
        output_dir (str): 输出目录
        organism (str): 生物体标识符，如mm10(小鼠)或hg19(人类)
        
    返回:
        str: motif注释文件路径
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 选择适当的URL
    if organism.startswith("mm"):  # 小鼠
        url = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
        output_file = os.path.join(output_dir, f"motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
    else:  # 默认使用人类hg19
        url = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        output_file = os.path.join(output_dir, f"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
    
    print(f"正在从{url}下载motif注释文件...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"motif注释文件已下载到：{output_file}")
        return output_file
    except Exception as e:
        print(f"下载失败：{str(e)}")
        return None

class SCENICAnalyzer:
    """基于pySCENIC命令行工具的基因调控网络分析类"""
    
    def __init__(self, 
                 input_file: str = None, 
                 output_dir: str = None,
                 tf_list_file: str = None,
                 db_folder: str = None,
                 num_workers: int = None,
                 normalize: bool = False,
                 log: bool = False,
                 scale: bool = False,
                 config: str = None,
                 auc_threshold: float = 0.05,
                 nes_threshold: float = 3.0,
                 species: str = "hg19"):
        """
        初始化SCENIC分析器
        
        参数:
            input_file (str): 输入文件路径，支持h5ad或csv/tsv格式
            output_dir (str): 输出目录
            tf_list_file (str): 转录因子列表文件
            db_folder (str): 包含RcisTarget数据库的文件夹
            num_workers (int): 并行计算的工作进程数
            auc_threshold (float): AUCell阈值
            nes_threshold (float): 富集得分阈值
            species (str): 物种，例如"hg19"(人类)或"mm10"(小鼠)
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.tf_list_file = tf_list_file
        self.db_folder = db_folder
        self.auc_threshold = auc_threshold
        self.nes_threshold = nes_threshold
        self.anndata = None
        self.tf_list = []
        self.num_workers = num_workers if num_workers else cpu_count()
        self.log = log
        self.normalize = normalize
        self.scale = scale
        self.species = species
        with open(config, "r") as f:
            config = json.load(f)
        self.config = config

        # 创建输出目录
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # 加载转录因子列表
        if tf_list_file:
            self._load_tf_list()
        
    def _load_tf_list(self):
        """加载转录因子列表"""
        if not os.path.exists(self.tf_list_file):
            print(f"警告：指定的转录因子文件{self.tf_list_file}不存在")
            return
        
        with open(self.tf_list_file, 'r') as f:
            self.tf_list = [line.strip() for line in f if line.strip()]
        print(f"加载了{len(self.tf_list)}个转录因子")
    
    def load_data(self):
        """加载输入数据"""
        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f"输入文件{self.input_file}不存在")
        
        # 根据文件扩展名确定文件类型
        if self.input_file.endswith('.h5ad'):
            self.anndata = sc.read_h5ad(self.input_file)
        elif self.input_file.endswith('.csv'):
            self.anndata = sc.read_csv(self.input_file)
        elif self.input_file.endswith('.tsv') or self.input_file.endswith('.txt'):
            self.anndata = sc.read_text(self.input_file)
        else:
            raise ValueError("不支持的文件格式，请使用h5ad/csv/tsv格式文件")
        
        print(f"成功加载数据: {self.anndata.shape[0]}个细胞, {self.anndata.shape[1]}个基因")
    
    def preprocess(self, normalize: bool = False, log_transform: bool = False, scale: bool = False):
        """数据预处理"""
        if self.anndata is None:
            raise ValueError("请先加载数据")
        if self.normalize:
            sc.pp.normalize_total(self.anndata)
            print("完成总量归一化")
        if self.log:
            sc.pp.log1p(self.anndata)
            print("完成log1p转换")
        if self.scale:
            sc.pp.scale(self.anndata)
            print("完成数据缩放")
        # 保存处理后的数据
        self.processed_file = os.path.join(self.output_dir, "processed_data.loom")
        self._export_to_loom(self.processed_file)
        print(f"预处理后的数据已保存至: {self.processed_file}")
    
    def _export_to_loom(self, loom_file: str):
        """
        将AnnData对象导出为loom文件
        
        参数:
            loom_file (str): 输出loom文件路径
        """
        try:
            import loompy as lp
        except ImportError:
            print("正在安装loompy包...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", "loompy"])
            import loompy as lp
        
        # 导出为loom文件
        if isinstance(self.anndata.X, np.ndarray):
            matrix = self.anndata.X
        else:
            matrix = self.anndata.X.toarray()
        
        # loompy期望行是基因，列是细胞，所以需要转置
        matrix = matrix.T 
        
        # 确保基因名称是pySCENIC期望的格式
        # 从anndata对象中获取基因名
        gene_names = self.anndata.var.index.tolist()
        
        # 检查是否有合适的gene_symbol列，如果有则优先使用
        if 'gene_symbol' in self.anndata.var.columns:
            print("使用gene_symbol列作为基因标识符")
            gene_names = self.anndata.var['gene_symbol'].tolist()
        elif 'symbol' in self.anndata.var.columns:
            print("使用symbol列作为基因标识符")
            gene_names = self.anndata.var['symbol'].tolist()
        elif 'gene_name' in self.anndata.var.columns:
            print("使用gene_name列作为基因标识符")
            gene_names = self.anndata.var['gene_name'].tolist()
        
        # 确保所有基因名都是字符串类型
        gene_names = [str(name) for name in gene_names]
        
        print(f"基因名示例: {gene_names[:5]}")
        
        # 创建正确的行属性（基因）和列属性（细胞）
        row_attrs = {"Gene": np.array(gene_names)}
        col_attrs = {"CellID": np.array(self.anndata.obs.index)}
        
        # 使用with语句创建loom文件
        lp.create(loom_file, matrix, row_attrs, col_attrs)
            
        print(f"数据已导出为loom格式: {loom_file}")
    
    def run_scenic(self):
        """运行SCENIC分析流程，通过调用pyscenic命令行工具"""
        if not os.path.exists(self.processed_file):
            raise ValueError("请先运行预处理步骤")
        
        print("步骤1: 推断转录因子与靶基因之间的关联...")
        
        # 检查GRNBoost2输出文件
        grn_output = os.path.join(self.output_dir, "adjacencies.csv")
        if os.path.exists(grn_output) and os.path.getsize(grn_output) > 0:
            # 验证文件内容
            try:
                test_df = pd.read_csv(grn_output, nrows=5)
                if len(test_df) > 0 and 'TF' in test_df.columns and 'target' in test_df.columns:
                    print(f"发现有效的共表达网络结果文件: {grn_output}，跳过GRNBoost2步骤")
                else:
                    print(f"共表达网络结果文件格式有误，将重新生成")
                    self._run_grn_step(grn_output)
            except Exception as e:
                print(f"读取共表达网络结果文件时出错: {str(e)}，将重新生成")
                self._run_grn_step(grn_output)
        else:
            self._run_grn_step(grn_output)
        
        # 步骤2: 运行cistrome-AUC分析
        print("\n步骤2: 进行顺式调控元件富集分析...")
        
        # 检查motif富集结果文件
        motif_output = os.path.join(self.output_dir, "motifs.csv")
        if os.path.exists(motif_output) and os.path.getsize(motif_output) > 0:
            # 验证文件内容
            try:
                test_df = pd.read_csv(motif_output, nrows=5)
                if len(test_df) > 0:
                    print(f"发现有效的顺式调控元件富集结果文件: {motif_output}，跳过ctx步骤")
                else:
                    print(f"顺式调控元件富集结果文件内容为空，将重新生成")
                    self._run_ctx_step(grn_output, motif_output)
            except Exception as e:
                print(f"读取顺式调控元件富集结果文件时出错: {str(e)}，将重新生成")
                self._run_ctx_step(grn_output, motif_output)
        else:
            self._run_ctx_step(grn_output, motif_output)
        
        # 步骤3: 进行AUCell分析
        print("\n步骤3: 计算调控子活性...")
        
        # 检查AUCell结果文件
        auc_output = os.path.join(self.output_dir, "scenic_result.loom")
        if os.path.exists(auc_output) and os.path.getsize(auc_output) > 0:
            # 验证文件内容
            try:
                import loompy as lp
                with lp.connect(auc_output) as ds:
                    if 'RegulonsAUC' in ds.ca and len(ds.ca.RegulonsAUC) > 0:
                        print(f"发现有效的AUCell结果文件: {auc_output}，跳过aucell步骤")
                    else:
                        print(f"AUCell结果文件内容不完整，将重新生成")
                        self._run_aucell_step(motif_output, auc_output)
            except Exception as e:
                print(f"读取AUCell结果文件时出错: {str(e)}，将重新生成")
                self._run_aucell_step(motif_output, auc_output)
        else:
            self._run_aucell_step(motif_output, auc_output)
        
        print("SCENIC分析完成！可以使用SCope进行可视化分析: http://scope.aertslab.org")
        
        # 加载结果到内存
        self._load_results()
        
        return {
            "adjacencies": os.path.join(self.output_dir, "adjacencies.csv"),
            "motifs": os.path.join(self.output_dir, "motifs.csv"),
            "auc_result": os.path.join(self.output_dir, "scenic_result.loom")
        }
    
    def _run_grn_step(self, output_file):
        """运行GRNBoost2步骤"""
        # 运行pyscenic grn命令
        grn_cmd = [
            'pyscenic', "grn",
            self.processed_file,
            self.tf_list_file,
            "-o", output_file,
            "--num_workers", str(self.num_workers),
            "--method", "grnboost2"
        ]
        
        print("运行命令:", " ".join(grn_cmd))
        subprocess.run(grn_cmd, check=True)
        print(f"共表达网络结果已保存到: {output_file}")
        
    def _run_ctx_step(self, grn_output, output_file):
        """运行ctx步骤"""
        # 根据物种选择数据库文件
        db_files = []
        motif_annotations_fname = None
        
        # 确定物种类型
        is_mouse = self.species.startswith("mm")
        is_human = self.species.startswith("hg") or self.species.startswith("GRCh")
        
        print(f"当前分析物种: {self.species} ({'小鼠' if is_mouse else '人类' if is_human else '未知'})")
        
        # 检查db_folder中的文件
        for f in os.listdir(self.db_folder):
            if f.endswith('.feather'):
                # 根据指定的物种选择匹配的数据库文件
                if is_human and any(name in f.lower() for name in ['human', 'hg19', 'hg38', 'grch']):
                    print(f"检测到人类数据库文件: {f}")
                    db_files.append(os.path.join(self.db_folder, f))
                elif is_mouse and any(name in f.lower() for name in ['mouse', 'mm9', 'mm10']):
                    print(f"检测到小鼠数据库文件: {f}")
                    db_files.append(os.path.join(self.db_folder, f))
        
        # 如果没有找到数据库文件，抛出错误
        if not db_files:
            raise ValueError(f"在{self.db_folder}中未找到匹配{self.species}的数据库文件")
        
        # 检查是否有匹配的motif注释文件
        for f in os.listdir(self.db_folder):
            if f.endswith('.tbl'):
                if is_human and 'hgnc' in f.lower():  # 人类
                    # 优先使用v10版本
                    if 'v10' in f:
                        print(f"使用人类motif注释文件(v10): {f}")
                        motif_annotations_fname = os.path.join(self.db_folder, f)
                        break
                    # 记住v9版本，但继续搜索v10
                    elif motif_annotations_fname is None:
                        print(f"找到人类motif注释文件(v9): {f}")
                        motif_annotations_fname = os.path.join(self.db_folder, f)
                elif is_mouse and 'mgi' in f.lower():  # 小鼠
                    # 优先使用v10版本
                    if 'v10' in f:
                        print(f"使用小鼠motif注释文件(v10): {f}")
                        motif_annotations_fname = os.path.join(self.db_folder, f)
                        break
                    # 记住v9版本，但继续搜索v10
                    elif motif_annotations_fname is None:
                        print(f"找到小鼠motif注释文件(v9): {f}")
                        motif_annotations_fname = os.path.join(self.db_folder, f)
        
        # 如果没有找到匹配的注释文件，尝试下载
        if not motif_annotations_fname:
            print(f"未找到匹配的motif注释文件，尝试下载...")
            motif_annotations_fname = download_motif_annotations(self.db_folder, self.species)
            if not motif_annotations_fname:
                raise ValueError(f"无法下载motif注释文件，请手动下载并放置在{self.db_folder}目录中")
        
        # 运行pyscenic ctx命令
        ctx_cmd = [
            'pyscenic', "ctx",
            grn_output,
            *db_files,
            "--annotations_fname", motif_annotations_fname,
            "--expression_mtx_fname", self.processed_file,
            "--output", output_file,
            "--num_workers", str(self.num_workers),
            "--nes_threshold", str(self.nes_threshold),
            "--all_modules",
            "--mode", "custom_multiprocessing"
        ]
        
        print("运行命令:", " ".join(ctx_cmd))
        subprocess.run(ctx_cmd, check=True)
        print(f"顺式调控元件富集结果已保存到: {output_file}")
        
    def _run_aucell_step(self, motif_output, output_file):
        """运行aucell步骤"""
        # 运行pyscenic aucell命令
        aucell_cmd = [
            'pyscenic', "aucell",
            self.processed_file,
            motif_output,
            "-o", output_file,
            "--num_workers", str(self.num_workers),
            "--auc_threshold", str(self.auc_threshold),
            "--cell_id_attribute", "CellID",
            "--gene_attribute", "Gene"
        ]
        
        print("运行命令:", " ".join(aucell_cmd))
        subprocess.run(aucell_cmd, check=True)
        print(f"AUCell结果已保存到: {output_file}")
    
    def run_rss(self):
        """运行RSS分析，计算调控子特异性得分"""
        try:
            print("开始计算调控子特异性得分(RSS)...")
            
            if not hasattr(self, 'auc_mtx') or self.auc_mtx.empty:
                print("错误：AUC矩阵不可用，无法计算RSS。请先运行SCENIC分析。")
                return
            
            # 检查auc_mtx是否有数据
            if self.auc_mtx.shape[0] == 0 or self.auc_mtx.shape[1] == 0:
                print("错误：AUC矩阵为空，无法计算RSS。")
                return
            
            # 检查是否已存在rss_results
            if hasattr(self, 'rss_results') and self.rss_results is not None and not self.rss_results.empty:
                print("RSS结果已存在，跳过计算。")
                return
            
            # 加载分组信息
            cell_type = None
            if self.anndata is not None:
                if 'cell_type' in self.anndata.obs:
                    cell_type = self.anndata.obs['cell_type']
                    print("使用'cell_type'列进行RSS分析")
                elif 'louvain' in self.anndata.obs:
                    cell_type = self.anndata.obs['louvain']
                    print("使用'louvain'列进行RSS分析")
                elif 'leiden' in self.anndata.obs:
                    cell_type = self.anndata.obs['leiden']
                    print("使用'leiden'列进行RSS分析")
                else:
                    print("未找到细胞类型注释，尝试使用leiden算法进行聚类...")
                    # 使用leiden算法进行聚类
                    try:
                        import scanpy as sc
                        
                        # 确保anndata已经经过预处理
                        if 'X_pca' not in self.anndata.obsm:
                            print("执行PCA降维...")
                            sc.pp.pca(self.anndata)
                        
                        if 'neighbors' not in self.anndata.uns:
                            print("计算邻居图...")
                            sc.pp.neighbors(self.anndata)
                        
                        print("执行leiden聚类...")
                        sc.tl.leiden(self.anndata)
                        cell_type = self.anndata.obs['leiden']
                        print(f"leiden聚类完成，识别出{len(cell_type.unique())}个细胞群体")
                    except Exception as e:
                        print(f"leiden聚类失败: {str(e)}")
                        print("错误：未找到细胞类型注释且聚类失败。请在AnnData对象中添加'cell_type'、'louvain'或'leiden'列。")
                        return
            else:
                print("错误：AnnData对象不可用，无法获取细胞分组信息。")
                return
            
            # 确保cell_type索引与auc_mtx索引相匹配
            try:
                cell_type = cell_type.loc[self.auc_mtx.index]
            except:
                print("警告：细胞ID不匹配，尝试使用位置索引...")
                if len(cell_type) == len(self.auc_mtx):
                    cell_type.index = self.auc_mtx.index
                else:
                    print(f"错误：细胞数量不匹配 (AUC矩阵: {len(self.auc_mtx)}，细胞类型注释: {len(cell_type)})。")
                    return
            
            # 计算每个细胞类型的平均AUC值
            cell_types = cell_type.unique()
            regulons = self.auc_mtx.columns
            
            # 创建一个DataFrame来存储RSS结果
            self.rss_results = pd.DataFrame(index=regulons, columns=cell_types)
            
            # 计算每个细胞类型的平均AUC值
            for ct in cell_types:
                cell_indices = cell_type == ct
                if sum(cell_indices) > 0:  # 确保该类型有细胞
                    self.rss_results[ct] = self.auc_mtx.loc[cell_indices].mean()
            
            # 计算RSS (Regulon Specificity Score)
            # 1. z-score标准化每个调控子在不同细胞类型中的表达
            self.rss_results = self.rss_results.subtract(self.rss_results.mean(axis=1), axis=0).divide(self.rss_results.std(axis=1), axis=0)
            
            # 2. 将负值设为0（只保留高于平均值的部分）
            self.rss_results[self.rss_results < 0] = 0
            
            # 3. 归一化使每行的和为1
            row_sums = self.rss_results.sum(axis=1)
            self.rss_results = self.rss_results.divide(row_sums, axis=0)
            
            # 4. 替换NaN为0
            self.rss_results.fillna(0, inplace=True)
            
            # 保存RSS结果到文件
            rss_file = os.path.join(self.output_dir, "rss_scores.csv")
            self.rss_results.to_csv(rss_file)
            print(f"RSS计算完成。结果已保存至: {rss_file}")
            
        except Exception as e:
            import traceback
            print(f"计算RSS时出错: {str(e)}")
            print(traceback.format_exc())
            self.rss_results = None
    
    def _load_results(self):
        """加载计算结果到内存"""
        if hasattr(self, 'results_loaded') and self.results_loaded:
            return

        print("加载分析结果以供可视化...")
        try:
            # 加载adjacencies.csv
            adj_file = os.path.join(self.output_dir, "adjacencies.csv")
            if os.path.exists(adj_file):
                self.adjacencies = pd.read_csv(adj_file)
                print(f"已加载{len(self.adjacencies)}个TF-目标基因关系")
            else:
                print(f"警告：找不到adjacencies文件: {adj_file}")
                self.adjacencies = None

            # 加载motifs.csv
            motif_file = os.path.join(self.output_dir, "motifs.csv")
            if os.path.exists(motif_file):
                self.motifs = pd.read_csv(motif_file)
                print(f"已加载{len(self.motifs)}个潜在调控子")
            else:
                print(f"警告：找不到motifs文件: {motif_file}")
                self.motifs = None

            # 加载scenic_result.loom
            loom_file = os.path.join(self.output_dir, "scenic_result.loom")
            if os.path.exists(loom_file):
                try:
                    # 使用h5py而不是loompy读取loom文件
                    import h5py
                    with h5py.File(loom_file, 'r') as f:
                        # 检查是否有AUC矩阵
                        if 'col_attrs' in f and 'RegulonsAUC' in f['col_attrs']:
                            # 提取AUC矩阵
                            auc_values = f['col_attrs']['RegulonsAUC'][()]
                            cell_ids = f['col_attrs']['CellID'][()] if 'CellID' in f['col_attrs'] else np.arange(len(auc_values))
                            
                            # 转换为pandas DataFrame
                            # 从auc_values中提取列名
                            auc_dtypes = auc_values.dtype
                            regulon_names = auc_dtypes.names if hasattr(auc_dtypes, 'names') else [f"Regulon_{i}" for i in range(auc_values.shape[1])]
                            
                            # 将structured array转换为普通数组
                            if hasattr(auc_dtypes, 'names'):
                                # 对于structured array
                                auc_data = np.column_stack([auc_values[name] for name in regulon_names])
                            else:
                                # 对于普通数组
                                auc_data = auc_values
                            
                            # 创建DataFrame
                            self.auc_mtx = pd.DataFrame(auc_data, index=cell_ids, columns=regulon_names)
                            print(f"已加载AUCell矩阵，包含{self.auc_mtx.shape[0]}个细胞和{self.auc_mtx.shape[1]}个调控子")
                            
                            # 加载RSS分数，如果存在
                            if 'row_attrs' in f and 'RegulonsSpecificityScores' in f['row_attrs']:
                                rss_values = f['row_attrs']['RegulonsSpecificityScores'][()]
                                if hasattr(rss_values.dtype, 'names'):
                                    rss_names = rss_values.dtype.names
                                    rss_data = {name: rss_values[name] for name in rss_names}
                                    self.rss_results = pd.DataFrame(rss_data)
                                    print(f"已加载调控子特异性得分")
                                else:
                                    print("警告：RSS数据格式不支持，需要重新计算")
                        else:
                            print("警告：loom文件中未找到RegulonsAUC数据")
                            self.auc_mtx = None
                except Exception as e:
                    print(f"读取loom文件时出错: {str(e)}")
                    import traceback
                    print(traceback.format_exc())
                    self.auc_mtx = None
            else:
                print(f"警告：找不到loom文件: {loom_file}")
                self.auc_mtx = None

            # 如果需要，加载RSS结果
            rss_file = os.path.join(self.output_dir, "rss_scores.csv")
            if os.path.exists(rss_file) and (not hasattr(self, 'rss_results') or self.rss_results is None):
                try:
                    self.rss_results = pd.read_csv(rss_file, index_col=0)
                    print(f"已加载RSS结果文件: {rss_file}")
                except Exception as e:
                    print(f"读取RSS文件时出错: {str(e)}")
                    self.rss_results = None
            
            self.results_loaded = True
        except Exception as e:
            print(f"加载结果时出错: {str(e)}")
            import traceback
            print(traceback.format_exc())
            self.results_loaded = False
    
    def plot_rss_heatmap(self, n_regulons: int = 20):
        """绘制调控子特异性得分热图"""
        try:
            if not hasattr(self, 'rss_results') or self.rss_results is None:
                print("错误：RSS结果不可用，无法绘制热图。这可能是因为pySCENIC未生成RegulonsSpecificityScores数据。")
                print("请检查SCENIC分析是否完整，或使用其他工具（如SCope）进行可视化。")
                return
            
            if self.rss_results.empty:
                print("错误：RSS结果数据为空，无法绘制热图")
                return
        
            # 获取前N个最具特异性的调控子
            top_regulons = self.rss_results.loc[self.rss_results.sum(axis=1).sort_values(ascending=False).index[:min(n_regulons, len(self.rss_results))]]
        
            # 绘制热图
            plt.figure(figsize=(10, n_regulons//2))
            sns.heatmap(top_regulons, cmap="viridis", linewidths=.5)
            plt.title(f"Top {n_regulons} Regulon Specificity Scores")
            plt.tight_layout()
        
            # 保存图像
            rss_plot_fname = os.path.join(self.output_dir, "rss_heatmap.png")
            plt.savefig(rss_plot_fname, dpi=300)
            plt.close()
            print(f"RSS热图已保存至: {rss_plot_fname}")
        
        except Exception as e:
            print(f"绘制RSS热图时出错: {str(e)}")
    
    def plot_network(self, top_n: int = 20, tf: str = None):
        """绘制基因调控网络可视化图"""
        try:
            if not hasattr(self, 'adjacencies') or self.adjacencies.empty:
                print("警告：调控关系数据不可用，无法绘制网络图")
                return
            
            try:
                import networkx as nx
            except ImportError:
                print("正在安装networkx包...")
                subprocess.check_call([sys.executable, "-m", "pip", "install", "networkx"])
                import networkx as nx
        
            # 根据权重过滤边
            if tf:
                # 仅显示指定TF的网络
                if tf not in set(self.adjacencies['TF']):
                    print(f"警告: 转录因子 {tf} 不在网络中")
                    return
                    
                edges_df = self.adjacencies[self.adjacencies['TF'] == tf].sort_values('importance', ascending=False).head(top_n)
            else:
                # 显示整体网络中权重最高的边
                edges_df = self.adjacencies.sort_values('importance', ascending=False).head(top_n)
                
            if edges_df.empty:
                print("警告：没有找到满足条件的边，无法绘制网络图")
                return
        
            # 创建网络
            G = nx.DiGraph()

            # 添加边
            for _, row in edges_df.iterrows():
                G.add_edge(row['TF'], row['target'], weight=row['importance'])

            # 计算节点大小 (转录因子节点更大)
            tfs = set(edges_df['TF'])
            node_sizes = [1000 if node in tfs else 300 for node in G.nodes()]

            # 绘制网络
            plt.figure(figsize=(14, 10))

            # 改进的布局选择策略
            node_count = len(G.nodes())

            tf_nodes = [node for node in G.nodes() if node in tfs]
            target_nodes = [node for node in G.nodes() if node not in tfs]
            pos = nx.bipartite_layout(G, tf_nodes, align='vertical') 

            # 绘制节点 - 改进颜色和透明度
            nx.draw_networkx_nodes(G, pos, 
                                node_size=node_sizes, 
                                node_color=['#FF6B6B' if node in tfs else '#4ECDC4' for node in G.nodes()],
                                alpha=0.9,
                                edgecolors='black',
                                linewidths=1)

            # 绘制边 - 改进样式
            edge_weights = [G[u][v]['weight'] * 0.02 for u, v in G.edges()]  # 调整权重系数
            nx.draw_networkx_edges(G, pos, 
                                width=edge_weights, 
                                alpha=0.7, 
                                arrowsize=10, 
                                arrowstyle='->',
                                connectionstyle='arc3,rad=0.1',
                                edge_color='#555555')

            # 绘制标签 - 改进样式
            nx.draw_networkx_labels(G, pos, 
                                font_size=9, 
                                font_family='sans-serif',
                                font_weight='bold')

            # 添加图例（可选）
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='#FF6B6B', alpha=0.9, label='Transcription Factors'),
                Patch(facecolor='#4ECDC4', alpha=0.9, label='Targets')
            ]
            plt.legend(handles=legend_elements, loc='upper right')
        
            # 设置标题
            if tf:
                plt.title(f"Gene Regulatory Network for {tf} (top {top_n} targets)",fontsize=14, fontweight='bold', pad=20)
            else:
                plt.title(f"Gene Regulatory Network (top {top_n} interactions)",fontsize=14, fontweight='bold', pad=20)
        
            plt.axis('off')
            plt.tight_layout()
        
            # 保存图像
            if tf:
                network_plot_fname = os.path.join(self.output_dir, f"grn_network_{tf}.png")
            else:
                network_plot_fname = os.path.join(self.output_dir, "grn_network_top.png")
                
            plt.savefig(network_plot_fname, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"网络图已保存至: {network_plot_fname}")
        except Exception as e:
            print(f"绘制网络图时出错: {str(e)}")

def main():
    """主函数，解析命令行参数并执行分析"""
    parser = argparse.ArgumentParser(description="GRN - 基于pySCENIC的基因调控网络分析工具")
    
    # 基本参数
    parser.add_argument("-i", "--input", required=True, help="输入文件路径，支持h5ad/csv/tsv格式")
    parser.add_argument("-o", "--output", default="./scenic_results", help="输出目录")
    parser.add_argument("--tf-list", default=None, help="转录因子列表文件")
    parser.add_argument("--db-folder", required=True, help="RcisTarget数据库文件夹")
    parser.add_argument("--species", default="hg19", help="分析物种，如hg19(人类)或mm10(小鼠)")
    # SCENIC参数
    parser.add_argument("--num-workers", type=int, default=None, help="并行计算的工作进程数")
    parser.add_argument("--auc-threshold", type=float, default=0.05, help="AUCell阈值")
    parser.add_argument("--nes-threshold", type=float, default=3.0, help="富集得分阈值")
    
    # 预处理参数
    parser.add_argument("--normalize", type=str,default='False', help="是否进行总量归一化")
    parser.add_argument("--log", type=str,default='False', help="是否进行log1p转换")
    parser.add_argument("--scale", type=str,default='False', help="是否进行缩放")
    
    # 可视化参数
    parser.add_argument("--top_regulons", type=int, default=20, help="显示前N个调控子")
    parser.add_argument("--top_targets", type=int, default=20, help="每个转录因子显示前N个靶基因")
    parser.add_argument("--plot_tf", default=None, help="为特定转录因子绘制网络图")
    
    args = parser.parse_args()
    from distutils.util import strtobool
    args.log = bool(strtobool(args.log))
    args.scale = bool(strtobool(args.scale))
    args.normalize = bool(strtobool(args.normalize))
    if args.plot_tf == 'None':
        args.plot_tf = None
    # 创建SCENIC分析器
    analyzer = SCENICAnalyzer(
        input_file=args.input,
        output_dir=args.output,
        tf_list_file=args.tf_list,
        db_folder=args.db_folder,
        num_workers=args.num_workers,
        auc_threshold=args.auc_threshold,
        nes_threshold=args.nes_threshold,
        species=args.species,
        normalize = bool(args.normalize),
        log = bool(args.log),
        scale = bool(args.scale),
        config = os.environ.get('CDesk_config')
    )
    
    # 加载数据
    analyzer.load_data()
    
    # 预处理
    analyzer.preprocess()
    
    # 运行SCENIC分析
    results = analyzer.run_scenic()
    
    # 对pyscenic的结果进行rss分析
    analyzer.run_rss()
    
    # 可视化
    try:
        # 先检查RSS分析是否成功
        has_rss = hasattr(analyzer, 'rss_results') and analyzer.rss_results is not None and not analyzer.rss_results.empty
        
        # 绘制RSS热图
        try:
            analyzer.plot_rss_heatmap(n_regulons=args.top_regulons)
        except Exception as e:
            print(f"RSS热图绘制失败: {str(e)}")
            print("由于pySCENIC未能生成调控子特异性得分，您可以尝试使用SCope进行可视化分析")
        
        # 绘制网络图
        try:
            analyzer.plot_network(top_n=args.top_targets, tf=args.plot_tf)
        except Exception as e:
            print(f"网络可视化失败: {str(e)}")
    except Exception as e:
        print(f"可视化过程中出现错误: {str(e)}")
    
    # 总结分析结果
    print("\nSCENIC分析完成！")
    print(f"- 共表达网络已保存至: {os.path.join(args.output, 'adjacencies.csv')}")
    print(f"- 富集分析结果已保存至: {os.path.join(args.output, 'motifs.csv')}")
    print(f"- AUCell结果已保存至: {os.path.join(args.output, 'scenic_result.loom')}")
    
    if has_rss:
        print(f"- RSS热图已保存至: {os.path.join(args.output, 'rss_heatmap.png')}")
    else:
        print("- RSS分析未成功，未生成热图")
    
    print("\n要进行进一步的交互式可视化，请考虑使用SCope工具: http://scope.aertslab.org")

if __name__ == "__main__":
    main() 

