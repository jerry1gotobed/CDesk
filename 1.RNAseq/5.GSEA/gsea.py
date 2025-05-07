import argparse
import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import pickle

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run GSEA analysis")
    subparser = parser.add_subparsers(dest='mode', required=True, help="Choose to analyze or plot")
    # analysis
    gsea_parser = subparser.add_parser('gsea', help='Perform the gsea analysis')
    gsea_parser.add_argument("--rnk",required=True,type=str, help="The ranked gene list file (e.g. gene_name(upper_case) log2FC(descending)) or the CDesk deg analysis result csv file")
    gsea_parser.add_argument("-o",required=True,  type=str, help="Directory to save the results")
    gsea_parser.add_argument("-s","--species",default="human",type=str, help="Specify the species",choices=['human','mouse'])
    gsea_parser.add_argument("--preprocess", action="store_true", help="Preprocess from CDesk DEG analysis result")   
    gsea_parser.add_argument("--permutation_num",default=1000 , type=int, help="Number of permutation, default: 1000")
    gsea_parser.add_argument("-t","--thread",default=4 , type=int, help="Number of threads, default: 4")
    # plot
    plot_parser = subparser.add_parser('plot', help='Specify the GSEA term plot')
    plot_parser.add_argument("-i","--input",required=True, help="The GSEA analysis output directory")
    plot_parser.add_argument("--term", help="Specify the terms to plot, separate by ,")
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()
    
    if args.mode == 'gsea':
        msig = gp.Msigdb()
        if args.species=='human':
            gmt = msig.get_gmt(category='h.all', dbver="2024.1.Hs")
        elif args.species=='mouse':
            gmt = msig.get_gmt(category='h.all', dbver="2024.1.Mm") # Mm

        if not os.path.exists(args.o):
            os.makedirs(args.o)
        rank = pd.read_csv(args.rnk, sep=None, engine='python',)
        if args.preprocess:
            rank = rank.loc[:, ['gene_name','log2FC']]
            rank = rank.sort_values(by='log2FC',ascending=False)
            rank['gene_name'] = rank['gene_name'].str.upper()
        gs_res = gp.prerank(rnk = rank, gene_sets = gmt, seed = 6, permutation_num = args.permutation_num,threads=args.thread)
        gs_res.res2d = gs_res.res2d.sort_values(by='FDR q-val',ascending=True)
        df = gs_res.res2d
        df.to_csv(f'{args.o}/gsea_result_sort.csv')
        terms = gs_res.res2d.Term
        axs = gs_res.plot(terms[:5], show_ranking=False, legend_kws={'loc': (1.05, 0)},)
        fig = axs.get_figure()
        fig.savefig(f"{args.o}/gsea_plot_top5_significant_path.pdf", format="pdf",bbox_inches="tight")
        with open(f"{args.o}/gs_res.pkl", "wb") as file:
            pickle.dump(gs_res, file)
        print(f"运行完成,对象已保存到{args.o}/gs_res.pkl")

    if args.mode == 'plot':
        with open(f"{args.input}/gs_res.pkl", "rb") as file:
            gs_res = pickle.load(file)
        terms = [x.strip() for x in args.term.split(",")]
        axs = gs_res.plot(terms, show_ranking=False, legend_kws={'loc': (1.05, 0)}, )
        fig = axs.get_figure()
        fig.savefig(f"{args.input}/gsea_plot_specified_terms.pdf", format="pdf", bbox_inches="tight")
        print(f"运行完成,图片已保存到{args.input}/gsea_plot_specified_terms.pdf")

if __name__ == "__main__":
    main()
