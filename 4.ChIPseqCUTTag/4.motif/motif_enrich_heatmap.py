# %%
import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=UserWarning)
from typing import List
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mpt
import matplotlib.cm as mpm
from matplotlib.patches import PathPatch
from matplotlib.font_manager import fontManager
from itertools import product, combinations, chain
from tqdm import tqdm
from math import ceil
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess
import threading
import traceback
import time
import random
import sys
import re
import os
import concurrent.futures as fu
import argparse
import logging
import json
mpl.rc('pdf', fonttype=42)

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

# Execute the command
def exec_cmd(cmd):
    global logger
    if isinstance(cmd, (list, tuple)):
        cmd = ' '.join(cmd)
    #logger.info('exec cmd: {}'.format(cmd))
    pipe = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        universal_newlines=True
    )
    output_lines = []
    while True:
        buff = pipe.stdout.readline()
        if (buff == '') and (pipe.poll() != None): # No output, and poll != 0 (0: successful, 1: error)
            break
        buff = buff.strip()
        if buff == '':
            continue
        line = buff.rstrip()
        output_lines.append(line)
        #logger.info(buff)
    retcode = pipe.poll()
    if retcode != 0:
        # Execute the command wrongly, print the output
        err_text = "\n".join(output_lines)
        logger.error(f"Command failed (exit code {retcode}): {cmd}")
        logger.error("=== Begin command output ===")
        logger.error(err_text)
        logger.error("===  End command output  ===")
        os._exit(retcode)

# Record the bed file row number
def region_num(file_path: str):
    num = 0
    with open(file_path) as f:
        for line in f:
            if line.strip() == '':
                continue
            num += 1
    return num

# Multi tasks
def multi_task(func_list: list, arg_list: List[dict], task_num: int = None):
    task_all = list()
    task_data = tuple(zip(func_list, arg_list))
    bar = tqdm(total=len(task_data))
    if task_num is None:
        task_num = len(func_list)
    
    def task_done(f: fu.Future):
        nonlocal bar
    
        bar.update(1)
        # Error thread
        e = f.exception()
        if e:
            logger.error(traceback.format_exc())
    
    with fu.ThreadPoolExecutor(max_workers=task_num) as executor:
        for f, args in task_data:
            # Submit task
            task_temp = executor.submit(f, **args)
            # Callback task
            task_temp.add_done_callback(task_done)
            task_all.append(task_temp)
        fu.wait(task_all)
    bar.close()

# XX vs XX.txt. Conduct a Permutation Test for region overlap.
# Perform statistical analysis on the overlap between two genomic region sets (A and B).
def enrich(file_bed_a: str, file_bed_b: str, dir_result: str, title: str, species: str):
    if os.path.exists(os.path.join(dir_result, "{}.txt".format(title))):
        return
    r_str = """
suppressMessages(library(regioneR))

setwd("{4}")

file_bed_a <- "{0}"
file_bed_b <- "{1}"

A <- toGRanges(
    file_bed_a,
    genome = "{2}"
)
B <- toGRanges(
    file_bed_b,
    genome = "{2}"
)
pt <- overlapPermTest(
    A = A,
    B = B,
    ntimes = 500,
    genome = "{2}",
    mc.cores = 1
)

res <- c(
    filea=file_bed_a,
    fileb=file_bed_b,
    pval=pt$numOverlaps$pval,
    zscore=pt$numOverlaps$zscore
)
write.table(
    res,
    "{3}.txt",
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
)
""".format(file_bed_a, file_bed_b, species, title, dir_result)
    file_r = os.path.join(dir_result, "{}.r".format(title))
    with open(file_r, "w") as f:
        f.write(r_str)
    exec_cmd("Rscript {}".format(file_r))
    exec_cmd("rm -rf {}".format(file_r))

def enrich_main(file_bed_dict: dict, dict_type_name: str,species:str,task_num: int = 12):
    global palette, dir_base, dir_cooc
    # now（MotifAll folder）and base folder
    beds = pd.read_csv(dir_cooc)
    if 'tag' not in beds.columns or 'bed' not in beds.columns or 'group' not in beds.columns:
        print('Need columns: bed,tag,group')
        sys.exit(1)
    if not beds['tag'].is_unique:
        print('Duplicates in tag column')
        sys.exit(1)
    dir_enrich_now = os.path.join(dir_base, dict_type_name)
    os.makedirs(dir_enrich_now, exist_ok=True)
    dir_enrich_base = os.path.join(dir_base, "base")
    os.makedirs(dir_enrich_base, exist_ok=True)
    # Merge each bed again and put in now folder, maybe can skip or merge at the step before. The bed file is large.
    for bed_name in file_bed_dict:
        file_bed_temp = os.path.join(dir_enrich_now, "{}.bed".format(bed_name))
        if not os.path.exists(file_bed_temp):
            exec_cmd("cat {} | cut -f 1-3 | sort -k1,1 -k2,2n | uniq | {} merge -i - > {}".format(
                file_bed_dict[bed_name],
                'bedtools',
                file_bed_temp
            ))
        file_bed_dict[bed_name] = file_bed_temp
    # Merge the CO and OC again and put in base folder (Maybe not necessary)
    bed_dict = dict(zip(beds['tag'], beds['bed']))
    for bed_name in bed_dict:
        file_bed_temp = os.path.join(dir_enrich_base, "{}.bed".format(bed_name))
        if not os.path.exists(file_bed_temp):
            exec_cmd("cat {} | cut -f 1-3 | sort -k1,1 -k2,2n | uniq | {} merge -i - > {}".format(
                bed_dict[bed_name],
                'bedtools',
                file_bed_temp
            ))
        bed_dict[bed_name] = file_bed_temp
    # species，threads #####################
    func_list = list()
    args_list = list()
    file_list = list()
    # Run enrich() in multiple tasks
    for bed_name_a in bed_dict:
        for bed_name_b in file_bed_dict:
            func_list.append(enrich)
            args_list.append(dict(
                file_bed_a=bed_dict[bed_name_a],
                file_bed_b=file_bed_dict[bed_name_b],
                dir_result=dir_enrich_now,
                title="{}_vs_{}".format(bed_name_a, bed_name_b),
                species=species
            ))
            file_list.append("{}_vs_{}.txt".format(bed_name_a, bed_name_b))
    multi_task(func_list, args_list, task_num)

    if not os.path.exists(os.path.join(dir_base, "enrich.{}.txt".format(dict_type_name))):
        # Record the result to a dataframe
        data_res = list()
        for file_single in file_list:
            if not os.path.exists(os.path.join(dir_enrich_now, file_single)):
                continue
            type_name_a, type_name_b = file_single[:-4].split("_vs_", 1)
            with open(os.path.join(dir_enrich_now, file_single)) as f:
                data_temp = {
                    "a": type_name_a,
                    "b": type_name_b
                }
                for line in f:
                    line = line.strip()
                    if line == "":
                        continue
                    line = line.split("\t", 1)
                    data_temp[line[0]] = line[1]
                data_res.append(data_temp)
        df = pd.DataFrame(data_res)
        df["pval"] = df["pval"].astype(float)
        df["-log10(pval)"] = -np.log10(df["pval"])
        df["zscore"] = df["zscore"].astype(float)
        # 
        def get_ratio(se: pd.Series):
            file_temp = "{}.{}.tmp.{}.bed".format(
                os.path.basename(se["filea"]),
                os.path.basename(se["fileb"]),
                random.randint(100000, 999999)
            )
            exec_cmd("{} -wa -a {} -b {} | sort | uniq > {}".format(
                'intersectBed',
                se["filea"],
                se['fileb'],
                file_temp
            ))
            se["ratio"] = round(region_num(file_temp) / region_num(se["filea"]) * 100, 1)
            exec_cmd("rm -rf {}".format(file_temp))
            return se

        df = df.apply(get_ratio, axis=1)
        df.to_csv(
            os.path.join(dir_base, "enrich.{}.txt".format(dict_type_name)),
            sep="\t",
            index=False
        )
        pass
    else:
        df = pd.read_table(os.path.join(dir_base, "enrich.{}.txt".format(dict_type_name)))
        df["zscore"] = df["zscore"].astype(float)

    # According to TF order
    df = df[df["b"].isin(tuple(file_bed_dict.keys()))]
    df["b_index"] = df["b"].map(lambda d: tuple(file_bed_dict.keys()).index(d))
    df.sort_values(["a", "b_index"], inplace=True)
    mapping = beds.set_index("tag")["group"].to_dict()
    df["tag"] = df["a"].map(mapping)

    df_all = df.copy()

    def get_size(d: float):
        if d < 0:
            return 10
        elif d < 5:
            return 100
        elif d < 10:
            return 200
        elif d < 20:
            return 350
        else:
            return 500
    # zscore size
    df_all["zscore_size"] = df_all["zscore"].map(get_size)
    # Plot
    for cooc_type, df in df_all.groupby("tag"):
        fig = plt.figure(figsize=(len(df["a"].unique()) * 0.5, len(df["b"].unique()) * 0.7))
        grid = plt.GridSpec(20, 7, hspace=0.1, wspace=0.1)
        vmin = 0
        vmax = 80
        cmap = mpl.colors.LinearSegmentedColormap.from_list("www", ("lightblue", "red"))
        norm = mpl.colors.Normalize(
            vmin=vmin,
            vmax=vmax
        )

        ax_color = fig.add_subplot(
            grid[-3:-2, :3], 
            xticklabels=[],
            yticklabels=[]
        )
        plt.colorbar(
            mpm.ScalarMappable(
                norm=norm,
                cmap=cmap
            ),
            cax=ax_color,
            orientation='horizontal',
            label="% of targets",
        )
        ax_color.set_xticks(
            [vmin, vmax],
            ["{:.1f}".format(vmin), "{:.1f}".format(vmax)]
        )
        ax_color.spines[:].set_visible(True)

        ax_size = fig.add_subplot(
            grid[-4:, 4:], 
            xticklabels=[],
            yticklabels=[]
        )
        
        ss = [
            (10, "0"),
            (100, ""),
            (200, ""),
            (350, ""),
            (500, "20"),
        ]
        ax_size.scatter(
            [(i/len(ss))**1.4 for i in range(len(ss))],
            [0 for _ in range(len(ss))],
            s=list(map(lambda d: d[0], ss)),
        )
        for i, si in enumerate(ss):
            ax_size.text(
                x=(i/len(ss))**1.4,
                y=-1,
                s=si[1],
                ha="center",
                va="top"
            )
        ax_size.set_xlim(-0.1, 0.9)
        ax_size.set_ylim(-2, 1.5)
        ax_size.set_xlabel("EnrichmentZscore")
        ax_size.set_xticks([])
        ax_size.set_yticks([])
        ax_size.spines[:].set_visible(False)

        ax = fig.add_subplot(
            grid[:-6, :],
        )
        ax.scatter(
            df["a"].to_numpy(),
            df["b"].to_numpy(),
            s=df["zscore_size"].to_numpy(),
            c=df["ratio"].to_numpy(),
            cmap=cmap,
            norm=norm
        )
        ax.set_title("Sensitive Enrich")
        ax.set_ylabel("")
        ax.set_xlabel("")
        
        xlim = ax.get_xlim()
        ax.set_xlim((xlim[0] - 0.5, xlim[1] + 0.5))
        ylim = ax.get_ylim()
        ax.set_ylim((ylim[0] - 0.5, ylim[1] + 0.5))

        plt.savefig(
            os.path.join(dir_base, "enrich.zscore.{}.{}.pdf".format(dict_type_name, cooc_type)),
            bbox_inches="tight",
            dpi=200
        )
        plt.close()

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ATAC sample replicate correlation")
    parser.add_argument("--motif_bed",required=True, type=str, help="The motif bed files directory")
    parser.add_argument("--bed",required=True, type=str, help="The input bed files directory")
    parser.add_argument("-o",required=True,type=str, help="The output directory")
    parser.add_argument("-t",default=12,type=int, help="The number of threads")
    parser.add_argument("--species",default='mm10',type=str, help="Species")

    return parser.parse_args()

args = parse_arguments()

dir_base = os.path.abspath(args.o)
dir_cooc = args.bed
dir_motif = args.motif_bed
threads = args.t
species = args.species

os.makedirs(dir_base, exist_ok=True)
os.chdir(dir_base)

# Color palette
global palette
palette = "#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75\n#1E7CAF\n#46489A\n#0F231F\n#1187CD\n#b7704f\n#596032\n#2570a1\n#281f1d\n#de773f\n#525f42\n#2585a6\n#2f271d\n#c99979\n#5f5d46\n#1b315e\n#1d1626\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75".split("\n")

df = pd.read_csv(dir_motif, sep=",")
if 'motif' not in df.columns or 'bed' not in df.columns:
    print('Need columns: motif,bed')
    sys.exit(1)

filtered_bed_dict = dict(zip(df['motif'], df['bed']))

if not filtered_bed_dict:
    print("Error: Can not find any specified motif corresponding bed file")
    sys.exit(1)

# Run
enrich_main(
    file_bed_dict=filtered_bed_dict,
    dict_type_name="MotifAll",
    task_num=threads,
    species = species
)

print('Done, you can see the results now.')
