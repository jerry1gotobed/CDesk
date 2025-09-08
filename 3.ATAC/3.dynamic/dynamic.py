import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=UserWarning)
# %%
from sklearn import preprocessing
import concurrent.futures as fu
#from pyg2plot import Plot, JS
from matplotlib.font_manager import fontManager
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mpt
import matplotlib.cm as mpm
from itertools import product
from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess
import traceback
import argparse
import sys
import os
from typing import List
import chardet
import json

mpl.use("Agg")
mpl.rc('pdf', fonttype=42)
#mpl.rcParams['font.sans-serif'] = ["Arial"]

# 得到样本信息，生成需要的系列信息
def get_samples(file_meta: str,important_index: int):
    """ 获取样本的信息
    """
    # 检查是否存在，格式，重复，两个分组
    head_tail_names = set(("head", "tail"))
    # 检查是否存在，格式，重复，两个分组
    assert os.path.exists(file_meta), "meta not exists: {}".format(file_meta)
    with open(file_meta, 'rb') as f:
        result = chardet.detect(f.read(1000))  # 读取前 1000 字节检测
    encoding = result['encoding']
    # 读取文件
    df = pd.read_csv(file_meta, sep=None, engine='python', encoding=encoding)
    df.fillna("", inplace=True)
    assert tuple(df.columns) == ("sample", "tag_time", "tag_group"), "meta file columns error!"
    assert len(df) == len(df["sample"].unique()), "sample dumplates !"
    assert len(set(df["tag_group"].unique()) - head_tail_names) <= 2, "Only less than 2 tag groups accepted"

    # 两个组的time一致，记录每个组的sample是什么
    samples_group = dict()
    samples_head_tail = dict()
    samples_time = list(filter(
        lambda d: d.strip() != "",
        df[~df["tag_group"].isin(head_tail_names)]["tag_time"].unique()
    ))
    for group_name, df_temp in df.groupby("tag_group"):
        if group_name not in head_tail_names:
            assert tuple(df_temp["tag_time"]) == tuple(samples_time), "tag time order not same !"
            samples_group[group_name] = df_temp["sample"].to_list()
            continue
        samples_head_tail[group_name] = df_temp["sample"].to_list()

    # '0111111111'生成这样的
    def padding_false(add_head: bool, len_samples: int):
        list_now = [[ 1 for _ in range(i)] for i in range(1, len_samples+1)]
        len_final = max(map(
            lambda d: len(d),
            list_now
        ))
        for index in range(len(list_now)):
            gap = [0 for _ in range(len_final - len(list_now[index]))]
            if add_head:
                list_now[index] = tuple(gap + list_now[index])
            else:
                list_now[index] = tuple(list_now[index] + gap)
            list_now[index] = "".join(map(str, list_now[index]))
        list_now = list_now[:-1]
        if add_head:
            list_now.reverse()
        return list_now

    # need
    samples_need = df["sample"].to_list()

    # pat_samples
    samples_pat_samples = {
        k: samples_head_tail.get("head", list()) + samples_group[k] + samples_head_tail.get("tail", list()) for k in samples_group
    }

    # pat_last
    samples_pat_last = {
        k: samples_group[k][-1] for k in samples_group
    }

    # pat
    ks = list(samples_group.keys())
    if len(ks) > 1:
        assert len(samples_pat_samples[ks[0]]) == len(samples_pat_samples[ks[1]]), "sample has not same length: {}, {}".format(*samples_pat_samples)
    samples_len = len(samples_pat_samples[ks[0]])
    samples_pat = [
        padding_false(True, samples_len),
        padding_false(False, samples_len)
    ]
    samples_pat_type = ["CO", "OC"]

    return samples_pat_type, samples_pat, samples_pat_samples, samples_need, samples_pat_last, samples_time, samples_head_tail, samples_group, list(samples_group.keys())[important_index]

# 查找对应样本的结尾文件，缺失或重复则报错
def get_ext_files(dir_path: str, ext: str,samples_need):
    #nonlocal samples_need
    res = dict()
    end_str = "{}".format(ext)
    for root, _, files in os.walk(dir_path):
        for file_single in files:
            if not file_single.endswith(end_str):
                continue
            file_name = file_single[:-len(end_str)]
            if file_name not in samples_need:
                continue
            assert file_name not in res, "{} {} duplicate: {}, {}".format(file_name, end_str, res[file_name], os.path.join(dir_path, root, file_single))
            res[file_name] = os.path.join(root, file_single)

    sample_gap = set(samples_need) - set(res.keys())
    assert len(sample_gap) == 0, "sample missing: {}".format(sample_gap)
    return list(res.values())

# 得到一个merge.ov.bed，记录每个样本文件在合并bed的交集次数
def bed_combine(bed_file: List[str], ext_len: int, result_dir: str, tmp_dir: str, config:str,need_head: bool = True, file_name_tag: str = "", force: bool = True):
    """ 将bed合并
    """
    global samples_head_tail, samples_group, important_group
    res_txt_file = os.path.join(result_dir, "{}merge.ov.bed".format(file_name_tag))
    tmp_cat_file = os.path.join(tmp_dir, "cat.bed")
    tmp_merge_file = os.path.join(tmp_dir, "merge.bed")

    # 已运行过就不再运行
    if os.path.exists(res_txt_file) and (not force):
        return res_txt_file

    # merge
    os.system("cat {} > {}".format(
        " ".join(bed_file),
        tmp_cat_file
    ))
    os.system("cut -f 1-3 {} | sort -k1,1 -k2,2n | {} merge -i - > {}".format(
        tmp_cat_file,
        config['software']['bedtools'],
        tmp_merge_file
    ))

    #得到每个文件与merge的交集次数
    cmd = ""
    header = list()
    for bed_file_temp in bed_file:
        cmd = "{}| {} -c -a - -b {}".format(
            cmd,
            config['software']['intersectBed'],
            bed_file_temp
        )
        if need_head:
            header.append(os.path.basename(bed_file_temp)[:-ext_len].split(".",1)[-1])
    cmd = "cat {} {} > {}".format(
        tmp_merge_file,
        cmd,
        res_txt_file
    )
    os.system(" ".join(cmd) if isinstance(cmd, list) else cmd)

    #
    with open(tmp_cat_file, "w") as f:
        if need_head:
            f.write("\t".join(["chr", "start", "end"] + header))
        f.write("\n")

    os.system("cat {} {} > {}".format(
        tmp_cat_file,
        res_txt_file,
        tmp_merge_file
    ))
    os.system("cat {} > {}".format(
        tmp_merge_file,
        res_txt_file
    ))

    os.system(" ".join(["rm", "-rf", tmp_cat_file]))
    os.system(" ".join(["rm", "-rf", tmp_merge_file]))
    if len(list(samples_group.keys())) > 1:
        less_important = list(set(samples_group.keys()) - {important_group})[0]
        column_order = ["chr", "start", "end"] + samples_head_tail['head'] +  samples_group[important_group] + samples_group[less_important] + samples_head_tail['tail']
    else:
        column_order = ["chr", "start", "end"] + samples_head_tail['head'] +  samples_group[important_group] + samples_head_tail['tail']
    # 读取 res_txt_file
    df = pd.read_csv(res_txt_file, sep="\t")
    # 按指定顺序重新排列列
    df = df[column_order]
    # 保存为新的 res_txt_file
    df.to_csv(res_txt_file, sep="\t", index=False)

    return res_txt_file

# 计算bw文件bed区域的平均信号值
def calc_bed_signal(bed_file: str, bw_file: str, signal_dir: str, tmp_dir: str,config: str, force: bool = True):
    """ 计算bw文件bed区域的平均信号值
    """
    res_signal_file = os.path.join(signal_dir, os.path.basename(bw_file).replace(".bw", ".txt"))
    tmp_signal_file = os.path.join(tmp_dir, os.path.basename(bw_file).replace(".bw", ".signal.txt"))
    tmp_bed_file = os.path.join(tmp_dir, os.path.basename(bw_file).replace(".bw", ".calc.bed"))

    # 已运行过就不再运行
    if os.path.exists(res_signal_file) and (not force):
        return res_signal_file

    os.system("cat -n " + bed_file + """ | awk -F "\t" '{if(NR>1){print $2"\t"$3"\t"$4"\t"$1}}' > """ + tmp_bed_file)

    os.system(" ".join([config['software']['bigWigAverageOverBed'],bw_file, tmp_bed_file, tmp_signal_file, "> /dev/null 2>&1"]))
    os.system("""cat """+tmp_signal_file+""" | sort -k1,1n | awk -F "\t" '{print $5}' > """ + res_signal_file)

    with open(res_signal_file) as f:
        res_data = np.asarray(list(map(
            lambda dd: float(dd.strip()),
            f
        )))

    os.system(" ".join(["rm", "-rf", tmp_signal_file]))
    os.system(" ".join(["rm", "-rf", tmp_bed_file]))
    return res_data

# 根据是否和000011,1110000一致来判断所属CO和OC
def get_cluster(se: pd.Series):
    global samples_pat_type, samples_pat, samples_pat_samples, samples_pat_last

    for ii in samples_pat_samples:
        pat_now = "".join(map(
            lambda d: "0" if d == 0 else "1",
            se[samples_pat_samples[ii]]
        ))
        se['last_{}'.format(ii)] = "0" if se[samples_pat_last[ii]] == 0 else "1"
        se['cluster_{}'.format(ii)] = "Other"

        for i in range(len(samples_pat)):
            for index, pat_temp in enumerate(samples_pat[i], start=1):
                if pat_temp != pat_now:
                    continue
                se['cluster_{}'.format(ii)] = "{}{}".format(samples_pat_type[i], index)
                break
            if se['cluster_{}'.format(ii)] != "Other":
                break
    return se

# 分类
def cluster_samples(signal_file: str, samples_pat_type: list, samples_pat: list, samples_pat_samples: dict, samples_need: list, samples_pat_last: list, bw_dir: str, result_dir: str, signal_dir: str, tmp_dir: str,config:str, force: bool = True):
    """ 将样本分类
    """
    res_num_file = os.path.join(result_dir, "cluster.num.{}".format(os.path.basename(signal_file)))
    res_signal_file = os.path.join(result_dir, "cluster.signal.{}".format(os.path.basename(signal_file)))
    res_signal_file_ori = os.path.join(result_dir, "cluster.signal_ori.{}".format(os.path.basename(signal_file))) # 没有标准化

    # 已运行过就不再运行
    if os.path.exists(res_signal_file) and (not force):
        return res_signal_file, res_signal_file_ori

    # 读取信息,分类
    def get_cluster(se: pd.Series):
        nonlocal samples_pat_type, samples_pat, samples_pat_samples, samples_pat_last
        for ii in samples_pat_samples:
            pat_now = "".join(map(
                lambda d: "0" if d == 0 else "1",
                se[samples_pat_samples[ii]]
            ))
            se['last_{}'.format(ii)] = "0" if se[samples_pat_last[ii]] == 0 else "1"
            se['cluster_{}'.format(ii)] = "Other"
            for i in range(len(samples_pat)):
                for index, pat_temp in enumerate(samples_pat[i], start=1):
                    if pat_temp != pat_now:
                        continue
                    se['cluster_{}'.format(ii)] = "{}{}".format(samples_pat_type[i], index)
                    break
                if se['cluster_{}'.format(ii)] != "Other":
                    break
        return se

    if not os.path.exists(res_num_file) or force:
        df_peak_signal = pd.read_table(signal_file)
        df_peak_signal = df_peak_signal.apply(get_cluster, axis=1)
        df_peak_signal.to_csv(
            res_num_file,
            index=False,
            sep="\t"
        )
    else:
        df_peak_signal = pd.read_table(res_num_file)

    # bw file
    files_bw_dict = dict()
    for i in samples_need:
        file_temp = os.path.join(bw_dir, "{}.bw".format(i))
        if not os.path.exists(file_temp):
            print("No bw file: {}, {}".format(i, bw_dir))
            sys.exit(1)
        files_bw_dict[i] = file_temp


    # 计算每一列的signal
    # original
    for i in samples_need:
        df_peak_signal[i] = calc_bed_signal(
            bed_file=res_num_file,
            bw_file=files_bw_dict[i],
            signal_dir=signal_dir,
            tmp_dir=tmp_dir,
            force=force,
            config = config
        ).tolist()
    df_peak_signal.to_csv(
        res_signal_file_ori,
        index=False,
        sep="\t"
    )
    # normalized
    for i in samples_need:
        df_peak_signal[i] = preprocessing.scale(df_peak_signal[i]).tolist()
        df_peak_signal[i] = df_peak_signal[i] - df_peak_signal[i].min()
    df_peak_signal.to_csv(
        res_signal_file,
        index=False,
        sep="\t"
    )

    return res_signal_file, res_signal_file_ori

# 画图，单个
def draw_image(signal_file: str, sample_group_now: str, min_percent: float, max_percent: float, samples_pat_samples: dict, result_dir: str, remove_cluster: list, sample_time: list, force: bool = True):
    """ 画一张图
    """
    title_name = sample_group_now
    res_file = os.path.join(result_dir, "{}{}.pdf".format(title_name, ".remove_cluster" if len(remove_cluster) > 0 else ""))

    # 已运行过就不再运行
    if os.path.exists(res_file) and (not force):
        return res_file

    df = pd.read_table(signal_file)
    cluster_name = "cluster_{}".format(sample_group_now)
    df = df[df[cluster_name] != "Other"]

    # 去除不要的类
    if len(remove_cluster) > 0:
        remove_cluster = set(remove_cluster)
        df = df[~df[cluster_name].isin(remove_cluster)]

    df.sort_values([cluster_name], inplace=True)
    df.set_index(cluster_name, inplace=True)
    df_temp = df[samples_pat_samples[sample_group_now]]

    df_temp.reset_index(inplace=True)
    df_temp_co = df_temp[df_temp[cluster_name].str.startswith("CO")]
    df_temp_oc = df_temp[df_temp[cluster_name].str.startswith("OC")]
    df_temp = pd.concat((
        df_temp_co.sort_values(cluster_name),
        df_temp_oc.sort_values(cluster_name),
    ))
    # df_temp.drop("score", inplace=True, axis=1)
    df_temp.set_index(cluster_name, inplace=True)
    cutoff_max = -100
    cutoff_min = 100
    for i in samples_pat_samples[sample_group_now]:
        cutoff_max = max(cutoff_max, df_temp[i].quantile(max_percent))
        cutoff_min = min(cutoff_max, df_temp[i].quantile(min_percent))


    x_n = len(samples_pat_samples[sample_group_now])
    fig = plt.figure(figsize=(x_n*0.5, 14))
    grid = plt.GridSpec(16, 24, hspace=3, wspace=0)
    main_ax = fig.add_subplot(grid[:-1, 1:])
    y_ax = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
    x_ax = fig.add_subplot(grid[-1, 5:19], yticklabels=[], sharex=y_ax)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("www", ("#030104", "#d57237", "#ecec9e", "#ffffcc"))
    sns.heatmap(
        df_temp,
        cmap=cmap,
        vmax=cutoff_max,
        vmin=cutoff_min,
        ax=main_ax,
        # cbar=None,
        cbar_ax=x_ax,
        cbar_kws={
            "orientation": "horizontal",
        },
        rasterized=True
    )
    #x_ax.spines[:].set_visible(True)
    for spine in x_ax.spines.values():
        spine.set_visible(True)
    x_ax.set_xticks(
        [cutoff_min, cutoff_max],
        ["{:.1f}".format(cutoff_min), "{:.1f}".format(cutoff_max)]
    )
    x_ax.set_xlabel("Normalized Signal")
    labels = dict(enumerate(df_temp.index))
    main_ax.set_title(title_name)
    main_ax.xaxis.set_ticklabels(ticklabels=sample_time, rotation=0)
    _ = main_ax.set_ylabel("")
    _ = main_ax.yaxis.set_tick_params(left=False)

    label_dict = dict()
    label_max = 0
    for lp in labels:
        lt = labels[lp]
        if lt not in label_dict:
            label_dict[lt] = list()
        label_dict[lt].append(lp)
        label_max = max(label_max, lp)
    label_dict_labels = dict(map(
        lambda d: ((max(label_dict[d]) - min(label_dict[d]))/2 + min(label_dict[d]), d),
        label_dict
    ))
    # g = label_max * 0.005
    g = min(map(lambda d: len(d), label_dict.values())) / 3
    l = g
    label_dict_labels_posi = list()
    label_dict_labels_lab = list()
    y_ax_x_posi = sum(y_ax.get_xlim())/2
    for d in label_dict_labels.keys():
        t = l+(d-l)*2
        if l > t:
            continue
        d_t = [l, t]
        # print(d_t)
        l = t + g*2
        if (d_t[1]-d_t[0]) <= g:
            continue
        label_dict_labels_posi.append(d)
        label_dict_labels_lab.append(label_dict_labels[d])
        y_ax.plot(
            [y_ax_x_posi, y_ax_x_posi],
            d_t,
            c='black'
        )
    y_ax.yaxis.set_ticks([])
    for p, l in zip(label_dict_labels_posi, label_dict_labels_lab):
        y_ax.text(-0.5, p, l, ha="right", va="center")
    _ = y_ax.yaxis.set_tick_params(left=False)
    y_ax.axis("off")
    y_ax.invert_yaxis()
    plt.savefig(res_file, dpi=500, bbox_inches='tight')
    plt.close()

# 画图，组合
def draw_image_by_one(signal_file: str, important_group: str, min_percent: float, max_percent: float, samples_pat_samples: dict, samples_group: dict, result_dir: str, remove_cluster: list,sample_time: list, samples_head_tail: dict, head_tail:dict,force: bool = True):
    """ 画一张图
    """
    title_name = "Align_by_{}".format(important_group)
    res_file = os.path.join(result_dir, "{}{}.pdf".format(title_name, ".remove_cluster" if len(remove_cluster) > 0 else ""))


    # 已运行过就不再运行
    if os.path.exists(res_file) and (not force):
        return res_file

    df = pd.read_table(signal_file)
    cluster_name = "cluster_{}".format(important_group)
    df = df[df[cluster_name] != "Other"]

    # 去除不要的类
    if len(remove_cluster) > 0:
        remove_cluster = set(remove_cluster)
        df = df[~df[cluster_name].isin(remove_cluster)]

    df.sort_values([cluster_name], inplace=True)
    # df.set_index(cluster_name, inplace=True)

    important_less_group = list(samples_pat_samples.keys())
    important_less_group.remove(important_group)
    important_less_group = important_less_group[0]
    samples_now_need = set(samples_pat_samples[important_group]).union(set(samples_pat_samples[important_less_group]))

    df_temp = df
    # 将斜率缩放处理（除以 100，取负值）后返回，作为一个分数
    def get_score1(se: pd.Series):
        y = np.array(se.to_list())
        x = np.array(range(1, y.shape[0]+1))
        theta, _, _, _ = np.linalg.lstsq(np.vstack([x, np.ones(len(x))]).T, y, rcond=None)
        return -round(theta.flatten()[0]/100, 1)
    df_temp['score1'] = df_temp[samples_pat_samples[important_group]].apply(get_score1, axis=1)
    last_key = "last_{}".format(important_less_group)
    df_temp[last_key] = df_temp[last_key].astype("str")
    df_temp['score2'] = (df_temp[cluster_name].map(lambda d: "0" if d.startswith("OC") else "1") != df_temp[last_key]).astype("int")
    # CO是1，less important last 是 0'; OC 是 0， 1；相同的更大

    # df_temp.reset_index(inplace=True)
    df_temp_co = df_temp[df_temp[cluster_name].str.startswith("CO")]
    df_temp_oc = df_temp[df_temp[cluster_name].str.startswith("OC")]
    # df_temp_oc['score'] = -df_temp_oc['score']
    df_temp = pd.concat((
        df_temp_co.sort_values([cluster_name, "score1", "score2"]),
        df_temp_oc.sort_values([cluster_name, "score1", "score2"]),
    ))
    df_temp.set_index(cluster_name, inplace=True)
    df_temp = df_temp[list(samples_now_need)]
    # df_temp.drop(["score1", "score2"], inplace=True, axis=1)
    cutoff_max = -10000
    cutoff_min = 10000
    for i in samples_now_need:
        cutoff_max = max(cutoff_max, df_temp[i].quantile(max_percent))
        cutoff_min = min(cutoff_max, df_temp[i].quantile(min_percent))

    cluster_list = [
        samples_head_tail.get("head", list()),
        samples_group[important_group],
        samples_group[important_less_group],
        samples_head_tail.get("tail", list())
    ]
    cluster_list_name = [
        "head",
        important_group,
        important_less_group,
        "tail"
    ]
    cluster_list_len_sum = sum(map(
        lambda d: len(d),
        cluster_list
    ))
    x_n = cluster_list_len_sum + len(cluster_list)
    grid_min = 2
    fig = plt.figure(figsize=(x_n*0.5, 12))
    grid = plt.GridSpec(16, cluster_list_len_sum*2+grid_min, hspace=1, wspace=0.8)
    main_ax_list = list()
    is_first = True
    for i in range(len(cluster_list)):
        if len(cluster_list[i]) == 0:
            continue
        grid_min_temp = grid_min + len(cluster_list[i]) * 2
        if is_first:
            ax_t = fig.add_subplot(grid[:-2, grid_min:grid_min_temp])
            is_first = False
        else:
            ax_t = fig.add_subplot(grid[:-2, grid_min:grid_min_temp], xticklabels=[], sharey=main_ax_list[0])
        main_ax_list.append(ax_t)
        grid_min = grid_min_temp
    y_ax = fig.add_subplot(grid[:-2, 0:2], xticklabels=[], sharey=main_ax_list[0])
    x_ax = fig.add_subplot(grid[-1, 11:27], yticklabels=[], sharex=y_ax)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("www", ("#030104", "#d57237", "#ecec9e", "#ffffcc"))

    def main_ax_single(i: int, name: str, sample_list_now: list):
        nonlocal df_temp, cmap, cutoff_min, cutoff_max, main_ax_list, x_ax, sample_time, samples_pat_samples

        main_ax = main_ax_list[i]
        df_temp_now = df_temp[sample_list_now]

        if i == 0:
            sns.heatmap(
                df_temp_now,
                cmap=cmap,
                vmax=cutoff_max,
                vmin=cutoff_min,
                ax=main_ax,
                cbar_ax=x_ax,
                cbar_kws={
                    "orientation": "horizontal",
                },
                rasterized=True
            )
        else:
            sns.heatmap(
                df_temp_now,
                cmap=cmap,
                vmax=cutoff_max,
                vmin=cutoff_min,
                ax=main_ax,
                cbar=None,
                rasterized=True
            )

        if name not in ("head", "tail"):
            main_ax.xaxis.set_ticklabels(ticklabels=sample_time, rotation=0)
            _ = main_ax.set_xlabel(name)
        elif name in ("head"):
            main_ax.xaxis.set_ticklabels(ticklabels=head_tail['head'], rotation=0)
            _ = main_ax.set_xlabel("")
        elif name in ("tail"):
            main_ax.xaxis.set_ticklabels(ticklabels=head_tail['tail'], rotation=0)
            _ = main_ax.set_xlabel("")
        else:
            _ = main_ax.set_xlabel("")
            # _ = main_ax.set_xlabel(name)
        _ = main_ax.xaxis.set_tick_params(left=False)
        _ = main_ax.set_ylabel("")
        _ = main_ax.yaxis.set_tick_params(left=False)

    ii = 0
    for i in range(len(cluster_list)):
        if len(cluster_list[i]) == 0:
            continue
        main_ax_single(ii, cluster_list_name[i], cluster_list[i])
        ii += 1

    #x_ax.spines[:].set_visible(True)
    for spine in x_ax.spines.values():
        spine.set_visible(True)
    x_ax.set_xticks(
        [cutoff_min, cutoff_max],
        ["{:.1f}".format(cutoff_min), "{:.1f}".format(cutoff_max)]
    )
    x_ax.set_xlabel("Normalized Signal")
    labels = dict(enumerate(df_temp.index))
    label_dict = dict()
    label_max = 0
    for lp in labels:
        lt = labels[lp]
        if lt not in label_dict:
            label_dict[lt] = list()
        label_dict[lt].append(lp)
        label_max = max(label_max, lp)
    label_dict_labels = dict(map(
        lambda d: ((max(label_dict[d]) - min(label_dict[d]))/2 + min(label_dict[d]), d),
        label_dict
    ))
    # g = label_max * 0.005
    g = min(map(lambda d: len(d), label_dict.values())) / 4
    l = g
    label_dict_labels_posi = list()
    label_dict_labels_lab = list()
    y_ax_x_posi = max(y_ax.get_xlim())
    for d in label_dict_labels.keys():
        t = l+(d-l)*2
        if l > t:
            continue
        d_t = [l, t]
        # print(d_t)
        l = t + g*2
        if (d_t[1]-d_t[0]) <= g:
            continue
        label_dict_labels_posi.append(d)
        label_dict_labels_lab.append(label_dict_labels[d])
        y_ax.plot(
            [y_ax_x_posi, y_ax_x_posi],
            d_t,
            c='black'
        )
    y_ax.yaxis.set_ticks([])
    for p, l in zip(label_dict_labels_posi, label_dict_labels_lab):
        y_ax.text(0, p, l, va="center", ha="right",
                #   fontsize=4
        )
    _ = y_ax.yaxis.set_tick_params(left=False)
    y_ax.axis("off")
    y_ax.invert_yaxis()
    plt.savefig(res_file,dpi=500,bbox_inches='tight')
    # plt.savefig(res_file.replace(".pdf", ".png"), dpi=500, bbox_inches='tight')
    plt.close()

# 生成每个bed的基因和一个总结文件,fail不理解
def get_genes_from_bed(bed_all_file: str, samples_pat_samples: dict, promoter_file: str, tf_file: str,result_dir: str,config:str):

    result_dir_bed = os.path.join(result_dir, "cluster")

    if not os.path.exists(result_dir_bed):
        os.makedirs(result_dir_bed)

    # 获得每个cluster的bed文件
    df = pd.read_table(bed_all_file)
    for i in samples_pat_samples:
        n = i
        i = "cluster_{}".format(i)
        for t, df_temp in df.groupby(i):
            if t == "Other":
                continue
            file_name_temp = os.path.join(result_dir_bed, "{}_{}.bed".format(n, t))
            df_temp[["chr", "start", "end"]].to_csv(
                file_name_temp,
                header=False,
                index=False,
                sep="\t"
            )

    result_dir_gene = bed_to_gene(
        bed_dir=result_dir_bed,
        promoter_file=promoter_file,
        tf_file=tf_file,
        config=config
    )

    return result_dir_bed, result_dir_gene

def bed_to_gene(bed_dir: str, promoter_file: str, tf_file: str,config:str):
    result_dir_gene = "{}_genes".format(bed_dir)
    gene_summary_file = "{}_gene_summary.txt".format(bed_dir)
    if not os.path.exists(result_dir_gene):
        os.makedirs(result_dir_gene)

    # 打开各个bed
    data_summary = list()
    df_tfs = pd.read_table(tf_file)
    tfs = set(df_tfs['Symbol'].to_list())
    del df_tfs
    get_sample_gene = lambda f: set(filter(
        lambda dd: dd != "",
        map(
            lambda d: d.strip(),
            f
        )
    ))

    for file_single in os.listdir(bed_dir):
        file_single = os.path.join(bed_dir, file_single)
        region_num = 0
        with open(file_single) as f:
            for line in f:
                if line.strip() == "":
                    break
                region_num += 1
        tmp = os.path.join(result_dir_gene, "temp.bed")
        type_temp = os.path.basename(file_single)[:-4]
        gene_file_temp = os.path.join(result_dir_gene, "{}.txt".format(type_temp))
        os.system("{} -wa -a {} -b {} | sort | uniq > {}".format(
            config['software']['intersectBed'],
            promoter_file,
            file_single,
            tmp
        ))
        os.system("""cat """ + tmp + """ | awk -F '\t' '{print $5}' | sort | uniq > """ + gene_file_temp)
        os.system("rm -rf {}".format(tmp))
        # tf
        with open(gene_file_temp) as f:
            gene_temp = get_sample_gene(f)
        tf_temp = tfs & gene_temp
        with open(gene_file_temp.replace(".txt", ".tf.txt"), "w") as f:
            f.write("\n".join(tf_temp))
        data_summary.append({
            "type": type_temp,
            "region": region_num,
            "gene": len(gene_temp),
            "tf": len(tf_temp)
        })
    df = pd.DataFrame(data_summary)
    df.sort_values("type", inplace=True)
    df.to_csv(gene_summary_file, index=False, sep="\t")

    return result_dir_gene

warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=UserWarning)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ATAC sample replicate correlation")
    parser.add_argument("--file_meta",required=True, type=str, help="The meta file")
    parser.add_argument("--bed_dir", required=True,type=str, help="Specify the bed file directory")
    parser.add_argument("--bed_summit",default='', type=str, help="Specify the bed file summit")
    parser.add_argument("--bw_dir", required=True,type=str, help="Specify the bw file directory")
    parser.add_argument("--remove_cluster",default='',type=str, help="Input peak file directory")
    parser.add_argument("-o",required=True,type=str, help="The output directory")
    parser.add_argument("--important_index",default=0,type=int, help="Specify the group to align, 1: the latter one in meta file, 0(default): the former one in meta file")
    parser.add_argument("--species", required=True,type=str, help="Species")
    parser.add_argument("--config_path", type=str,help="Path to configuration file")
    return parser.parse_args()

args = parse_arguments()
file_meta=args.file_meta
bed_file_dir=args.bed_dir
bed_file_ext=args.bed_summit
bw_dir=args.bw_dir
remove_cluster=args.remove_cluster
if remove_cluster != '':
    remove_cluster=remove_cluster.split(',')
important_index=args.important_index
species = args.species
config_path = args.config_path
result_dir = args.o

with open(args.config_path, "r") as f:
    config = json.load(f)
tf_file = config['data'][species]['tf_file']
promoter_file = config['data'][species]['promoter_file']

# 控制颜色
min_percent = 0.1
max_percent = 0.7
remove_cluster_null = list()

# 生成结果文件夹
tmp_dir = os.path.join(result_dir, "tmp")
dir_all = [result_dir, tmp_dir]
for i in dir_all:
    if os.path.exists(i):
        continue
    os.makedirs(i)

# 获取样本信息
samples_pat_type, samples_pat, samples_pat_samples, samples_need, samples_pat_last, sample_time, samples_head_tail, samples_group, important_group = get_samples(file_meta, important_index)
sample_time_all = samples_head_tail.get("head", list()) + sample_time + samples_head_tail.get("tail", list())

# 生成文件夹
signal_dir = os.path.join(result_dir, "signal")
if not os.path.exists(signal_dir):
    os.makedirs(signal_dir)
bed_dir_merge = os.path.join(result_dir, "beds_merge")
if not os.path.exists(bed_dir_merge):
    os.makedirs(bed_dir_merge)

# 查找对应样本的结尾文件，缺失或重复则报错
# 获得对应文件路径
print('1. Get sample imformation')
bed_file_ori = get_ext_files(
    dir_path=bed_file_dir,
    ext=bed_file_ext,
    samples_need=samples_need
)

assert len(bed_file_ori) > 0, "No bed files !"
bed_file = bed_file_ori

# 得到一个merge.ov.bed，记录每个样本文件在合并bed的交集次数
print('2. Merge the samples and get the OC/CO state')
signal_num_file_ori = bed_combine(
    bed_file=bed_file,
    ext_len=len(bed_file_ext),
    result_dir=bed_dir_merge,
    tmp_dir=tmp_dir,
    force=True,
    config = config
)

# 将数据分类
print('3. Get the cluster type and calculate the signal strength')
signal_file, signal_file_ori = cluster_samples(
    signal_file=signal_num_file_ori,
    samples_pat_type=samples_pat_type,
    samples_pat=samples_pat,
    samples_pat_samples=samples_pat_samples,
    samples_need=samples_need,
    samples_pat_last=samples_pat_last,
    bw_dir=bw_dir,
    result_dir=bed_dir_merge,
    signal_dir=signal_dir,
    tmp_dir=tmp_dir,
    force=True,
    config = config
)

# 画图(dynamic)
meta = pd.read_csv(file_meta)
head_tag = meta[meta['tag_group'] == 'head']['tag_time'].values[0] if meta[meta['tag_group'] == 'head']['tag_time'].values[0] != '' else meta[meta['tag_group'] == 'head']['sample'].values[0]
tail_tag = meta[meta['tag_group'] == 'tail']['tag_time'].values[0] if meta[meta['tag_group'] == 'head']['tag_time'].values[0] != '' else meta[meta['tag_group'] == 'tail']['sample'].values[0]
sample_time_tags = [head_tag] + sample_time + [tail_tag]
print('4. Plot heatmap')
image_dir = os.path.join(result_dir, "img")
if not os.path.exists(image_dir):
    os.makedirs(image_dir)
for sample_group_now in samples_pat_samples:
    draw_image(
        signal_file=signal_file,
        sample_group_now=sample_group_now,
        min_percent=min_percent,
        max_percent=max_percent,
        samples_pat_samples=samples_pat_samples,
        result_dir=image_dir,
        remove_cluster=remove_cluster_null,
        sample_time=sample_time_tags,
        force=True
    )
    if remove_cluster != '':
        draw_image(
            signal_file=signal_file,
            sample_group_now=sample_group_now,
            min_percent=min_percent,
            max_percent=max_percent,
            samples_pat_samples=samples_pat_samples,
            result_dir=image_dir,
            remove_cluster=remove_cluster,
            sample_time=sample_time_tags,
            force=True
        )

head_tail = meta[meta['tag_group'].isin(['head', 'tail'])].groupby("tag_group")["tag_time"].apply(list).to_dict()
# align one
if len(samples_pat_samples) > 1:
    #print(f"Plot 2 samples, align by {important_group}")
    draw_image_by_one(
        signal_file=signal_file,
        important_group=important_group,
        min_percent=min_percent,
        max_percent=max_percent,
        samples_pat_samples=samples_pat_samples,
        samples_group=samples_group,
        result_dir=image_dir,
        remove_cluster=remove_cluster_null,
        sample_time=sample_time,
        samples_head_tail=samples_head_tail,
        head_tail = head_tail,
        force=True
    )
    if remove_cluster != '':
        #print(f"Plot 2 samples, align by {important_group}, remove unwanted clusters")
        draw_image_by_one(
            signal_file=signal_file,
            important_group=important_group,
            min_percent=min_percent,
            max_percent=max_percent,
            samples_pat_samples=samples_pat_samples,
            samples_group=samples_group,
            result_dir=image_dir,
            remove_cluster=remove_cluster,
            sample_time=sample_time,
            samples_head_tail=samples_head_tail,
            head_tail = head_tail,
            force=True
        )

# 修改输出目录为result_dir参数指定的路径
print('5. Get gene information from bed files')
output_dir = os.path.join(result_dir, 'cluster')
os.makedirs(output_dir, exist_ok=True)
dir_bed, dir_gene = get_genes_from_bed(
    bed_all_file=signal_file,
    samples_pat_samples=samples_pat_samples,
    promoter_file=promoter_file,
    tf_file=tf_file,
    result_dir=result_dir,
    config=config
)

print("运行完成！您可以查看结果了！")
