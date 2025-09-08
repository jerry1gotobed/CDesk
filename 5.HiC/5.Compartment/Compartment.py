import cooler
import argparse
import pandas as pd
import os
import sys
import subprocess
import shutil
import json
import numpy as np
import cooltools
import bioframe
import warnings
from cytoolz import merge
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
from datetime import datetime
import warnings

def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['cooler','cooltools']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

current_dir = os.path.dirname(os.path.abspath(__file__))
# Parameter
parser = argparse.ArgumentParser(description="HiC compartment analysis")
# Required parameter
parser.add_argument('-i', '--input', type=str, required=True, help="The information file")
parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution to specify")
parser.add_argument('-o', '--output', type=str, required=True, help="The output directory")
# Optional parameter
parser.add_argument('-t', '--thread', type=int, default=20, help="The number of threads, default: 20")
parser.add_argument('--GC', type=str, default='no', help="The GC file to provide, default: no")
parser.add_argument('-s', '--species', type=str, default='no', help="The species, need to specify if no GC file provided, default: no")
parser.add_argument('--width', type=float, default=10, help="Default: 10")
parser.add_argument('--height', type=float, default=8, help="Default: 8")

# 解析参数
args = parser.parse_args()

meta = args.input
meta = pd.read_csv(meta)
if ('file' not in meta.columns) or ('tag' not in meta.columns):
    print('The meta file has to have columns: file tag')
    sys.exit(1)
input_files = list(meta['file'])
resolution = args.resolution
species = args.species
thread = args.thread
ATAC = args.GC
if ATAC != 'no':
    if (not os.path.exists(ATAC)):
        print(f'{ATAC} does not exit!')
        sys.exit(1)
    gc_df = pd.read_csv(ATAC,sep='\t',header=None)
    if gc_df.shape[1] != 4:
        print("Error: The file does not have exactly 4 columns: chrom,start,end,GC")
    gc_df.columns = ['chrom', 'start', 'end', 'GC']

config = os.getenv("CDesk_config")
with open(config, "r") as f:
    config = json.load(f)

if species != 'no':
    fasta = config['data'][species]['fasta']
    if (not os.path.exists(fasta)):
        print(f'{fasta} does not exit!')
        sys.exit(1)

if (ATAC == 'no') and (species == 'no'):
    print('Must provide a parameter: track/species')
    sys.exit(1)

output_dir = args.output
cool_dir = os.path.join(output_dir,'cool')
expected_cis_dir = os.path.join(output_dir,'expected_cis')
tmp_dir = os.path.join(output_dir,'tmp')
interaction_dir =  os.path.join(output_dir,'interaction')
img_dir = os.path.join(output_dir,'img')
AB_dir = os.path.join(output_dir,'AB')
os.makedirs(output_dir, exist_ok=True)
os.makedirs(cool_dir, exist_ok=True)
os.makedirs(expected_cis_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)
os.makedirs(interaction_dir, exist_ok=True)
os.makedirs(img_dir, exist_ok=True)
os.makedirs(AB_dir, exist_ok=True)

width = args.width
height = args.height

def saddleplot(
    track,
    saddledata,
    n_bins,
    width,
    height,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Generate a saddle plot.
    Parameters
    ----------
    track : pd.DataFrame
        See cooltools.digitize() for details.
    saddledata : 2D array-like
        Saddle matrix produced by `make_saddle`. It will include 2 flanking
        rows/columns for outlier signal values, thus the shape should be
        `(n+2, n+2)`.
    cmap : str or matplotlib colormap
        Colormap to use for plotting the saddle heatmap
    scale : str
        Color scaling to use for plotting the saddle heatmap: log or linear
    vmin, vmax : float
        Value limits for coloring the saddle heatmap
    color : matplotlib color value
        Face color for margin bar plots
    fig : matplotlib Figure, optional
        Specified figure to plot on. A new figure is created if none is
        provided.
    fig_kws : dict, optional
        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional
        Passed on to `ax.imshow()`
    margin_kws : dict, optional
        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional
        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object
        Specify a subregion of a figure to using a GridSpec.
    Returns
    -------
    Dictionary of axes objects.
    """

#     warnings.warn(
#         "Generating a saddleplot will be deprecated in future versions, "
#         + "please see https://github.com/open2c_examples for examples on how to plot saddles.",
#         DeprecationWarning,
#     )

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    # Old version
    # hist = np.bincount(x, minlength=len(binedges) + 1)

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Layout
    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    # Figure
    if fig is None:
        fig_kws_default = dict(figsize=(width, height))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    # Heatmap
    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    # Margins
    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
    # left margin hist
    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])

    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
    )

    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    # top margin hist
    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])

    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
    )

    plt.xlim(lo, hi)
    # plt.ylim(plt.ylim())  # correct
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

#     # Colorbar
    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        # cb.set_ticks(np.arange(vmin, vmax + 0.001, 0.5))
        # # do linspace between vmin and vmax of 5 segments and trunc to 1 decimal:
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        print('cbar')

        cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    # extra settings
    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)
    return grid

# Check the input files and whether the resolution is OK
for file in input_files:
    if (not os.path.exists(file)):
        print(f'{file} does not exit!')
        sys.exit(1)
    # if the resolution os OK
    cooler_paths = cooler.fileops.list_coolers(file)
    resolutions = [int(path.split('/')[-1]) for path in cooler_paths] if len(cooler_paths) > 1 else [cooler.Cooler(file).binsize]
    hic_resolution = min([int(path.split('/')[-1]) for path in cooler_paths]) if len(cooler_paths) > 1 else cooler.Cooler(file).binsize
    if (resolution % hic_resolution) != 0:
        print(f'The {resolution} cannot be evenly divided by the minimal resolution {hic_resolution}, please change the resolution')
        sys.exit(1)
    if (len(resolutions)>1) and (resolution not in resolutions):
        print(f'The {resolution} not in mcool resolutions {resolutions}, please change the resolution')
        sys.exit(1)

samples = []
for file in input_files:
    cooler_paths = cooler.fileops.list_coolers(file)
    resolutions = [int(path.split('/')[-1]) for path in cooler_paths] if len(cooler_paths) > 1 else [cooler.Cooler(file).binsize]
    hic_resolution = min([int(path.split('/')[-1]) for path in cooler_paths]) if len(cooler_paths) > 1 else cooler.Cooler(file).binsize
    # Balance if needed
    cool = file if len(resolutions) == 1 else f"{file}::resolutions/{resolution}"
    cool_temp = cooler.Cooler(cool)
    is_balanced = "weight" in cool_temp.bins().columns
    # Transfer to the cool with specified resolution
    # cool -> norm_cool ; cool -> mcool(norm)
    sample_prefix = os.path.splitext(os.path.basename(file))[0]
    print('________________________________________________________________________________')
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Analyze {file}")
    samples.append(sample_prefix)
    if resolution in resolutions:
        if not is_balanced:
            print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Balance the matrix")
            tmp_cool = os.path.join(cool_dir,os.path.basename(file))
            shutil.copy(file, tmp_cool)
            cool = tmp_cool if len(resolutions) == 1 else f"{tmp_cool}::resolutions/{resolution}"
            command = ["cooler", "balance", cool]
            try:
                subprocess.run(command, check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                print(f"Error：{e.stderr.decode()}")
            print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Matrix balance done")
    else: # cool - > mcool
        #command = f'cooler zoomify -n {thread} --balance -r {resolution}N -o {cool_dir}/{sample_prefix}_{resolution}N.mcool {file}'
        print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Modify the resolution")
        command = ['cooler', 'zoomify','-n', str(thread),'--balance','-r',str(resolution),'-o',os.path.join(cool_dir, f'{sample_prefix}_{resolution}N.mcool'),file]
        try:
            subprocess.run(command, check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
        cool = os.path.join(cool_dir, f'{sample_prefix}_{resolution}N.mcool::resolutions/{resolution}')
        print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Resolution modification done")
    
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Compartment analysis")
    if ATAC == 'no':
        command = f"cooler dump --header -t bins {cool} | cut -f1-3 > {os.path.join(tmp_dir,sample_prefix+'.bins.'+str(resolution)+'.tsv')}"
        try:
            subprocess.run(command, check=True,shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
        command = f"cooltools genome gc {os.path.join(tmp_dir,sample_prefix+'.bins.'+str(resolution)+'.tsv')} {fasta} > {os.path.join(tmp_dir,sample_prefix+'.gc.'+str(resolution)+'.tsv')}"
        try:
            subprocess.run(command, check=True,shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
        command = f"cooltools eigs-cis -o {os.path.join(expected_cis_dir,sample_prefix+'.eigs.'+str(resolution))} --phasing-track {os.path.join(tmp_dir,sample_prefix+'.gc.'+str(resolution)+'.tsv')} --n-eigs 1 {cool}"
        try:
            subprocess.run(command, check=True,shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
        gc_df = pd.read_csv(os.path.join(tmp_dir,sample_prefix+'.gc.'+str(resolution)+'.tsv'),sep='\t')
    else:
        command = f"cooltools eigs-cis -o {os.path.join(expected_cis_dir,sample_prefix+'.eigs.'+str(resolution))} --phasing-track {ATAC} --n-eigs 1 {cool}"
        try:
            subprocess.run(command, check=True,shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Compartment analysis done")
    
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Saddle plot")
    clr = cooler.Cooler(cool)

    # PCA eigenvector
    view_df = pd.DataFrame({'chrom': clr.chromnames, 'start': 0, 'end': clr.chromsizes.values, 'name': clr.chromnames})
    cis_eigs = cooltools.eigs_cis(clr, gc_df, view_df=view_df, n_eigs=3)
    eigenvector_track = cis_eigs[1][['chrom','start','end','E1']].dropna()
    valid_chroms = eigenvector_track['chrom'].unique()
    view_df_filtered = view_df[view_df['chrom'].isin(valid_chroms)].copy()
    # expect-cis
    cvd = cooltools.expected_cis(clr=clr, view_df=view_df_filtered)

    # Data for lot
    Q_LO = 0
    Q_HI = 1
    N_GROUPS = 50
    interaction_sum, interaction_count = cooltools.saddle(clr, cvd, eigenvector_track, 'cis', n_bins=N_GROUPS, qrange=(Q_LO,Q_HI), view_df=view_df_filtered)
    summatrix = pd.DataFrame(interaction_sum)
    countmatrix = pd.DataFrame(interaction_count)

    summatrix.to_csv(os.path.join(interaction_dir, sample_prefix) + '.s_matrix.csv',index=False,sep="\t")
    countmatrix.to_csv(os.path.join(interaction_dir, sample_prefix) + '.c_matrix.csv',index=False,sep="\t")
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        saddleplot(eigenvector_track, interaction_sum/interaction_count, N_GROUPS, width=width,height=height,qrange=(Q_LO,Q_HI), cbar_kws={'label':'average observed/expected contact frequency'})
    plt.savefig(os.path.join(img_dir, sample_prefix) + '.saddle_plot.pdf')

tag = list(meta['tag'].astype(str))
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Compute the compartment strength and analyze A/B transit for all samples")
command = ['Rscript',os.path.join(current_dir,'TransferAB.R'),output_dir,str(resolution),','.join(samples),','.join(tag),str(width),str(height)]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")

command = ['Rscript',os.path.join(current_dir,'Strength.R'),output_dir,','.join(samples),','.join(tag),str(width),str(height)]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done, you can check the results now")
