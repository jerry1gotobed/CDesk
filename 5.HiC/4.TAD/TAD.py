import fanc
import subprocess
import os
from datetime import datetime
import sys
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import pybedtools
from matplotlib_venn import venn2
import argparse
import shutil

def is_tool_available(name):
    """检查命令是否存在于 PATH 中"""
    return shutil.which(name) is not None

tools = ['fanc','bedtools']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

current_dir = os.path.dirname(os.path.abspath(__file__))
# Parameter
parser = argparse.ArgumentParser(description="HiC compartment analysis")
# Required parameter
parser.add_argument('-i', '--input', type=str, required=True, help="The information file")
parser.add_argument('-o', '--output', type=str, required=True, help="The output directory")
# Optional parameter
#parser.add_argument('-t', '--thread', type=int, default=20, help="The number of threads, default: 20")
parser.add_argument('--insulation_window', type=int, default=1000000, help="Window sizes to calculate insulation score in base pairs, default: 1000000")
parser.add_argument('--directionality_window', type=int, default=1000000, help="Window sizes to calculate directionality index in base pairs, default: 1000000")
parser.add_argument('--expand', type=int, default=1000000, help="Window sizes of boundary expansion in base pairs to calulate insulation score and DI near/on boundary, default: 1000000")

args = parser.parse_args()
input_file = args.input # "/mnt/linzejie/CDesk_test/data/5.HiC/4.TAD/test.csv"
output_dir = args.output # '/mnt/linzejie/tmp/hic/tad'
insulation_window = args.insulation_window # 1000000
directionality_window = args.directionality_window # 1000000
expand = args.expand # 1000000

input_file = pd.read_csv(input_file)
required_columns = ['file', 'tag']
assert all(col in input_file.columns for col in required_columns), f"Require column {required_columns} for the input dataframe"
tags = list(input_file['tag'].fillna('').astype(str).str.strip())
hic_samples = list(input_file['file'].fillna('').str.strip())

def human_format(num, precision=2, lowercase=False):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0

    num = round(num, precision)

    if num - int(num) == 0:
        num = int(num)
    num_str = '{}{}'.format(num, ['', 'k', 'M', 'G', 'T', 'P'][magnitude])
    return num_str if not lowercase else num_str.lower()

# Check hic samples resolutions to be the same
resolution_check = []
for hic_sample in hic_samples:
    temp = fanc.load(hic_sample)
    resolution_check.append(temp.bin_size)
if len(set(resolution_check)) > 1:
    print('The hic samples resolutions are not the same')
    sys.exit(1)

resolution = resolution_check[0]

tmp_dir = os.path.join(output_dir,'insulation_boundary_di')
contact_dir = os.path.join(output_dir,'contact')
img_dir = os.path.join(output_dir,'img')
tad_dir = os.path.join(output_dir,'tad')
os.makedirs(output_dir,exist_ok=True)
os.makedirs(tmp_dir,exist_ok=True)
os.makedirs(contact_dir,exist_ok=True)
os.makedirs(img_dir,exist_ok=True)
os.makedirs(tad_dir,exist_ok=True)

for hic_sample, tag in zip(hic_samples, tags):
    # fanc insulation score
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Calculate {tag} insulation score")
    command = f"fanc insulation {hic_sample} {os.path.join(tmp_dir,tag+'.insulation')} -i -g -w {insulation_window} -o bed"
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {tag} Insulation score calculation done")

    # fanc boundary, ignore the error message
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Analyze {tag} boundaries")
    command = f'''
                fanc boundaries {os.path.join(tmp_dir,tag+'.insulation_'+human_format(insulation_window).lower() + 'b'+'.bed')} {os.path.join(tmp_dir,tag+'.insulation_boundaries')} 2>&1 | \
                awk '
                /Traceback/ {{in_traceback=1}}
                in_traceback && /AttributeError/ {{in_traceback=0; next}}
                !in_traceback
                '
                '''
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {tag} Boundaries analysis done")

    # fanc directionality index
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Calculate {tag} directionality index")
    command = f"fanc directionality {hic_sample} {os.path.join(tmp_dir,tag+'.directionality')}  -w {directionality_window} -o bed "
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {tag} Directionality index calculation done")

    # Get TAD between boundary
    command = f'''
        sed -e "1d" {os.path.join(tmp_dir, tag + ".insulation_boundaries")} \
        | paste {os.path.join(tmp_dir, tag + ".insulation_boundaries")} - \
        | awk '$1 == $7 {{print $1 "\\t" $2 "\\t" $8}}' \
        > {os.path.join(tad_dir, tag + ".tad.bed")}
    '''
    subprocess.run(command,check=True,shell=True)

    # Expand boundary
    command = f'''
        awk 'BEGIN{{OFS="\\t"}} {{
            $2 = ($2 - {expand} < 0) ? 0 : $2 - {expand};
            $3 = $3 + {expand};
            print $0, ($2+$3)/2
        }}' {os.path.join(tmp_dir, tag + ".insulation_boundaries")} > {os.path.join(tmp_dir, tag + ".insulation_boundaries_expand")}
    '''
    subprocess.run(command, check=True, shell=True)

    # Get insulation score for each boundary (on or near)
    command = f'''
        bedtools window -w 0 -a {os.path.join(tmp_dir, tag + ".insulation_boundaries_expand")} -b {os.path.join(tmp_dir,tag+'.insulation_'+human_format(insulation_window).lower() + 'b'+'.bed')}|awk '{{bdr=($2+$3)/2;pos=($9+$10)/2;print (pos-bdr)"\t"$12}}' > {os.path.join(tmp_dir,tag+'.insulation_distance2boundary.txt')}
    '''
    subprocess.run(command, check=True, shell=True)

    # Get directionality index for each boundary (on or near)
    command = f'''
        bedtools window -w 0 -a {os.path.join(tmp_dir, tag + ".insulation_boundaries_expand")} -b {os.path.join(tmp_dir,tag+'.directionality_'+human_format(insulation_window).lower() + 'b'+'.bed')}|awk '{{bdr=($2+$3)/2;pos=($9+$10)/2;print (pos-bdr)"\t"$12}}' > {os.path.join(tmp_dir,tag+'.direction_distance2boundary.txt')}
    '''
    subprocess.run(command, check=True, shell=True)

    # Get contact count matrix
    # fanc dump Hi-C_early_2-cell.100000N_fanc.hic --only-intra -u > 1.txt
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Get {tag} contact matrix of each chromosom")
    command = f"fanc dump {hic_sample} --only-intra -u > {os.path.join(contact_dir,tag+'_contact.txt')}"
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")

    command = f"""awk '
        $1 == $4 {{
            chr = $1
            filename = "{contact_dir}/{tag}_" chr "_contact.txt"
            print $3 "\t" $6 "\t" $7 >> filename
            close(filename)
        }}' "{os.path.join(contact_dir, tag+'_contact.txt')}"
    """
    subprocess.run(command, check=True, shell=True)
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {tag} contact matrix extraction done")

    # Aggregate plot
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Aggregate tad for {tag} ")
    command = f"""
        fanc aggregate {hic_sample} {os.path.join(tad_dir, tag + ".tad.bed")} \
        -p {os.path.join(img_dir,tag+'_tad_aggregate.pdf')} --tads --tad-strength {os.path.join(tad_dir, tag + "_tad_strength.bed")}
    """
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {tag} Tad aggregate finished")

    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {tag} Sample finished")

# Compare TAD
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Tad compare")
command = ['Rscript',os.path.join(current_dir,'tad.R'),output_dir,','.join(tags),str(resolution)]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")

# Plot
for i in os.listdir(tad_dir):
    if i.endswith('.tadcompare.csv'):
        # 创建新图形
        plt.figure(figsize=(10, 8))
        prefix = i.replace('.tadcompare.csv','')
        output_pdf = os.path.join(img_dir,prefix+'.tadcompare.pdf')
        temp = pd.read_csv(os.path.join(tad_dir,i))
        sizes = list(temp['Type'].value_counts(sort=False).sort_index())
        labels = list(temp['Type'].value_counts(sort=False).sort_index().index)
        explode = [0] * len(labels)
        if 'strength-down' in labels:
            explode[labels.index('strength-down')] = 0.1
        if 'strength-up' in labels:
            explode[labels.index('strength-up')] = 0.15
        with PdfPages(output_pdf) as pdf:
            # Pie chart
            wedges, texts, autotexts = plt.pie(
                sizes,
                explode=explode,
                shadow=True,
                startangle=90,
                counterclock=False,
                wedgeprops={'linewidth': 0.7, 'edgecolor': 'black'},
                pctdistance=0.8,
                autopct='%1.1f%%'
            )
            
            plt.title(prefix.replace('_vs_',' vs '), fontsize=20, color="black")
            plt.legend(
                wedges,
                labels,
                loc="center left",
                bbox_to_anchor=(1, 0.5),
                fontsize=12
            )
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
        
            # ------------------------- Venn -------------------------
            plt.figure(figsize=(10,8))
            prefix.split('_vs_')
            # 加载BED文件
            f1 = pybedtools.BedTool(os.path.join(tad_dir,prefix.split('_vs_')[0]+".tad.bed"))
            f2 = pybedtools.BedTool(os.path.join(tad_dir,prefix.split('_vs_')[1]+".tad.bed"))
            
            # 计算重叠区域
            n_f1 = len(f1)
            n_f2 = len(f2)
            n_common = len(f1.intersect(f2, u=True))
            n_common_1 = len(f2.intersect(f1, u=True))
            n_common = min(n_common,n_common_1)

            # 绘制韦恩图
            v = venn2(subsets=(n_f1-n_common, n_f2-n_common, n_common),
                    set_labels=(None, None),
                    set_colors=("#DE582B", "#1868B2"),
                    alpha=0.7)
            
            # 设置标题
            plt.title("Overlap of TADs", fontsize=16)
            plt.text(-0.7, 0.3, prefix.split('_vs_')[0], fontsize=14, ha='center', va='center')
            plt.text(0.7, -0.3, prefix.split('_vs_')[1], fontsize=14, ha='center', va='center')
            
            # 保存韦恩图到PDF
            pdf.savefig(bbox_inches='tight')
            plt.close()

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Tad compare finished") 
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Finished, you can check the results now")

