import argparse
import os
import subprocess
import fanc
import cooler
import sys
import json
import shutil
from datetime import datetime

config_path = os.environ.get('CDesk_config')

def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['fanc','hicConvertFormat','cooler','java','awk','sort']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

def fanc_to_cools(fanC,outputFormat,output_directory,tmp_directory,thread,sample,resolution):
    resolutions = ','.join(map(str, resolution))
    # 先转成cool格式
    command = ["fanc","to-cooler","-t",str(thread),fanC,os.path.join(tmp_directory, sample+".cool"),'-M']
    try:
        subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    
    # cool -> mcool,h5
    if outputFormat== 'cool':
        src = os.path.join(tmp_directory, sample + ".cool")
        dst = os.path.join(output_directory, os.path.basename(src)) 
        shutil.move(src, dst)
        return 'cool',os.path.join(output_directory, sample+".cool")
    if outputFormat== 'h5':
        command = ["hicConvertFormat","--matrices",os.path.join(tmp_directory, sample+".cool"),"--inputFormat","cool","--outputFormat","h5","-o",os.path.join(output_directory, sample+".h5")]
        try:
            subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
        return 'cool',os.path.join(tmp_directory, sample+".cool")
    if outputFormat== 'mcool': # mcool/juicerhic
        # if (len(resolution) <= 1) and (outputFormat == 'mcool'):
        #   print('Need multiple resolutions ')
        #   sys.exit(1)
        command = ["cooler",'zoomify',os.path.join(tmp_directory, sample+".cool"),"-n",str(thread),"-r",resolutions,"-o",os.path.join(output_directory,sample+'.mcool')]
        try:
            subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
        return 'mcool',os.path.join(output_directory, sample+".mcool")
    # cool -> juicerhic
    if outputFormat== 'juicerhic':
        if len(resolution) <= 1:
            intermediate = 'cool'
            intermediate_cool = os.path.join(tmp_directory, sample+".cool")
            return intermediate,intermediate_cool
        else:
            command = ["cooler",'zoomify',os.path.join(tmp_directory, sample+".cool"),"-n",str(thread),"-r",resolutions,"-o",os.path.join(tmp_directory,sample+'.mcool')]
            try:
                subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                print(f"Error：{e.stderr.decode()}")
            intermediate = 'mcool'
            intermediate_cool = os.path.join(tmp_directory,sample+'.mcool')
    return intermediate,intermediate_cool

def cools_to_juicerhic(input_cool,sample,output_directory,chrom_sizes,juicer_tools_jar,tmp_directory,resolution,thread,input_format):
    output_hic = os.path.join(output_directory,sample+'.hic')
    output_bedpe = os.path.join(tmp_directory,sample+".bedpe")
    output_short = os.path.join(tmp_directory,output_bedpe+".short")
    output_sorted = os.path.join(tmp_directory,output_short+".sorted")
    
    if input_format == 'mcool':
        cooler_paths = cooler.fileops.list_coolers(input_cool)
        hic_resolution = [int(path.split('/')[-1]) for path in cooler_paths]
        highest_res = min(hic_resolution)
        command = ["cooler", "dump", "--join",f"{input_cool}::/resolutions/{highest_res}"] 
    elif input_format == 'cool':
        command = ["cooler", "dump", "--join",input_cool]
    with open(output_bedpe, "w") as fout:
        subprocess.run(command, stdout=fout, check=True)
    
    command = ["awk", "-F", "\t","{print 0, $1, $2, 0, 0, $4, $5, 1, $7}",output_bedpe]
    with open(output_short, "w") as fout:
        subprocess.run(command, stdout=fout, check=True)

    command = ["sort", "-k2,2d", "-k6,6d",output_short]
    with open(output_sorted, "w") as fout:
        subprocess.run(command, stdout=fout, check=True)

    res_list = ','.join([str(i) for i in resolution])
    command = [
        "java","-jar", juicer_tools_jar,
        "pre", "-r", res_list,
        output_sorted, output_hic, chrom_sizes,'-j',str(thread),'--threads',str(thread)
    ] 
    try:
        subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    return output_hic

##########################################################################################################################################
# Parameters
parser = argparse.ArgumentParser(description="HiC format transfer")

# Required parameters
parser.add_argument('--inputFormat', type=str, required=True, help="Input format", choices=['pairs','cool','mcool','juicerhic','fanchic','h5'])
parser.add_argument('--outputFormat', type=str, required=True, help="Output format", choices=['cool','mcool','juicerhic','fanchic','h5'])
parser.add_argument('-i', '--input', type=str, required=True, help="Input file")
parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")
# Optional parameters
parser.add_argument('-t','--thread', type=int, default=10, help="Number of threads, default: 10")
parser.add_argument('-r','--resolution', type=str, default='', help="Resolutions")
parser.add_argument('-n','--norm', type=str, default='None', help="Norm method",choices=['KR','ICE','VC','VC-SQRT','None'])
parser.add_argument('--oe',type=str, default='',help="The resolution of O/E matrix to get")
parser.add_argument('--sparse',type=str, default='',help="The resolution of sparse matrix to get")
parser.add_argument('--species',type=str, default='',help="Set the species to transfer to juicer.hic")

args = parser.parse_args()

# pairs input needs fasta
if args.inputFormat == 'pairs':
    if args.species == '':
        print('No species set')
        sys.exit(1)
    try:
        with open(config_path,'r') as f:
            config = json.load(f)
        fasta = config['data'][args.species]['fasta']
    except Exception as e:
        print(f"Error reading the configuration file: {e}")
        sys.exit(1)
# Need tools to transfer to juicer hic, check
if args.outputFormat == 'juicerhic':
    if args.species == '':
        print("Please set the species to transfer to juicer hic format")
        sys.exit(1)
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        chrom_sizes = config['data'][args.species]['chromInfo']
        juicer_tools_jar = config['software']['juicer_tools_jar']
    except Exception as e:
        print(f"Error reading the configuration file: {e}")
        sys.exit(1)

inputFormat = args.inputFormat
outputFormat = args.outputFormat
thread = int(args.thread)

output_directory = args.output
input_file = args.input
tmp_directory = os.path.join(args.output,'tmp')
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
if not os.path.exists(tmp_directory):
    os.makedirs(tmp_directory)

file_name = os.path.basename(input_file)
last_dot_index = file_name.rfind('.')
sample = file_name[:last_dot_index]

# h5 -> cool, treat as cool file
if inputFormat == 'h5':
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Convert h5 -> cool")
    command = ["hicConvertFormat","--matrices",input_file,"--inputFormat","h5","--outputFormat","cool","-o",os.path.join(tmp_directory, sample+".cool")]
    try:
        subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} h5 -> cool done")
    input_file = os.path.join(tmp_directory, sample+".cool")
    inputFormat = 'cool'

if (inputFormat == 'pairs') and (args.resolution == ''):
    print('Need to specify the resolution if input format is pairs')
    sys.exit(1)
# Get the resolutions
if inputFormat != 'pairs':
    # Get resolutions
    if inputFormat == 'cool':
        hic_sample = fanc.load(input_file)
        hic_resolution = [hic_sample.bin_size] # one resolution
    elif inputFormat == 'mcool':
        cooler_paths = cooler.fileops.list_coolers(input_file)
        hic_resolution = [int(path.split('/')[-1]) for path in cooler_paths] # multiple
    elif inputFormat == 'juicerhic':
        hic_sample = fanc.load(input_file)
        hic_resolution = hic_sample.resolutions()[0] # one / multiple resolutions
    elif inputFormat == 'fanchic':
        hic_sample = fanc.load(input_file)
        hic_resolution = [hic_sample.bin_size] # one resolution

def can_other_numbers_be_divided_by_min(nums):
    if not nums:
        return False  
    min_num = min(nums)
    for num in nums:
        if num != min_num and num % min_num != 0:
            return False
    return True

if args.resolution != '':
    resolutions = args.resolution
    resolution = [int(x) for x in resolutions.split(',')]
    resolution = sorted(resolution)
    resolution_min = min(resolution)
else:
    resolution = hic_resolution
    resolution = sorted(resolution)
    resolution_min = min(resolution)

norm = args.norm
if args.oe != '':
    oe_res = int(args.oe)
if args.sparse != '':
    sparse_res = int(args.sparse)
##########################################################################################################################################
# Check resolutions
if inputFormat == 'pairs':
    hic_resolution = [resolution_min]
if min(hic_resolution) > resolution_min:
    print(f'Can not transfer to higher resolution {resolution_min}, highest available resolution: {min(hic_resolution)}')
    sys.exit(1)
if not can_other_numbers_be_divided_by_min(resolution+hic_resolution):
    print(f'Can not transform the resolutions:{resolution}, the resolution can not be divisible by the highest resolution')
    sys.exit(1)
highest_res = min(hic_resolution) 

if (args.inputFormat == args.outputFormat) and (args.norm == 'None'):
    if hic_resolution == resolution:
        print("No format change / normalization / resolution change is set")   
        sys.exit(1)
if (args.oe != '') and (highest_res > oe_res):
    print(f'OE resolution {oe_res} higher than the highest resolution: {highest_res}')
    sys.exit(1)
if (args.oe != '') and (oe_res not in resolution):
    print(f'OE resolution {oe_res} not in resolution list {resolution}')
    sys.exit(1)
if (args.sparse != '') and (highest_res > sparse_res):
    print(f'Sparse resolution {sparse_res} higher than the highest resolution: {highest_res}')
    sys.exit(1)
if (args.sparse != '') and (sparse_res not in resolution):
    print(f'Sparse resolution {sparse_res} not in resolution list {resolution}')
    sys.exit(1)

if (len(resolution) <= 1) and (outputFormat == 'mcool'):
    print('Need to set multiple resolutions for mcool')
    sys.exit(1)

if (len(resolution) > 1) and ((outputFormat == 'cool') or (outputFormat == 'fanchic') or (outputFormat == 'h5')):
    print('Cool and fanchic and h5 only accept one resolution')
    sys.exit(1)
##########################################################################################################################################
# 1. Transfer to highest resolution, normalized (or not) FanC .hic format first
# pairs/cool/mcool/juicerhic/fanchic -> fanchic

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Convert to fanc format")
if (inputFormat == 'mcool') or (inputFormat == 'juicerhic'):
    command = ["fanc","hic",input_file+'@'+str(highest_res),os.path.join(tmp_directory, sample+"_fanc.hic"),"-t", str(thread), "-f",'--deepcopy','-b',str(resolution_min)]
if (inputFormat == 'cool') or (inputFormat == 'fanchic'):
    command = ["fanc","hic",input_file,os.path.join(tmp_directory, sample+"_fanc.hic"),"-t", str(thread), "-f",'--deepcopy','-b',str(resolution_min)]
if inputFormat == 'pairs':
    command = ["fanc","fragments",fasta,str(highest_res),os.path.join(tmp_directory,"frag.bed")]
    try:
        subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    command = ["fanc","pairs",input_file,os.path.join(tmp_directory,sample+"_fanc.pairs"),"-g",os.path.join(tmp_directory,"frag.bed"),'-t',str(thread),'-f']
    try:
        subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    command = ["fanc","hic",os.path.join(tmp_directory,sample+"_fanc.pairs"),os.path.join(tmp_directory, sample+"_fanc.hic"),"-t", str(thread), "-f",'--deepcopy']

if norm != 'None':
    command = command + ["-n","-m",norm]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} fanc format conversion done")

if outputFormat== 'fanchic':
    src = os.path.join(tmp_directory, sample + "_fanc.hic")
    dst = os.path.join(output_directory, os.path.basename(src))
    shutil.move(src, dst)
    Output = os.path.join(output_directory, sample+"_fanc.hic")
else:
# 2. fanchic -> cool/mcool
    fanC = os.path.join(tmp_directory, sample+"_fanc.hic")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Convert to cool format")
    intermediate,intermediate_cool = fanc_to_cools(fanC,outputFormat,output_directory,tmp_directory,thread,sample,resolution)
    Output = intermediate_cool
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} cool format conversion done")
# 3. cool/mcool -> juicerhic
if outputFormat== 'juicerhic':
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Convert to juicer hic format")
    Output = cools_to_juicerhic(intermediate_cool,sample,output_directory,chrom_sizes,juicer_tools_jar,tmp_directory,resolution,thread,intermediate)
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} juicer hic format conversion done")

# Get the oe matrix if set
if args.oe != '':
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extract O/E matrix")
    if outputFormat in ['mcool','juicerhic']:
        command = ["fanc","dump","-e","--only-intra",Output+'@'+str(oe_res),os.path.join(output_directory,sample+"_oe.txt")]
        try:
            subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
    else:
        command = ["fanc","dump","-e","--only-intra",Output,os.path.join(output_directory,sample+"_oe.txt")]
        try:
            subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} O/E matrix extraction done")

# Get the sparse matrix if set
if args.sparse != '':
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extract sparse contact matrix")
    if outputFormat in ['mcool','juicerhic']:
        command = ["fanc","dump","--only-intra","-u",Output+'@'+str(sparse_res),os.path.join(output_directory,sample+"_sparse.txt")]
        try:
            subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
    else:
        command = ["fanc","dump","--only-intra","-u",Output,os.path.join(output_directory,sample+"_sparse.txt")]
        try:
            subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error：{e.stderr.decode()}")
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extract sparse contact matrix")

print('\nDone, you can check the results now')
