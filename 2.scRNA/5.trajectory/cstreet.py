import argparse
import warnings
import subprocess
import os
import json
warnings.filterwarnings('ignore')

def flatten_list(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            # 如果元素是列表，递归调用
            flat_list.extend(flatten_list(item))
        else:
            # 如果元素不是列表，直接添加到结果列表
            flat_list.append(item)
    return flat_list


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run CStreet analysis")
    
    parser.add_argument("-i",required=True,type=str, help="Input expression matrixes list file")
    parser.add_argument("-s",required=True,type=str, help="Input cell state list file")
    parser.add_argument("-n",required=True, type=str,default='',help="Output directory name")
    parser.add_argument("-o",required=True, type=str, help="Output directory")
    parser.add_argument("--config",required=True, type=str, help="The config file")
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()
    
    with open(args.i, 'r') as file:
        expMatrix = [line.strip() for line in file if line.strip()]

    with open(args.s, 'r') as file:
        cellState = [line.strip() for line in file if line.strip()]

    os.makedirs(args.o, exist_ok=True)
    
    with open(args.config, "r") as f:
        config = json.load(f)

    cmd = [config['software']['CStreet'],'-i',expMatrix,'-s',cellState,'-n',args.n,'-o',args.o]
    cmd = flatten_list(cmd)
    subprocess.run(cmd, check=True)


if __name__ == "__main__":

    main()
