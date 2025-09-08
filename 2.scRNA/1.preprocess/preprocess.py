import os
import argparse
import sys
import subprocess
import resource
import json

def cellranger(args):
    files = os.listdir(args.i)
    fastqgz_samples = []; fq_samples = []

    for file in files:
        if file.endswith('_1.fastq.gz'):
            prefix = file.replace('_1.fastq.gz', '')
            if os.path.exists(os.path.join(args.i, prefix+'_2.fastq.gz')):
                fastqgz_samples.append(prefix)
            else:
                print(f'{prefix}_1.fastq.gz has no corresponding {prefix}_2.fastq.gz')
        if file.endswith('_1.fq.gz'):
            prefix = file.replace('_1.fq.gz', '')
            if os.path.exists(os.path.join(args.i, prefix+'_2.fq.gz')):
                fq_samples.append(prefix)
            else:
                print(f'{prefix}_1.fq.gz has no corresponding {prefix}_2.fq.gz')

    if not fastqgz_samples and not fq_samples:
        print('No fastq.gz or fq.gz files found')
        sys.exit(1)
    
    os.makedirs(args.o,exist_ok=True)
    os.makedirs(f'{args.o}/cellranger',exist_ok=True)
    os.makedirs(f'{args.o}/cellranger/input',exist_ok=True)

    with open(args.config, "r") as f:
        config = json.load(f)

    # cellranger
    # transcriptome = {    
    #     "pig": "/mnt/liudong/data/Genome/susScr11Self",
    #     "mm10": "/mnt/liudong/data/Genome/mm10Self",
    #     "human": "/mnt/linzejie/CDesk/data/cellranger/human" 
    # }

    if not os.path.exists(config['data'][args.s]['scRef10x']):
        print(f'No {args.s} cellranger database found in config file')
        sys.exit(1)

    # 软链接，修改文件名
    for sample in fastqgz_samples:
        symlink_name = f"{sample}_S1_L001_R1_001.fastq.gz"
        source_file = os.path.join(args.i,sample+'_1.fastq.gz')
        target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
        os.symlink(source_file, target_file)

        symlink_name = f"{sample}_S1_L001_R2_001.fastq.gz"
        source_file = os.path.join(args.i,sample+'_2.fastq.gz')
        target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
        os.symlink(source_file, target_file)

    for sample in fq_samples:
        symlink_name = f"{sample}_S1_L001_R1_001.fastq.gz"
        source_file = os.path.join(args.i,sample+'_1.fq.gz')
        target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
        os.symlink(source_file, target_file)

        symlink_name = f"{sample}_S1_L001_R2_001.fastq.gz"
        source_file = os.path.join(args.i,sample+'_2.fq.gz')
        target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
        os.symlink(source_file, target_file)

    samples = fastqgz_samples + fq_samples
    # 运行cellranger
    for sample in samples:
        #cmd = (f"{config['software']['cellranger']} count --id={sample} --fastqs={args.o}/cellranger/input --sample={sample} --transcriptome={config['data'][args.s]['scRef10x']} --localcores {args.t}")
        cmd = [config['software']['cellranger'],"count",f"--id={sample}",f"--fastqs={args.o}/cellranger/input",f"--sample={sample}",f"--transcriptome={config['data'][args.s]['scRef10x']}","--localcores",str(args.t)]
        subprocess.run(cmd, check=True)
        os.system(f'mv {sample} {args.o}/cellranger')

def drseq(args):
    files = os.listdir(args.i)
    fastq_samples = []

    for file in files:
        if file.endswith('_1.fastq'):
            prefix = file.replace('_1.fastq', '')
            if os.path.exists(os.path.join(args.i, prefix+'_2.fastq')):
                fastq_samples.append(prefix)
            else:
                print(f'{prefix}_1.fastq has no corresponding {prefix}_2.fastq')

    if not fastq_samples:
        print('No fastq files found')
        sys.exit(1)
    
    os.makedirs(args.o,exist_ok=True)

    # DrSeq
    # refgene = {    
    #     "pig": "/mnt/liudong/data/Genome/susScr11/susScr11.DrSeq.refgenes.txt"
    # }
    # mapindex = {
    #     "pig":"/mnt/linzejie/CDesk/data/drseq/susScr11_bowtie2/susScr11"
    # }

    with open(args.config, "r") as f:
        config = json.load(f)

    if not os.path.exists(f'{args.o}/DrSeq'):
        os.makedirs(f'{args.o}/DrSeq')

        
    # 运行DrSeq
    for sample in fastq_samples:
        #cmd = (f"{config['software']['DrSeq']} simple -b {args.i}/{sample}_1.fastq -r {args.i}/{sample}_2.fastq -n {sample} -g {config['data'][args.s]['refgenes']} --maptool bowtie2 --mapindex {config['data'][args.s]['mapping_index']} --thread {args.t} -f")
        cmd = [config['software']['DrSeq'],"simple","-b", f"{args.i}/{sample}_1.fastq","-r", f"{args.i}/{sample}_2.fastq","-n", sample,"-g", config['data'][args.s]['refgenes'],"--maptool", "bowtie2","--mapindex", config['data'][args.s]['bowtie2_mapindex'],"--thread", str(args.t),"-f","--clean"]
        subprocess.run(cmd, check=True)
        os.system(f'mv {sample} {args.o}/DrSeq')
        os.system(f'mv {sample}.conf {args.o}/DrSeq')

def singleron(args):
    #command = f"conda run -n celescope python -c 'import sys; print(sys.executable)'"
    #soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
    #resource.setrlimit(resource.RLIMIT_NOFILE, (1000000, 1000000))
    #try:
    #    subprocess.run(command, shell=True, check=True)
    #    print("激活celescope conda 环境成功")
    #except Exception as e:
    #    print(str(e))
    #    print('激活celescope conda环境失败，请检查')
    #    sys.exit(1)

    files = os.listdir(args.i)
    fastqgz_samples = []; fq_samples = []
    for file in files:
        if file.endswith('_1.fastq.gz'):
            prefix = file.replace('_1.fastq.gz', '')
            if os.path.exists(os.path.join(args.i, prefix+'_2.fastq.gz')):
                fastqgz_samples.append(prefix)
            else:
                print(f'{prefix}_1.fastq.gz has no corresponding {prefix}_2.fastq.gz')
        if file.endswith('_1.fq.gz'):
            prefix = file.replace('_1.fq.gz', '')
            if os.path.exists(os.path.join(args.i, prefix+'_2.fq.gz')):
                fq_samples.append(prefix)
            else:
                print(f'{prefix}_1.fq.gz has no corresponding {prefix}_2.fq.gz')
    
    if not fastqgz_samples and not fq_samples:
        print('No fastq.gz or fq.gz files found')
        sys.exit(1)

    os.makedirs(args.o, exist_ok=True)

    with open(args.config, "r") as f:
        config = json.load(f)
    # singleron
    # mapindex = {    
    #     "galGal": "/mnt/linzejie/temp/galGal6_celescope"
    # }
    #if not os.path.exists(f'{args.o}/singleron'):
    os.makedirs(f'{args.o}/singleron', exist_ok=True)

    samples = fastqgz_samples + fq_samples

    filename = f'{args.o}/singleron/rna.mapfile'
    with open(filename, 'w') as f:
        lines = [f'{x}\t{args.i}\t{x}' for x in samples]
        f.write('\n'.join(lines))

    #cmd = (f"multi_rna --mapfile {args.o}/singleron/rna.mapfile --genomeDir {mapindex[args.s]} --thread {args.t} --mod shell --outdir {args.o}/singleron")
    #os.system(cmd)
    cmd = [config['conda_env']['celescope']+'/bin/multi_rna', '--mapfile', f"{args.o}/singleron/rna.mapfile",
       '--genomeDir', config['data'][args.s]['singleron_mapindex'],
       '--thread', str(args.t),
       '--mod', 'shell',
       '--outdir', f"{args.o}/singleron"]
    subprocess.run(cmd, check=True)
    os.system(f'mv shell {args.o}/singleron')

    for file_name in os.listdir(f'{args.o}/singleron/shell'):
        #base_name = file_name.replace('_.sh', '')
        #cmd = (f'bash {args.o}/singleron/shell/{file_name}')
        #cmd = ['bash', f"{args.o}/singleron/shell/{file_name}"]
        os.environ['PATH'] = f"{config['conda_env']['celescope']+'/bin'}:{os.environ['PATH']}"
        cmd = ['bash', f"{args.o}/singleron/shell/{file_name}"]
        subprocess.run(cmd, check=True)

def dnbc(args):
    with open(args.config, "r") as f:
        config = json.load(f)
    os.makedirs(args.o, exist_ok=True)
    os.makedirs(f'{args.o}/dnbc', exist_ok=True)
    args.i = os.path.abspath(args.i)
    args.o = os.path.abspath(args.o)
    os.chdir(f'{args.o}/dnbc')
    # maps = {
    #     "pig": "/mnt/linzejie/CDesk_data/data/scRNA/dnbc4tools/pig",
    #     "mm10": "/mnt/linzejie/CDesk_data/data/scRNA/dnbc4tools/mm10",
    #     "human": "/mnt/linzejie/CDesk_data/data/scRNA/dnbc4tools/hg38"
    # }
    cmd = [
    config['software']['dnbc4tools'],'rna','multi',
    '--list',args.i, # 一个文件定义
    '--genomeDir',config['data'][args.s]['dnbc_mapindex'],
    '--threads', str(args.t),
    '--outdir',f'{args.o}/dnbc'
    ]
    subprocess.run(cmd, check=True)
    sh_files = [f for f in os.listdir(f'{args.o}/dnbc') if f.endswith('.sh')]

    for sh_file in sh_files:
        full_path = os.path.join(f'{args.o}/dnbc',sh_file) 
        # 获取文件前缀（去掉扩展名）
        file_prefix = os.path.splitext(sh_file)[0]
        cmd = ['bash',full_path]
        subprocess.run(cmd, check=True)
        print(f"{file_prefix} finished \n")

def main():
    # CDesk ArgumentParser Object
    parser = argparse.ArgumentParser(prog='scRNA', description="scRNA pipeline")
    # 创建二级子命令
    subparsers = parser.add_subparsers(dest='type', required=True, help="Choose an omic to run the script")

    # 10x
    # 创建cellranger子命令的子解析器
    cellranger_parser = subparsers.add_parser('cellranger', help='cellranger 10x preprocess')
    cellranger_parser.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq.gz, X_2.fastq.gz')
    cellranger_parser.add_argument('-s', required=True, help='Specify the species')
    cellranger_parser.add_argument('-o', required=True, help='Output directory')
    cellranger_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)
    cellranger_parser.add_argument('--config',required=True, help='The config file')

    # 创建DrSeq子命令的子解析器
    drseq_parser = subparsers.add_parser('drseq', help='DrSeq 10x preprocess')
    drseq_parser.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq, X_2.fastq')
    drseq_parser.add_argument('-s', required=True, help='Specify the species')
    drseq_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    drseq_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)
    drseq_parser.add_argument('--config',required=True, help='The config file')

    # singleron
    singleron_parser = subparsers.add_parser('singleron',help='singleron pipeline')
    singleron_parser.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq.gz, X_2.fastq.gz')
    singleron_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    singleron_parser.add_argument('-s', required=True, help='Specify the species')
    singleron_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)
    singleron_parser.add_argument('--config',required=True, help='The config file')
    
    # dnbc
    dnbc_parser = subparsers.add_parser('dnbc',help='dnbc pipeline')
    dnbc_parser.add_argument('-i', required=True, help='Specify the input tsv file')
    dnbc_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    dnbc_parser.add_argument('-s', required=True, help='Specify the species')
    dnbc_parser.add_argument('-t', help='Specify number of threads (default is 10)',default=10,type=int)
    dnbc_parser.add_argument('--config',required=True, help='The config file')

    # 解析命令行参数
    args = parser.parse_args()

    # 判断并执行 bulkRNA 相关子命令
    if args.type == 'cellranger':
        cellranger(args)
    elif args.type == 'drseq':
        drseq(args)
    elif args.type == 'singleron':
        singleron(args)
    elif args.type=='dnbc':
        dnbc(args)

if __name__ == "__main__":
    main()

