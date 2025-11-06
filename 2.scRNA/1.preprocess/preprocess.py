import os
import argparse
import sys
import subprocess
import resource
import json
import shutil
import pandas as pd

def cellranger(args):
    input_file = args.i
    test = pd.read_csv(input_file)

    # Check columns
    for col in ("sample","fq1","fq2"):
        if col not in test.columns:
            print(f"Error: Can not find column {col} in the input csv")
            sys.exit(1)
    
    for check in ['sample','fq1','fq2']:
        if test.duplicated(check).any():
            dup = test[test.duplicated(check)][check].iloc[0]
            print(f"Find duplications: {dup}")
            sys.exit(1)

    check_cols = ["fq1", "fq2"]
    for idx, row in test.iterrows():
        sample = row["sample"]
        # Check files exist or not
        for col in check_cols:
            raw = row[col]
            src = raw
            if not os.path.exists(src):
                print(f"Can not find {sample} {col}: {src}")
                sys.exit(1)

    with open(os.environ.get('CDesk_config'), "r") as f:
        config = json.load(f)

    if not os.path.exists(config['data'][args.s]['scRef10x']):
        print(f'No {args.s} cellranger database found in config file')
        sys.exit(1)

    os.makedirs(args.o,exist_ok=True)
    os.makedirs(f'{args.o}/cellranger',exist_ok=True)
    os.makedirs(f'{args.o}/cellranger/input',exist_ok=True)

    # Link files
    for idx, row in test.iterrows():
        # fq1
        source_file = os.path.abspath(row['fq1'])
        symlink_name = f"{row['sample']}_S1_L001_R1_001.fastq.gz"
        target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
        os.symlink(source_file, target_file)
        # fq2
        source_file = os.path.abspath(row['fq2'])
        symlink_name = f"{row['sample']}_S1_L001_R2_001.fastq.gz"
        target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
        os.symlink(source_file, target_file)

    # Run cellranger
    for sample in list(test['sample']):
        cmd = [config['software']['cellranger'],"count",f"--id={sample}",f"--fastqs={args.o}/cellranger/input",f"--sample={sample}",f"--transcriptome={config['data'][args.s]['scRef10x']}","--localcores",str(args.t)]
        subprocess.run(cmd, check=True)
        os.system(f'mv {sample} {args.o}/cellranger')

def drseq(args):
    input_file = args.i
    test = pd.read_csv(input_file)

    # Check columns
    for col in ("sample","fq1","fq2"):
        if col not in test.columns:
            print(f"Error: Can not find column {col} in the input csv")
            sys.exit(1)
    
    for check in ['sample','fq1','fq2']:
        if test.duplicated(check).any():
            dup = test[test.duplicated(check)][check].iloc[0]
            print(f"Find duplications: {dup}")
            sys.exit(1)

    check_cols = ["fq1", "fq2"]
    for idx, row in test.iterrows():
        sample = row["sample"]
        # Check files exist or not
        for col in check_cols:
            raw = row[col]
            src = raw
            if not os.path.exists(src):
                print(f"Can not find {sample} {col}: {src}")
                sys.exit(1)
    
    os.makedirs(args.o,exist_ok=True)

    with open(os.environ.get('CDesk_config'), "r") as f:
        config = json.load(f)

    if not os.path.exists(f'{args.o}/DrSeq'):
        os.makedirs(f'{args.o}/DrSeq')
          
    os.path.abspath(row['fq2'])  
    # Run Drseq
    for idx, row in test.iterrows():
        cmd = [config['software']['DrSeq'],"simple","-b",os.path.abspath(row['fq1']),"-r",os.path.abspath(row['fq2']),"-n",row["sample"],"-g",config['data'][args.s]['refgenes'],"--maptool", "bowtie2","--mapindex", config['data'][args.s]['bowtie2_mapindex'],"--thread", str(args.t),"-f","--clean"]
        subprocess.run(cmd, check=True)
        os.system(f'mv {sample} {args.o}/DrSeq')
        os.system(f'mv {sample}.conf {args.o}/DrSeq')

def singleron(args):
    input_file = args.i
    test = pd.read_csv(input_file)

    # Check columns
    for col in ("sample","fq1","fq2"):
        if col not in test.columns:
            print(f"Error: Can not find column {col} in the input csv")
            sys.exit(1)
    
    for check in ['sample','fq1','fq2']:
        if test.duplicated(check).any():
            dup = test[test.duplicated(check)][check].iloc[0]
            print(f"Find duplications: {dup}")
            sys.exit(1)

    check_cols = ["fq1", "fq2"]
    for idx, row in test.iterrows():
        sample = row["sample"]
        # Check files exist or not
        for col in check_cols:
            raw = row[col]
            src = raw
            if not os.path.exists(src):
                print(f"Can not find {sample} {col}: {src}")
                sys.exit(1)

    os.makedirs(args.o, exist_ok=True)

    with open(os.environ.get('CDesk_config'), "r") as f:
        config = json.load(f)

    os.makedirs(f'{args.o}/singleron', exist_ok=True)
    os.makedirs(f'{args.o}/singleron/data', exist_ok=True)

    # Link files
    for idx, row in test.iterrows():
        # fq1
        source_file = os.path.abspath(row['fq1'])
        symlink_name = f"{row['sample']}_1.fq.gz"
        target_file = os.path.join(f'{args.o}/singleron/data',symlink_name)
        os.symlink(source_file, target_file)
        # fq2
        source_file = os.path.abspath(row['fq2'])
        symlink_name = f"{row['sample']}_2.fq.gz"
        target_file = os.path.join(f'{args.o}/singleron/data',symlink_name)
        os.symlink(source_file, target_file)

    filename = f'{args.o}/singleron/rna.mapfile'
    with open(filename, 'w') as f:
        lines = [f'{x}\t{args.o}/singleron/data\t{x}' for x in list(test['sample'])]
        f.write('\n'.join(lines))

    # Run celescope
    cmd = [config['software']['celescope_conda']+'/bin/multi_rna', '--mapfile', f"{args.o}/singleron/rna.mapfile",
       '--genomeDir', config['data'][args.s]['singleron_mapindex'],
       '--thread', str(args.t),
       '--mod', 'shell',
       '--limitBAMsortRAM','100000000000',
       '--outdir', f"{args.o}/singleron"]
    subprocess.run(cmd, check=True)

    if os.path.exists(os.path.join(args.o,'singleron','shell')):
        shutil.rmtree(os.path.join(args.o,'singleron','shell'))

    shutil.move("shell", os.path.join(args.o,'singleron'))

    for file_name in os.listdir(f'{args.o}/singleron/shell'):
        os.environ['PATH'] = f"{os.path.join(config['software']['celescope_conda'],'bin')}:{os.environ['PATH']}"
        cmd = ['bash', f"{args.o}/singleron/shell/{file_name}"]
        subprocess.run(cmd, check=True)

def dnbc(args):
    input_file = args.i
    test = pd.read_csv(input_file)

    # Check columns
    for col in ("sample","cDNAfq1","cDNAfq2","oligofq1","oligofq2"):
        if col not in test.columns:
            print(f"Error: Can not find column {col} in the input csv")
            sys.exit(1)
    
    for check in ["sample","cDNAfq1","cDNAfq2","oligofq1","oligofq2"]:
        if test.duplicated(check).any():
            dup = test[test.duplicated(check)][check].iloc[0]
            print(f"Find duplications: {dup}")
            sys.exit(1)

    check_cols = ["cDNAfq1","cDNAfq2","oligofq1","oligofq2"]
    for idx, row in test.iterrows():
        sample = row["sample"]
        # Check files exist or not
        for col in check_cols:
            raw = row[col]
            src = raw
            if not os.path.exists(src):
                print(f"Can not find {sample} {col}: {src}")
                sys.exit(1)

    with open(os.environ.get('CDesk_config'), "r") as f:
        config = json.load(f)
    os.makedirs(args.o, exist_ok=True)
    os.makedirs(f'{args.o}/dnbc', exist_ok=True)
    os.chdir(f'{args.o}/dnbc')

    for idx, row in test.iterrows():
        cmd = [
            config['software']['dnbc4tools'],'rna','run',
            '--name',row['sample'],
            '--cDNAfastq1',row['cDNAfq1'],'--cDNAfastq2',row['cDNAfq2'],
            '--oligofastq1',row['oligofq1'],'--oligofastq2',row['oligofq2'],
            '--genomeDir',config['data'][args.s]['dnbc_mapindex'],
            '--threads', str(args.t),
            '--outdir',f'{args.o}/dnbc'
        ]
        subprocess.run(cmd, check=True)


def main():
    # CDesk ArgumentParser Object
    parser = argparse.ArgumentParser(prog='scRNA', description="scRNA pipeline")
    subparsers = parser.add_subparsers(dest='type', required=True, help="Choose an omic to run the script")

    # 10x
    # cellranger
    cellranger_parser = subparsers.add_parser('cellranger', help='cellranger 10x preprocess')
    cellranger_parser.add_argument('-i', required=True, help='Specify the input fq information file')
    cellranger_parser.add_argument('-s', required=True, help='Specify the species')
    cellranger_parser.add_argument('-o', required=True, help='Output directory')
    cellranger_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)

    # DrSeq
    drseq_parser = subparsers.add_parser('drseq', help='DrSeq 10x preprocess')
    drseq_parser.add_argument('-i', required=True, help='Specify the input fq information file')
    drseq_parser.add_argument('-s', required=True, help='Specify the species')
    drseq_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    drseq_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)

    # singleron
    singleron_parser = subparsers.add_parser('singleron',help='singleron pipeline')
    singleron_parser.add_argument('-i', required=True, help='Specify the input fq information file')
    singleron_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    singleron_parser.add_argument('-s', required=True, help='Specify the species')
    singleron_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)
    
    # dnbc
    dnbc_parser = subparsers.add_parser('dnbc',help='dnbc pipeline')
    dnbc_parser.add_argument('-i', required=True, help='Specify the input fq information file')
    dnbc_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    dnbc_parser.add_argument('-s', required=True, help='Specify the species')
    dnbc_parser.add_argument('-t', help='Specify number of threads (default is 10)',default=10,type=int)

    # Parameters
    args = parser.parse_args()

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
