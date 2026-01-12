import argparse
import subprocess
import os
import tempfile
import sys
import json
from xvfbwrapper import Xvfb

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run IGV snapshots visualization")
    parser.add_argument("-i",required=True, type=str, help="The input list file")
    parser.add_argument("-o",required=True,type=str, help="Output directory")
    parser.add_argument("--genome", required=True, type=str, help="Select the genome id")
    parser.add_argument("--display_mode",default='collapse',help="Display mode", choices=["collapse","expand","squish"])
    parser.add_argument("--region",required=True,help="The region file")
    parser.add_argument("--type", default='png', type=str,choices=['svg','png'] ,help="svg or png output")
    return parser.parse_args()

def main():
    args = parse_arguments()
    os.makedirs(args.o, exist_ok=True)

    igv_script = f"""
    new
    genome {args.genome}
    snapshotDirectory {args.o}
    """

    with open(args.i, 'r') as file:
        files = [line.strip() for line in file if line.strip()]

    for input_file in files:
        input_file = input_file.strip()
        if os.path.exists(input_file):
            if input_file.endswith('.bam'):
                bam_index_file = input_file + ".bai"
                if not os.path.exists(bam_index_file):
                    print('Can not find corresponding bai file, please check')
                    sys.exit(1)
            igv_script += f"load {input_file}\n    "
        else:
            print(f"{input_file} not exist")
            sys.exit(1)

    with open(args.region, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]
    regions = ','.join(lines)  

    for region in regions.split(','):
        region_save = region.replace(':','_')
        region_save = region_save.replace('-','_')
        igv_script += f"""
    goto {region}
    {args.display_mode}
    sort position
    snapshot {region_save+'.'+args.type}
    """
    igv_script += """
    exit
    """
    
    script_file = os.path.join(f'{args.o}/igv_batch_script.txt')
    try:
        with open(script_file, "w") as f:
            f.write(igv_script)
    except Exception as e:
        print(f"Error writing IGV script to file: {e}")
        return

    igv_cmd = [
        'igv',
        "-b", script_file
    ]

    with Xvfb(width=1024, height=768, colordepth=24):
        try:
            subprocess.run(igv_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Execute IGV wrong: {e.returncode}")

    print("Snapshot generation completed successfully.")

if __name__ == "__main__":
    main()
