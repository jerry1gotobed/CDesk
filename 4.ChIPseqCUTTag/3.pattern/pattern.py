import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Run ChIPseqCUTTag specific pattern analysis")

parser.add_argument("-bw",required=True, type=str, help="The input bw files list file")
parser.add_argument("-bed",required=True, type=str, help="The bed genomic regions to be analyzed list file")
parser.add_argument("-o",required=True,type=str, help="Output directory")

parser.add_argument("--mode",default='reference-point',help="Analysis mode",choices=["reference-point","scale-regions"])

parser.add_argument("--height",help="Heatmap height",type=int,default=10)
parser.add_argument("-a",help="Length to extend upstream of the genomic region start site (bp)",type=int,default=2000)
parser.add_argument("-b",help="Length to extend downstream of the genomic region end site (bp)",type=int,default=2000)
parser.add_argument("-m",help="Scaled genomic region length (bp)",type=int,default=3000)
parser.add_argument("-t",help="Number of processors",type=int,default=10)
parser.add_argument("--region",help="The reference point for reference-point",default='center',choices=['TSS','TES','center'])

args = parser.parse_args()

os.makedirs(args.o, exist_ok=True)

with open(args.bw, 'r') as file:
    bws = [line.strip() for line in file.readlines()]

with open(args.bed, 'r') as file:
    beds = [line.strip() for line in file.readlines()]

# Step 1: computeMatrix scale-regions
compute_matrix_cmd = [
    "computeMatrix", args.mode,
    "-a", str(args.a),
    "-b", str(args.b),
    "-p",str(args.t),
    "--skipZeros","--missingDataAsZero",#'--binSize','5',
    "-o", os.path.join(args.o,"matrix.mat.gz")] + ['-S'] + bws + ['-R'] + beds

if args.mode == "scale-regions":
    compute_matrix_cmd = compute_matrix_cmd + ['-m',str(args.m)]
else:
    compute_matrix_cmd = compute_matrix_cmd + ['--referencePoint',args.region]

subprocess.run(compute_matrix_cmd, check=True)

# Step 2: plotHeatmap
plot_heatmap_cmd = [
    "plotHeatmap",'--heatmapHeight',str(args.height),
    "-m", os.path.join(args.o,"matrix.mat.gz"),
    "-out", os.path.join(args.o,"Result.pdf"),
    "--whatToShow", "plot, heatmap and colorbar",
    #"--colorList","blue,white,red",
    "--colorMap","Reds","--missingDataColor","white",
    "--legendLocation","best"
]
subprocess.run(plot_heatmap_cmd, check=True)

print('Finished, you can see the results now')
