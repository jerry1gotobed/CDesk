import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import pandas as pd

def process_sample(sample_path, sample_name,output_dir, genome,mode,start,end):
    """处理单个样本"""
    print(f'Processing {sample_name} ......')
    try:
        # 创建样本专属目录
        sample_output_dir = os.path.join(output_dir, f"{sample_name}_motifDir")
        os.makedirs(sample_output_dir, exist_ok=True)

        # 运行 findMotifsGenome.pl（阻塞直到完成）
        motif_log = os.path.join(output_dir, f"{sample_name}_motif.log")

        if mode == 'peak':
            # 生成临时文件路径
            tmp_peak_file = os.path.join(output_dir, f"{sample_name}.homer_peaks.tmp")
            # 转换 BED 文件格式
            with open(tmp_peak_file, 'w') as f:
                subprocess.run(
                    ["awk", '{print $4"\t"$1"\t"$2"\t"$3"\t+"}', sample_path],
                    stdout=f,
                    check=True
                )

            with open(motif_log, 'w') as log:
                subprocess.run(
                    [
                        "findMotifsGenome.pl",
                        tmp_peak_file,
                        genome,
                        sample_output_dir,
                        "-len", "8,10,12",
                        "-preparsedDir", output_dir
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True
                )

            # 注释 peaks
            annotate_xls = os.path.join(output_dir, f"{sample_name}.peakAnn.xls")
            annotate_log = os.path.join(output_dir, f"{sample_name}.annLog.txt")
            subprocess.run(
                f"annotatePeaks.pl {tmp_peak_file} {genome} > {annotate_xls} 2> {annotate_log}",
                shell=True,  # 需要 shell 处理重定向
                check=True
            )
            os.remove(tmp_peak_file)
            
        elif mode == 'gene':
            with open(motif_log, 'w') as log:
                subprocess.run(
                    [
                        "findMotifs.pl",
                        sample_path,
                        genome,
                        sample_output_dir,
                        "-start", str(start),
                        "-end", str(end)
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True
                )



        # 压缩结果目录
        output_basename = os.path.basename(sample_output_dir)
        tar_file = f"{output_basename}.tar.gz"
        subprocess.run(
            [
                "tar", "-czf",
                os.path.join(output_dir, tar_file),
                "-C", output_dir,
                output_basename
            ],
            check=True
        )

        # 清理临时文件
        return f"✅ {sample_name} done"

    except subprocess.CalledProcessError as e:
        return f"❌ {sample_name} 失败：{e}"
    except Exception as e:
        return f"❌ {sample_name} 异常：{str(e)}"

def main():
    input = sys.argv[1]
    output_dir = sys.argv[2]
    genome = sys.argv[3]
    max_workers = int(sys.argv[4])
    mode = sys.argv[5]
    start = int(sys.argv[6])
    end = int(sys.argv[7])

    os.makedirs(output_dir, exist_ok=True)
    if mode == 'peak':
        input = pd.read_csv(input)
        if 'peak' not in input.columns or 'sample' not in input.columns:
            print('Need columns: peak,sample')
            sys.exit(1)
        files = list(input['peak'])
        sample_names = list(input['sample'])
    
    if mode == 'gene':
        input = pd.read_csv(input)
        if 'gene' not in input.columns or 'sample' not in input.columns:
            print('Need columns: gene,sample')
            sys.exit(1)
        files = list(input['gene'])
        sample_names = list(input['sample'])

    # 并发执行任务
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process_sample,
                sample_path,
                sample_name,
                output_dir,
                genome,
                mode,
                start,
                end
                ): (sample_path, sample_name) 
            for sample_path, sample_name in zip(files, sample_names)
        }

        # 跟踪进度
        for future in as_completed(futures):
            result = future.result()
            print(result)

    print("Done, you can see the result now")

if __name__ == "__main__":
    main()

