import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import json

def process_sample(sample_path, input_dir, output_dir, genome,config):
    """处理单个样本"""
    try:
        # 提取样本名
        base_name = os.path.basename(sample_path)
        sample_name = base_name.split('_peaks.bed')[0]

        # 创建样本专属目录
        sample_output_dir = os.path.join(output_dir, f"{sample_name}_motifDir")
        os.makedirs(sample_output_dir, exist_ok=True)

        # 生成临时文件路径
        tmp_peak_file = os.path.join(output_dir, f"{sample_name}.homer_peaks.tmp")

        # 转换 BED 文件格式
        with open(tmp_peak_file, 'w') as f:
            subprocess.run(
                ["awk", '{print $4"\t"$1"\t"$2"\t"$3"\t+"}', sample_path],
                stdout=f,
                check=True
            )

        # 运行 findMotifsGenome.pl（阻塞直到完成）
        motif_log = os.path.join(output_dir, f"{sample_name}_motif.log")
        with open(motif_log, 'w') as log:
            subprocess.run(
                [
                    #"findMotifsGenome.pl",
                    config['software']['findMotifsGenome'],
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
            f"{config['software']['annotatePeaks']} {tmp_peak_file} {genome} > {annotate_xls} 2> {annotate_log}",
            shell=True,  # 需要 shell 处理重定向
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
        os.remove(tmp_peak_file)
        return f"✅ {sample_name} 完成"

    except subprocess.CalledProcessError as e:
        return f"❌ {sample_name} 失败：{e}"
    except Exception as e:
        return f"❌ {sample_name} 异常：{str(e)}"

def main():
    if len(sys.argv) != 7:
        print("用法：python script.py <输入目录> <输出目录> <基因组> <最大并发数>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    genome = sys.argv[3]
    max_workers = int(sys.argv[4])
    peak_summit = sys.argv[5] #_peaks.bed
    config_path = sys.argv[6]
    
    with open(config_path, "r") as f:
        config = json.load(f)

    os.makedirs(output_dir, exist_ok=True)
    # 获取所有 BED 文件
    bed_files = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(peak_summit)
    ]

    # 并发执行任务
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process_sample,
                sample_path,
                input_dir,
                output_dir,
                genome,
                config
            ): sample_path
            for sample_path in bed_files
        }

        # 跟踪进度
        for future in as_completed(futures):
            result = future.result()
            print(result)

    print("所有样本处理完成！")

if __name__ == "__main__":
    main()
