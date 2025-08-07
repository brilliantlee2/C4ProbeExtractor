import os
import subprocess
import pysam
import gzip
import shutil
import sys
import logging
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--fastq", help="FASTQ of Cyclone reads", type=str)
    parser.add_argument("--probe", help="FASTA of Cyclone probe squence(s)", type=str)
    # Optional arguments
    parser.add_argument(
        "--output_fastq",
        help="Output file name for stranded FASTQ entries",
        type=str,
        default="output.fastq.gz",
    )
    parser.add_argument(
        "-t", "--threads", help="Threads to use [8]", type=int, default=8)
    
    args = parser.parse_args()

    return args
    


def run_subprocess(cmd):
    """
    Run OS command and return stdout & stderr
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return str(stdout), str(stderr)

def check_bwa():
    if shutil.which("bwa") is None:
        logging.error("Could not find 'bwa' executable in PATH. Please install it.")
        sys.exit(1)
    else:
        print("bwa 检测通过")

def bwa_index(reference):
    """构建bwa index"""
    if not all(os.path.exists(reference + ext) for ext in [".bwt", ".pac", ".ann", ".amb", ".sa"]):
        print("正在构建 BWA 索引...")
        bwa_idx_cmd = "bwa index {ref}".format(ref=reference)
        stdout, stderr = run_subprocess(bwa_idx_cmd)
    else:
        print("索引文件已存在，跳过 index")

def bwa_mem(reference, fastq, output_bam=None, threads=8):
    """运行bwa mem比对并排序输出为bam"""
    print("正在运行 BWA MEM 比对...")
    if output_bam is None:
        # 获取 basename（去掉路径），然后最多去掉两个扩展名（.gz 和 .fastq）
        base = os.path.basename(fastq)
        if base.endswith(".gz"):
            base = os.path.splitext(base)[0]  # 去掉 .gz
        base = os.path.splitext(base)[0]      # 再去掉 .fastq 或 .fq
        output_bam = base + ".sorted.bam"
    
    bwa_mem_cmd = "bwa mem -k 5 -B 1 -O 1,1 -E 1,1 -L 0,0 -T 10 -a -t {threads} {ref} {fq} | samtools sort -o {output_bam}".format(
    ref=reference, fq=fastq,output_bam=output_bam,threads=threads
    )

    stdout, stderr = run_subprocess(bwa_mem_cmd)
    print(f"比对完成，输出文件为：{output_bam}")
    return output_bam

def get_mapped_read_names(bam_path):
    """获取bam中比对成功的read名"""
    mapped_reads = set()
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if not read.is_unmapped:
                mapped_reads.add(read.query_name)
    return mapped_reads

def extract_reads_from_fastq(fastq_path, read_names_set, output_path):
    """从fastq中提取比对到的read"""
    open_func = gzip.open if fastq_path.endswith(".gz") else open
    with open_func(fastq_path, "rt") as fin,  gzip.open(output_path, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            read_name = header.strip().split()[0][1:]  # 去掉开头的 @

            if read_name in read_names_set:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
def main(agrs):
    #os.environ["PATH"] = "/home/liyy/anaconda3/envs/Py310/bin:" + os.environ["PATH"]
    check_bwa()
    # ====== 参数 ======
    reference = args.probe
    fastq = args.fastq
    extracted_fastq = args.output_fastq
    threads = args.threads

    # ====== 执行流程 ======
    bwa_index(reference)
    bam_file = bwa_mem(reference=reference, fastq=fastq, threads=threads)

    print("提取 BAM 中比对成功的 reads...")
    mapped_reads = get_mapped_read_names(bam_file)
    print(f"比对成功 reads 数量: {len(mapped_reads)}")

    print("提取 reads 到 FASTQ...")
    extract_reads_from_fastq(fastq, mapped_reads, extracted_fastq)
    print(f"输出提取结果：{extracted_fastq}")

if __name__ == "__main__":
    args = parse_args()
    main(args)