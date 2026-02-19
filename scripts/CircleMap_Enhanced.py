import argparse
import os
import subprocess
import sys
import pandas as pd
from pathlib import Path
from typing import List

class PipelineError(Exception):
    """Custom exception for pipeline errors"""
    pass

def run_command(command, log_file=None):
    """执行命令并处理可能的错误"""
    try:
        if log_file:
            with open(log_file, 'w') as f:
                result = subprocess.run(command, shell=True, check=True,
                                     stdout=f, stderr=subprocess.STDOUT)
        else:
            result = subprocess.run(command, shell=True, check=True,
                                  capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}")
        print(f"Error message: {e.stderr if e.stderr else 'No error message'}")
        return False

def check_output_exists(*files):
    """检查所有指定的输出文件是否存在"""
    return all(os.path.exists(f) for f in files)

def round_float_columns(df, columns, decimals):
    """Round float columns to a specific number of decimals"""
    for col in columns:
        if col in df.columns:
            df[col] = df[col].round(decimals)
    return df

def filter_circlemap_bed(input_bed: str, output_bed: str, sample_name: str) -> bool:
    """
    读取Circle-Map输出的原始bed文件，应用过滤条件，输出过滤后的bed文件

    CircleMap 实际输出的列名:
    chrom, start, end, discordants, soft-clipped, score, mean, std, start_ratio, end_ratio, continuity
    """
    # CircleMap 实际输出的列名（文件带表头）
    float_columns = ['score', 'mean', 'std', 'start_ratio', 'end_ratio', 'continuity']

    try:
        # 读取原始bed文件（带表头）
        df = pd.read_csv(input_bed, sep='\t')

        if df.empty:
            print(f"Warning: Input bed file is empty: {input_bed}")
            open(output_bed, 'w').close()
            return True

        # 计算Length
        df['length'] = df['end'] - df['start'] + 1

        # 添加name列
        df['name'] = df.apply(lambda row: f"{row['chrom']}-{row['start']}-{row['end']}", axis=1)

        # 添加sample列
        df['sample'] = sample_name

        # Round float columns
        df = round_float_columns(df, float_columns, 2)

        # 应用过滤条件（使用实际的列名）
        filtered_df = df[
            (df['score'] > 50) &
            (df['soft-clipped'] >= 1) &
            (df['discordants'] >= 1) &
            (df['start_ratio'] > 0.33) &
            (df['end_ratio'] > 0.33) &
            (df['length'] < 1e7)
        ]

        print(f"Filtering results: {len(df):,} raw -> {len(filtered_df):,} filtered eccDNAs")

        # 重排列顺序：name, sample, length放前面，其余保持原顺序
        original_cols = ['chrom', 'start', 'end', 'discordants', 'soft-clipped',
                        'score', 'mean', 'std', 'start_ratio', 'end_ratio', 'continuity']
        primary_cols = ['name', 'sample', 'length']
        final_cols = primary_cols + [col for col in original_cols if col in filtered_df.columns]
        filtered_df = filtered_df[final_cols]

        # 保存过滤后的bed文件（带表头）
        filtered_df.to_csv(output_bed, sep='\t', header=True, index=False)
        print(f"Filtered BED saved: {output_bed}")

        return True

    except Exception as e:
        print(f"Error filtering bed file: {e}")
        import traceback
        traceback.print_exc()
        return False

class Pipeline:
    def __init__(self, args):
        self.threads = args.threads
        self.fastq1 = args.fastq1
        self.fastq2 = args.fastq2
        self.reference = args.reference
        self.output_dir = args.output_dir
        self.sample_name = args.sample_name
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        os.chdir(self.output_dir)

    def run_fastp(self):
        """运行fastp步骤"""
        output_files = [
            f"{self.sample_name}.R1.fp.fastq.gz",
            f"{self.sample_name}.R2.fp.fastq.gz",
            f"{self.sample_name}.report.html",
            f"{self.sample_name}.report.json"
        ]
        
        if check_output_exists(*output_files):
            print("Fastp outputs already exist, skipping...")
            return True
            
        cmd = (f"fastp -w 16 -i {self.fastq1} -I {self.fastq2} "
               f"-o {self.sample_name}.R1.fp.fastq.gz "
               f"-O {self.sample_name}.R2.fp.fastq.gz "
               f"-h {self.sample_name}.report.html "
               f"-j {self.sample_name}.report.json")
        
        return run_command(cmd, f"fastp_{self.sample_name}.log")

    def run_bwa(self):
        """运行BWA MEM步骤"""
        output_file = f"{self.sample_name}.sam"
        
        if check_output_exists(output_file):
            print("BWA outputs already exist, skipping...")
            return True
            
        cmd = (f"bwa mem -q -t {self.threads} {self.reference} "
               f"{self.sample_name}.R1.fp.fastq.gz "
               f"{self.sample_name}.R2.fp.fastq.gz > {output_file}")
        
        return run_command(cmd)

    def run_samtools_sort_qname(self):
        """运行samtools按查询名排序"""
        output_file = f"qname_eccDNA_{self.sample_name}.bam"
        
        if check_output_exists(output_file):
            print("Samtools qname sort outputs already exist, skipping...")
            return True
            
        cmd = (f"samtools sort -@ {self.threads} -n "
               f"-o {output_file} {self.sample_name}.sam")
        
        return run_command(cmd)

    def run_samtools_sort_coordinate(self):
        """运行samtools按坐标排序"""
        output_files = [
            f"sorted_eccDNA_{self.sample_name}.bam",
            f"sorted_eccDNA_{self.sample_name}.bam.bai"
        ]
        
        if check_output_exists(*output_files):
            print("Samtools coordinate sort outputs already exist, skipping...")
            return True
            
        cmd1 = (f"samtools sort -@ {self.threads} "
                f"-o sorted_eccDNA_{self.sample_name}.bam {self.sample_name}.sam")
        cmd2 = (f"samtools index -@ {self.threads} "
                f"sorted_eccDNA_{self.sample_name}.bam")
        
        return run_command(cmd1) and run_command(cmd2)

    def run_circle_map_extractor(self):
        """运行Circle-Map ReadExtractor"""
        output_file = f"eccDNA_{self.sample_name}_candidates.bam"
        
        if check_output_exists(output_file):
            print("Circle-Map ReadExtractor outputs already exist, skipping...")
            return True
            
        cmd = (f"Circle-Map ReadExtractor -i qname_eccDNA_{self.sample_name}.bam "
               f"-o {output_file}")
        
        return run_command(cmd, f"circle_map_read_extractor_{self.sample_name}.log")

    def run_samtools_sort_candidates(self):
        """对candidates进行排序"""
        output_files = [
            f"sort_eccDNA_{self.sample_name}_candidates.bam",
            f"sort_eccDNA_{self.sample_name}_candidates.bam.bai"
        ]
        
        if check_output_exists(*output_files):
            print("Samtools candidates sort outputs already exist, skipping...")
            return True
            
        cmd1 = (f"samtools sort -@ {self.threads} "
                f"eccDNA_{self.sample_name}_candidates.bam "
                f"-o sort_eccDNA_{self.sample_name}_candidates.bam")
        cmd2 = (f"samtools index -@ {self.threads} "
                f"sort_eccDNA_{self.sample_name}_candidates.bam")
        
        return run_command(cmd1) and run_command(cmd2)

    def run_circle_map_realign(self):
        """运行Circle-Map++ Realign，输出原始bed和过滤后的bed"""
        raw_bed = f"eccDNA_{self.sample_name}_CM.bed"
        filtered_bed = f"eccDNA_{self.sample_name}_filtered.bed"
        log_file = f"circle_map_realign_{self.sample_name}.log"
        
        # 如果两个输出文件都存在，跳过
        if check_output_exists(raw_bed, filtered_bed):
            print("Circle-Map++ Realign outputs already exist, skipping...")
            return True
        
        # Step 1: 运行 Circle-Map++ Realign (同步执行，不用nohup)
        cmd = (f"circle_map++ Realign -t {self.threads} "
               f"-i sort_eccDNA_{self.sample_name}_candidates.bam "
               f"-qbam qname_eccDNA_{self.sample_name}.bam "
               f"-sbam sorted_eccDNA_{self.sample_name}.bam "
               f"-fasta {self.reference} "
               f"-o {raw_bed}")
        
        print(f"Running Circle-Map++ Realign...")
        if not run_command(cmd, log_file):
            print("Error executing Circle-Map++ Realign")
            return False
        
        print(f"Raw BED saved: {raw_bed}")
        
        # Step 2: 过滤原始bed文件，生成filtered bed
        if not filter_circlemap_bed(raw_bed, filtered_bed, self.sample_name):
            print("Error filtering bed file")
            return False
        
        return True

    def run_pipeline(self):
        """运行完整的流程"""
        steps = [
            (self.run_fastp, "Fastp"),
            (self.run_bwa, "BWA MEM"),
            (self.run_samtools_sort_qname, "Samtools Sort by Query Name"),
            (self.run_samtools_sort_coordinate, "Samtools Sort by Coordinate"),
            (self.run_circle_map_extractor, "Circle-Map ReadExtractor"),
            (self.run_samtools_sort_candidates, "Samtools Sort Candidates"),
            (self.run_circle_map_realign, "Circle-Map++ Realign")
        ]

        print(f"Starting pipeline for sample {self.sample_name}")
        print("="*60)
        
        for step_func, step_name in steps:
            print(f"\n[Step] {step_name}...")
            if not step_func():
                print(f"[Error] {step_name} step failed")
                return False
            print(f"[Done] {step_name} completed")
                
        print("\n" + "="*60)
        print(f"Pipeline completed successfully for sample {self.sample_name}")
        print(f"Output files:")
        print(f"  - Raw BED: eccDNA_{self.sample_name}_CM.bed")
        print(f"  - Filtered BED: eccDNA_{self.sample_name}_filtered.bed")
        return True

def main():
    parser = argparse.ArgumentParser(
        description="Circle-seq Pipeline with Circle-Map++ for eccDNA detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python circle_seq_pipeline_v2.py -t 8 \\
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \\
    -r reference.fa -o ./output -s sample_name

Output files:
  - eccDNA_{sample}_CM.bed        : Raw Circle-Map++ output
  - eccDNA_{sample}_filtered.bed  : Filtered eccDNAs (score>50, splits>2, etc.)
        """
    )
    parser.add_argument("-t", "--threads", required=True, type=int, help="Number of threads")
    parser.add_argument("-1", "--fastq1", required=True, help="Input fastq1 file")
    parser.add_argument("-2", "--fastq2", required=True, help="Input fastq2 file")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")

    args = parser.parse_args()

    # 验证输入文件是否存在
    for file_path in [args.fastq1, args.fastq2, args.reference]:
        if not os.path.exists(file_path):
            print(f"Error: Input file {file_path} does not exist")
            sys.exit(1)

    pipeline = Pipeline(args)
    success = pipeline.run_pipeline()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

