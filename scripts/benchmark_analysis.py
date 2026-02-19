#!/usr/bin/env python3
"""
eccDNA Benchmark Analysis
比较 CircleMap Enhanced, CircleSeeker, CReSIL, eccDNA_RCA_nanopore 四个工具
与模拟 truth BED 的检测性能

匹配规则: reciprocal overlap >= 90% (Uecc/Mecc), similarity >= 90% (Cecc ALL fragments)
评估维度: Overall + 按类型 (Uecc/Mecc/Cecc)

用法: python benchmark_analysis.py [benchmark_collect_dir]
"""

import os
import re
import csv
import sys
import statistics
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Tuple, Set


# ============================================================
# 配置
# ============================================================

GENOMES = ["ColCEN_5200", "human_10000_simple", "human_23000"]
REPS = ["rep1", "rep2", "rep3"]
DEPTHS = ["sequencing_10X", "sequencing_30X", "sequencing_50X"]
TOOLS = ["CircleMap", "CircleSeeker", "CReSIL", "eccDNA_RCA"]
OVERLAP_THRESHOLD = 0.9


# ============================================================
# 数据结构
# ============================================================

@dataclass
class Region:
    chrom: str
    start: int
    end: int
    name: str = ""

    @property
    def length(self):
        return self.end - self.start

    def overlap(self, other):
        if self.chrom != other.chrom:
            return 0
        return max(0, min(self.end, other.end) - max(self.start, other.start))

    def reciprocal_overlap(self, other, threshold=0.5):
        ovl = self.overlap(other)
        if ovl == 0:
            return False
        return (ovl / self.length >= threshold) and (ovl / other.length >= threshold)

    def similarity(self, other):
        """ovl / min(len1, len2), 用于 Cecc 片段匹配"""
        ovl = self.overlap(other)
        if ovl == 0:
            return 0.0
        return ovl / min(self.length, other.length)


# ============================================================
# Truth 解析
# ============================================================

def parse_truth(filepath):
    """
    解析 truth BED 文件
    返回 {name: {'type': 'U'/'M'/'C', 'fragments': [Region, ...]}}
    对 U/M 类型: fragments 只含 primary region
    对 C 类型: fragments 包含所有片段 (从第 8 列解析)
    """
    if not os.path.exists(filepath):
        return {}
    truth = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.strip().split('\t')
            name = p[3]
            ecc_type = p[6] if len(p) > 6 else 'U'
            primary = Region(p[0], int(p[1]), int(p[2]), name)

            fragments = [primary]
            if len(p) > 7 and p[7] != '.':
                fragments = []
                for s in p[7].split(';'):
                    m = re.match(r'([^:]+):(\d+)-(\d+)', s)
                    if m:
                        fragments.append(Region(
                            m.group(1), int(m.group(2)), int(m.group(3)), name
                        ))

            truth[name] = {'type': ecc_type, 'fragments': fragments}
    return truth


# ============================================================
# 工具输出解析
# ============================================================

def parse_circlemap(filepath):
    """
    解析 CircleMap_filtered.bed (TSV, 有 header)
    列: name, sample, length, chrom, start, end, ...
    """
    if not os.path.exists(filepath):
        return []
    regions = []
    with open(filepath) as f:
        f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue
            p = line.strip().split('\t')
            regions.append(Region(p[3], int(p[4]), int(p[5]), p[0]))
    return regions


def parse_circleseeker(filepath):
    """
    解析 CircleSeeker_summary.csv
    - Uecc: 直接取 chr, start, end
    - Mecc: 从 location 解析所有 mapping 位点 (用 | 分隔)
    - Cecc: 从 location 解析所有片段 (用 ; 分隔)
    """
    if not os.path.exists(filepath):
        return []
    regions = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ecc_id = row['eccDNA_id']
            ecc_type = row['type']

            if ecc_type == 'Cecc':
                loc = row.get('location', '')
                for frag in loc.split(';'):
                    m = re.match(r'([^:]+):(\d+)-(\d+)', frag)
                    if m:
                        regions.append(Region(
                            m.group(1), int(m.group(2)), int(m.group(3)), ecc_id
                        ))
            elif ecc_type == 'Mecc':
                loc = row.get('location', '')
                for site in loc.split('|'):
                    m = re.match(r'([^:]+):(\d+)-(\d+)', site)
                    if m:
                        regions.append(Region(
                            m.group(1), int(m.group(2)), int(m.group(3)), ecc_id
                        ))
            else:  # Uecc
                regions.append(Region(
                    row['chr'], int(row['start']), int(row['end']), ecc_id
                ))
    return regions


def parse_cresil(filepath):
    """
    解析 CReSIL_eccDNA_final.txt (TSV, 有 header)
    merge_region 格式:
      单区域: Chr3:11369279-11370899_-
      多区域: Chr1:6693790-6697982_-,Chr5:5428789-5429254_-
    """
    if not os.path.exists(filepath):
        return []
    regions = []
    with open(filepath) as f:
        f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue
            p = line.strip().split('\t')
            ecc_id = p[0]
            for region_str in p[1].split(','):
                m = re.match(r'([^:]+):(\d+)-(\d+)', region_str)
                if m:
                    regions.append(Region(
                        m.group(1), int(m.group(2)), int(m.group(3)), ecc_id
                    ))
    return regions


ECCDNA_RCA_DEDUP_THRESHOLD = 0.99


def parse_eccdna_rca(filepath):
    """
    解析 eccDNA_RCA_info.tsv (TSV, 有 header, read 级别)
    1. 筛选 Nfullpass > 0 的 reads
    2. 提取所有 fragments
    3. 按 99% reciprocal overlap 去冗余 (保留不同 eccDNA 边界)
    4. 记录冗余度统计

    返回: (deduped_regions, stats_dict)
    """
    if not os.path.exists(filepath):
        return [], {}

    all_frags = []
    total_reads_detected = 0

    with open(filepath) as f:
        f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue
            p = line.strip().split('\t')
            nfullpass = int(p[1])
            if nfullpass == 0:
                continue
            total_reads_detected += 1
            if len(p) <= 5 or not p[5]:
                continue
            for frag in p[5].split('|'):
                m = re.match(r'([^:]+):(\d+)-(\d+)', frag)
                if m:
                    all_frags.append((m.group(1), int(m.group(2)), int(m.group(3))))

    # 排序后按 99% reciprocal overlap 去冗余
    all_frags.sort()
    deduped = []
    threshold = ECCDNA_RCA_DEDUP_THRESHOLD

    for chrom, start, end in all_frags:
        length = end - start
        if length == 0:
            continue
        matched = False
        # 反向遍历已有 clusters, 找 99% reciprocal overlap 匹配
        for i in range(len(deduped) - 1, -1, -1):
            d_chrom, d_start, d_end = deduped[i]
            if d_chrom != chrom:
                break
            if start - d_start > 100000:  # 超出可能范围, 提前终止
                break
            ovl = min(end, d_end) - max(start, d_start)
            if ovl <= 0:
                continue
            d_length = d_end - d_start
            if ovl / length >= threshold and ovl / d_length >= threshold:
                matched = True
                break
        if not matched:
            deduped.append((chrom, start, end))

    regions = [Region(c, s, e, f"rca_{i}") for i, (c, s, e) in enumerate(deduped)]

    stats = {
        'total_reads_detected': total_reads_detected,
        'total_fragments': len(all_frags),
        'deduped_regions': len(deduped),
        'redundancy': round(len(all_frags) / max(len(deduped), 1), 2)
    }
    return regions, stats


# ============================================================
# 评估
# ============================================================

def evaluate_unified(truth_all, detected_regions, threshold=OVERLAP_THRESHOLD):
    """
    统一评估: 一次匹配, 按类型拆分结果。

    匹配顺序:
    1. Cecc (严格: ALL fragments 必须匹配, similarity >= threshold)
    2. Uecc / Mecc (reciprocal overlap >= threshold, 任一 fragment 匹配即可)

    返回: {
        'Overall': {Truth, Detected, TP, FP, FN, Precision, Recall, F1},
        'Uecc':    {Truth, TP, FN, Recall},   (如有)
        'Mecc':    {Truth, TP, FN, Recall},   (如有)
        'Cecc':    {Truth, TP, FN, Recall},   (如有)
    }
    Per-type 只输出 Recall: Precision/FP 在类型级别无意义
    (无法确定 FP 检出属于哪个类型)
    """
    truth_u = {k: v for k, v in truth_all.items() if v['type'] == 'U'}
    truth_m = {k: v for k, v in truth_all.items() if v['type'] == 'M'}
    truth_c = {k: v for k, v in truth_all.items() if v['type'] == 'C'}

    # 按 eccDNA ID 分组
    det_groups = defaultdict(list)
    for r in detected_regions:
        det_groups[r.name].append(r)

    # 染色体索引
    det_by_chr = defaultdict(list)
    for r in detected_regions:
        det_by_chr[r.chrom].append(r)

    matched_truth = {}   # truth_name -> type_code
    matched_groups = set()

    # Phase 1: Cecc 严格匹配 (所有片段都必须匹配)
    if truth_c:
        for det_name, det_frags in det_groups.items():
            if det_name in matched_groups:
                continue
            for truth_name, tinfo in truth_c.items():
                if truth_name in matched_truth:
                    continue
                if match_cecc_fragments(tinfo['fragments'], det_frags, threshold):
                    matched_truth[truth_name] = 'C'
                    matched_groups.add(det_name)
                    break

    # Phase 2: Uecc / Mecc (reciprocal overlap)
    for tname, tinfo in truth_all.items():
        if tname in matched_truth:
            continue
        if tinfo['type'] == 'C':
            continue  # Cecc 已在 Phase 1 处理
        found = False
        for frag in tinfo['fragments']:
            for dr in det_by_chr.get(frag.chrom, []):
                if dr.name in matched_groups:
                    continue
                if frag.reciprocal_overlap(dr, threshold):
                    matched_truth[tname] = tinfo['type']
                    matched_groups.add(dr.name)
                    found = True
                    break
            if found:
                break

    # ---- Overall ----
    n_det = len(det_groups)
    n_tp = len(matched_truth)
    n_fp = n_det - len(matched_groups)
    n_fn = len(truth_all) - n_tp
    prec = n_tp / (n_tp + n_fp) if (n_tp + n_fp) > 0 else 0
    rec = n_tp / (n_tp + n_fn) if (n_tp + n_fn) > 0 else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0

    results = {
        'Overall': {
            'Truth': len(truth_all), 'Detected': n_det,
            'TP': n_tp, 'FP': n_fp, 'FN': n_fn,
            'Precision': round(prec, 4), 'Recall': round(rec, 4),
            'F1': round(f1, 4)
        }
    }

    # ---- Per-type: 只输出 Recall ----
    for code, name, truth_sub in [
        ('U', 'Uecc', truth_u), ('M', 'Mecc', truth_m), ('C', 'Cecc', truth_c)
    ]:
        if not truth_sub:
            continue
        tp_t = sum(1 for tc in matched_truth.values() if tc == code)
        fn_t = len(truth_sub) - tp_t
        rec_t = tp_t / len(truth_sub)
        results[name] = {
            'Truth': len(truth_sub), 'Detected': '',
            'TP': tp_t, 'FP': '', 'FN': fn_t,
            'Precision': '', 'Recall': round(rec_t, 4),
            'F1': ''
        }

    return results


def match_cecc_fragments(truth_frags, detected_frags, similarity_threshold=0.9):
    """
    Cecc 匹配: 所有 truth 片段都必须在 detected 中找到对应片段
    使用 similarity (ovl / min(len1, len2)) 而非 reciprocal_overlap
    """
    if not truth_frags or not detected_frags:
        return False
    truth_matched = 0
    used_detected = set()
    for truth_frag in truth_frags:
        for i, det_frag in enumerate(detected_frags):
            if i in used_detected:
                continue
            if truth_frag.similarity(det_frag) >= similarity_threshold:
                truth_matched += 1
                used_detected.add(i)
                break
    return truth_matched == len(truth_frags)




# ============================================================
# 汇总输出
# ============================================================

def print_summary(all_results):
    """逐行展示每个 Genome × EvalType × Depth × Rep × Tool 的结果"""

    for genome in GENOMES:
        print(f"\n{'=' * 130}")
        print(f"  Genome: {genome}")
        print(f"{'=' * 130}")

        for eval_type in ['Overall', 'Uecc', 'Mecc', 'Cecc']:
            subset = [r for r in all_results
                      if r['Genome'] == genome and r['EvalType'] == eval_type]
            if not subset:
                continue

            is_overall = (eval_type == 'Overall')
            print(f"\n  [{eval_type}]")

            if is_overall:
                has_redundancy = any(r['Redundancy'] != '' for r in subset)
                header = (f"  {'Depth':<6} {'Rep':<5} {'Tool':<14} {'Truth':>6} "
                          f"{'Detected':>9} {'TP':>6} {'FP':>6} {'FN':>6} "
                          f"{'Precision':>10} {'Recall':>10} {'F1':>10}")
                if has_redundancy:
                    header += f" {'Redund.':>8}"
            else:
                has_redundancy = False
                header = (f"  {'Depth':<6} {'Rep':<5} {'Tool':<14} {'Truth':>6} "
                          f"{'TP':>6} {'FN':>6} {'Recall':>10}")
            print(header)
            print(f"  {'-' * (len(header) - 2)}")

            for depth in DEPTHS:
                depth_short = depth.replace('sequencing_', '')
                for rep in REPS:
                    for tool in TOOLS:
                        row = [r for r in subset
                               if r['Depth'] == depth and r['Rep'] == rep
                               and r['Tool'] == tool]
                        if not row:
                            continue
                        r = row[0]
                        if is_overall:
                            line = (
                                f"  {depth_short:<6} {rep:<5} {tool:<14} "
                                f"{r['Truth']:>6} {r['Detected']:>9} "
                                f"{r['TP']:>6} {r['FP']:>6} {r['FN']:>6} "
                                f"{r['Precision']:>10.4f} {r['Recall']:>10.4f} "
                                f"{r['F1']:>10.4f}"
                            )
                            if has_redundancy:
                                if r['Redundancy'] != '':
                                    line += f" {r['Redundancy']:>8}"
                                else:
                                    line += f" {'':>8}"
                        else:
                            line = (
                                f"  {depth_short:<6} {rep:<5} {tool:<14} "
                                f"{r['Truth']:>6} {r['TP']:>6} {r['FN']:>6} "
                                f"{r['Recall']:>10.4f}"
                            )
                        print(line)
                # depth 之间空一行
                if depth != DEPTHS[-1]:
                    print()


def print_rca_stats(rca_stats_all):
    """输出 eccDNA_RCA 冗余度统计 (逐个 rep)"""
    print(f"\n{'=' * 100}")
    print("  eccDNA_RCA 冗余度统计")
    print(f"{'=' * 100}")
    print(f"  {'Genome':<22} {'Depth':<6} {'Rep':<5} {'Detected Reads':>15} "
          f"{'Fragments':>11} {'Deduped':>8} {'Redundancy':>11}")
    print(f"  {'-' * 84}")

    for s in sorted(rca_stats_all, key=lambda x: (x['genome'], x['depth'], x['rep'])):
        depth_short = s['depth'].replace('sequencing_', '')
        print(
            f"  {s['genome']:<22} {depth_short:<6} {s['rep']:<5} "
            f"{s['total_reads_detected']:>15} {s['total_fragments']:>11} "
            f"{s['deduped_regions']:>8} {s['redundancy']:>11.2f}"
        )


# ============================================================
# 主流程
# ============================================================

def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "benchmark_collect"

    if not os.path.isdir(base_dir):
        print(f"错误: 目录不存在 {base_dir}")
        sys.exit(1)

    output_csv = os.path.join(os.path.dirname(base_dir) or '.', "benchmark_results.csv")
    rca_csv = os.path.join(os.path.dirname(base_dir) or '.', "benchmark_rca_stats.csv")

    all_results = []
    rca_stats_all = []

    for genome in GENOMES:
        for rep in REPS:
            for depth in DEPTHS:
                d = os.path.join(base_dir, genome, rep, depth)
                if not os.path.isdir(d):
                    print(f"  [SKIP] {d}")
                    continue

                # ---- 解析 truth ----
                truth_all = parse_truth(os.path.join(d, 'truth_all.bed'))

                # ---- 解析各工具 ----
                tools_detected = {
                    'CircleMap': parse_circlemap(
                        os.path.join(d, 'CircleMap_filtered.bed')),
                    'CircleSeeker': parse_circleseeker(
                        os.path.join(d, 'CircleSeeker_summary.csv')),
                    'CReSIL': parse_cresil(
                        os.path.join(d, 'CReSIL_eccDNA_final.txt')),
                }

                rca_regions, rca_stats = parse_eccdna_rca(
                    os.path.join(d, 'eccDNA_RCA_info.tsv'))
                tools_detected['eccDNA_RCA'] = rca_regions
                rca_stats.update({
                    'genome': genome, 'rep': rep, 'depth': depth
                })
                rca_stats_all.append(rca_stats)

                # 构建 tool → redundancy 映射 (仅 eccDNA_RCA 有值)
                rca_redundancy = rca_stats.get('redundancy', '')

                # ---- 评估 (统一匹配, 按类型拆分) ----
                for tool, detected in tools_detected.items():
                    redundancy = rca_redundancy if tool == 'eccDNA_RCA' else ''
                    unified = evaluate_unified(truth_all, detected)
                    for eval_type, res in unified.items():
                        all_results.append({
                            'Genome': genome, 'Rep': rep, 'Depth': depth,
                            'Tool': tool, 'EvalType': eval_type,
                            'Redundancy': redundancy, **res
                        })

                print(f"  [OK] {genome}/{rep}/{depth}")

    # ---- 保存 CSV ----
    fieldnames = [
        'Genome', 'Rep', 'Depth', 'Tool', 'EvalType',
        'Truth', 'Detected', 'TP', 'FP', 'FN',
        'Precision', 'Recall', 'F1', 'Redundancy'
    ]
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    rca_fieldnames = [
        'genome', 'rep', 'depth',
        'total_reads_detected', 'total_fragments', 'deduped_regions', 'redundancy'
    ]
    with open(rca_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rca_fieldnames)
        writer.writeheader()
        writer.writerows(rca_stats_all)

    # ---- 输出汇总 ----
    print_summary(all_results)
    print_rca_stats(rca_stats_all)

    print(f"\n{'=' * 90}")
    print(f"  详细结果已保存: {output_csv}")
    print(f"  eccDNA_RCA 冗余度: {rca_csv}")
    print(f"{'=' * 90}")

    # ---- 注释 ----
    print("""
注:
  - 统一匹配: 一次匹配所有 truth (Cecc 优先严格匹配, 再匹配 Uecc/Mecc)
  - Overall: Precision / Recall / F1 在全部类型上计算
  - Per-type (Uecc/Mecc/Cecc): 只报告 Recall (TP/FN)
    Precision 在类型级别无意义 (无法确定 FP 检出属于哪个类型)
  - 匹配规则: reciprocal overlap >= 90% (Uecc/Mecc)
  - Cecc: ALL 片段数量+位置都匹配才计为 TP (similarity >= 90%)
  - eccDNA_RCA: 原始结果经 99% reciprocal overlap 去冗余后评估
""")


if __name__ == '__main__':
    main()
