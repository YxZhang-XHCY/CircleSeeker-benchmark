#!/bin/bash
# collect_benchmark_results.sh
# 将各工具的结果文件和模拟 truth BED 文件收集到统一的 benchmark 目录中
# 使用方法: cd 到数据根目录后运行 bash scripts/collect_benchmark_results.sh

BASE_DIR="."
OUT_DIR="benchmark_collect"

GENOMES=("ColCEN_5200" "human_10000_simple" "human_23000")
REPS=("rep1" "rep2" "rep3")
DEPTHS=("sequencing_10X" "sequencing_30X" "sequencing_50X")

total=0
skipped=0

for genome in "${GENOMES[@]}"; do
    for rep in "${REPS[@]}"; do
        for depth in "${DEPTHS[@]}"; do
            prefix="${genome}_${rep}_${depth}"
            src="${BASE_DIR}/${genome}/${rep}/${depth}"
            dst="${OUT_DIR}/${genome}/${rep}/${depth}"

            mkdir -p "${dst}"

            # ============================================================
            # 1. Truth BED 文件 (region/ 目录, 所有深度共享同一套 truth)
            # ============================================================
            truth_dir="${BASE_DIR}/${genome}/${rep}/region"
            for bed_type in all unique chimeric multi; do
                truth_file="${truth_dir}/sim_ecc_${rep}.${bed_type}.bed"
                if [ -f "${truth_file}" ]; then
                    cp -n "${truth_file}" "${dst}/truth_${bed_type}.bed"
                    total=$((total + 1))
                fi
            done
            # lengths.tsv 包含 eccDNA 长度分布信息
            lengths_file="${truth_dir}/sim_ecc_${rep}.lengths.tsv"
            if [ -f "${lengths_file}" ]; then
                cp -n "${lengths_file}" "${dst}/truth_lengths.tsv"
                total=$((total + 1))
            fi

            # ============================================================
            # 2. CircleMap Enhanced (NGS 短读长工具)
            #    主结果: *_filtered.bed (过滤后) 和 *_CM.bed (原始)
            # ============================================================
            cm_filtered="${src}/circlemap_output/eccDNA_${prefix}_filtered.bed"
            cm_raw="${src}/circlemap_output/eccDNA_${prefix}_CM.bed"
            if [ -f "${cm_filtered}" ]; then
                cp -n "${cm_filtered}" "${dst}/CircleMap_filtered.bed"
                total=$((total + 1))
            else
                echo "  [MISS] ${cm_filtered}"
                skipped=$((skipped + 1))
            fi
            if [ -f "${cm_raw}" ]; then
                cp -n "${cm_raw}" "${dst}/CircleMap_raw.bed"
                total=$((total + 1))
            fi

            # ============================================================
            # 3. CircleSeeker (HiFi 长读长工具)
            #    主结果: *_eccDNA_summary.csv (汇总)
            #    辅助:   *_eccDNA_regions.csv (区域级别)
            # ============================================================
            cs_summary="${src}/circleseeker_output/${prefix}_eccDNA_summary.csv"
            cs_regions="${src}/circleseeker_output/${prefix}_eccDNA_regions.csv"
            if [ -f "${cs_summary}" ]; then
                cp -n "${cs_summary}" "${dst}/CircleSeeker_summary.csv"
                total=$((total + 1))
            else
                echo "  [MISS] ${cs_summary}"
                skipped=$((skipped + 1))
            fi
            if [ -f "${cs_regions}" ]; then
                cp -n "${cs_regions}" "${dst}/CircleSeeker_regions.csv"
                total=$((total + 1))
            fi

            # ============================================================
            # 4. CReSIL (ONT 长读长工具)
            #    主结果: eccDNA_final.txt
            # ============================================================
            cr_file="${src}/cresil_output/eccDNA_final.txt"
            if [ -f "${cr_file}" ]; then
                cp -n "${cr_file}" "${dst}/CReSIL_eccDNA_final.txt"
                total=$((total + 1))
            else
                echo "  [MISS] ${cr_file}"
                skipped=$((skipped + 1))
            fi

            # ============================================================
            # 5. CReSIL-HiFi (HiFi 长读长工具)
            #    主结果: eccDNA_final.txt (与 CReSIL ONT 格式相同)
            # ============================================================
            cr_hifi_file="${src}/cresil_hifi_output/eccDNA_final.txt"
            if [ -f "${cr_hifi_file}" ]; then
                cp -n "${cr_hifi_file}" "${dst}/CReSIL_HiFi_eccDNA_final.txt"
                total=$((total + 1))
            else
                echo "  [MISS] ${cr_hifi_file}"
                skipped=$((skipped + 1))
            fi

            # ============================================================
            # 6. eccDNA_RCA_nanopore (ONT 长读长工具)
            #    主结果: *_info.tsv (检测信息)
            #    辅助:   *_eccDNA.fa (eccDNA 序列)
            # ============================================================
            rca_info="${src}/eccDNA_RCA_output/${prefix}_info.tsv"
            rca_fa="${src}/eccDNA_RCA_output/${prefix}_eccDNA.fa"
            if [ -f "${rca_info}" ]; then
                cp -n "${rca_info}" "${dst}/eccDNA_RCA_info.tsv"
                total=$((total + 1))
            else
                echo "  [MISS] ${rca_info}"
                skipped=$((skipped + 1))
            fi
            if [ -f "${rca_fa}" ]; then
                cp -n "${rca_fa}" "${dst}/eccDNA_RCA_eccDNA.fa"
                total=$((total + 1))
            fi

            # ============================================================
            # 7. ecc_finder (NGS 短读长工具)
            #    主结果: *.csv (eccDNA 检测结果)
            # ============================================================
            ef_file="${src}/ecc_finder_output/${prefix}.csv"
            if [ -f "${ef_file}" ]; then
                cp -n "${ef_file}" "${dst}/ecc_finder.csv"
                total=$((total + 1))
            else
                echo "  [MISS] ${ef_file}"
                skipped=$((skipped + 1))
            fi

            echo "[OK] ${genome}/${rep}/${depth}"
        done
    done
done

echo ""
echo "========================================="
echo "收集完成！"
echo "  输出目录: ${OUT_DIR}/"
echo "  收集文件数: ${total}"
echo "  缺失文件数: ${skipped}"
echo "========================================="
echo ""
echo "目录结构示例:"
echo "  ${OUT_DIR}/{genome}/{rep}/{depth}/"
echo "    truth_all.bed           <- 全部 eccDNA 区域 (ground truth)"
echo "    truth_unique.bed        <- unique eccDNA 区域"
echo "    truth_chimeric.bed      <- chimeric eccDNA 区域"
echo "    truth_multi.bed         <- multi-mapping eccDNA 区域"
echo "    truth_lengths.tsv       <- eccDNA 长度分布"
echo "    CircleMap_filtered.bed  <- CircleMap Enhanced 结果"
echo "    CircleMap_raw.bed       <- CircleMap Enhanced 原始结果"
echo "    CircleSeeker_summary.csv <- CircleSeeker 汇总结果"
echo "    CircleSeeker_regions.csv <- CircleSeeker 区域结果"
echo "    CReSIL_eccDNA_final.txt <- CReSIL (ONT) 结果"
echo "    CReSIL_HiFi_eccDNA_final.txt <- CReSIL-HiFi 结果"
echo "    eccDNA_RCA_info.tsv     <- eccDNA_RCA_nanopore 结果"
echo "    eccDNA_RCA_eccDNA.fa    <- eccDNA_RCA_nanopore 序列"
echo "    ecc_finder.csv          <- ecc_finder 结果"
