#!/bin/bash
#SBATCH --job-name=CReSIL_HiFi
#SBATCH --output=logs_cresil_hifi/%x_%A_%a.out
#SBATCH --error=logs_cresil_hifi/%x_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --array=0-26

# CReSIL-HiFi 批量运行脚本 (SLURM array job)
# 用于 PacBio HiFi 长读长数据的 eccDNA 检测
# 基于 CReSIL V1.2.0+hifi (添加 HiFi 平台支持)

eval "$(conda shell.bash hook)"
conda activate cresil-hifi

THREADS=16
# === 用户配置 ===
# 请修改为实际的数据根目录
BASE_DIR="/path/to/eccDNA_simulation"
REF_DIR="${BASE_DIR}/Ref"

# 日志目录
LOG_DIR="${BASE_DIR}/logs_cresil_hifi"
mkdir -p ${LOG_DIR}

# 汇总统计文件
STATS_FILE="${LOG_DIR}/cresil_hifi_benchmark_stats.tsv"
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    echo -e "Sample\tStatus\tWall_Time(s)\tMax_Memory(KB)\tMax_Memory(GB)\tStart_Time\tEnd_Time" > ${STATS_FILE}
fi

# 构建任务列表: 3 datasets x 3 reps x 3 depths = 27 tasks
DATASETS=("ColCEN_5200" "human_23000" "human_10000_simple")
REPS=("rep1" "rep2" "rep3")
DEPTHS=("sequencing_10X" "sequencing_30X" "sequencing_50X")

# 从 array index 计算对应的 dataset/rep/depth
IDX=${SLURM_ARRAY_TASK_ID}
D_IDX=$((IDX / 9))
R_IDX=$(( (IDX % 9) / 3 ))
P_IDX=$((IDX % 3))

DATASET=${DATASETS[$D_IDX]}
REP=${REPS[$R_IDX]}
DEPTH=${DEPTHS[$P_IDX]}

# 设置参考基因组
if [ "$DATASET" == "ColCEN_5200" ]; then
    REF_FA="${REF_DIR}/ColCEN.fasta"
else
    REF_FA="${REF_DIR}/chm13v2.0.fa"
fi

REF_FAI="${REF_FA}.fai"
REF_MMI="${REF_FA}.hifi.mmi"

# 生成 fasta 索引（加锁防止并发）
if [ ! -f "$REF_FAI" ]; then
    (
        flock -x 200
        if [ ! -f "$REF_FAI" ]; then
            echo "[INFO] 生成 fasta 索引: $REF_FAI"
            samtools faidx "$REF_FA"
        fi
    ) 200>${REF_FA}.fai.lock
fi

# 生成 minimap2 HiFi 索引（加锁防止并发）
if [ ! -f "$REF_MMI" ]; then
    (
        flock -x 200
        if [ ! -f "$REF_MMI" ]; then
            echo "[INFO] 生成 minimap2 HiFi 索引: $REF_MMI"
            minimap2 -x map-hifi -d "$REF_MMI" "$REF_FA"
        fi
    ) 200>${REF_MMI}.lock
fi

# 构建路径 (HiFi 数据使用 .HiFi.fasta 文件)
INPUT_DIR="${BASE_DIR}/${DATASET}/${REP}/${DEPTH}"
INPUT_FILE="${INPUT_DIR}/sim_ecc_${REP}.HiFi.fasta"
OUTPUT_DIR="${INPUT_DIR}/cresil_hifi_output"
PREFIX="${DATASET}_${REP}_${DEPTH}"
LOG_FILE="${LOG_DIR}/${PREFIX}.log"
TIME_FILE="${LOG_DIR}/${PREFIX}.time"

echo "=== CReSIL-HiFi Task ${SLURM_ARRAY_TASK_ID} ==="
echo "Sample: ${PREFIX}"
echo "Input:  ${INPUT_FILE}"
echo "Ref:    ${REF_FA}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# 检查输入文件
if [ ! -f "$INPUT_FILE" ]; then
    echo "[跳过] 输入文件不存在: $INPUT_FILE"
    exit 0
fi

# 检查是否已完成
if [ -f "${OUTPUT_DIR}/eccDNA_final.txt" ]; then
    echo "[已完成] ${PREFIX}"
    exit 0
fi

START_TIME=$(date "+%Y/%m/%d %H:%M")

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 运行 CReSIL-HiFi (使用 -platform hifi 参数)
/usr/bin/time -v -o "${TIME_FILE}" bash -c "
    cd ${OUTPUT_DIR}

    # Step 1: Trim - 预处理 reads (HiFi 模式)
    echo '[Step 1/2] Running cresil trim (HiFi mode)...'
    cresil trim -t ${THREADS} \
        -platform hifi \
        -fq ${INPUT_FILE} \
        -r ${REF_MMI} \
        -o ${OUTPUT_DIR}

    # Step 2: Identify - 识别 eccDNA (HiFi 模式, 跳过 Medaka 校正)
    echo '[Step 2/2] Running cresil identify (HiFi mode)...'
    cresil identify -t ${THREADS} \
        -platform hifi \
        -fa ${REF_FA} \
        -fai ${REF_FAI} \
        -fq ${INPUT_FILE} \
        -trim ${OUTPUT_DIR}/trim.txt
" > "${LOG_FILE}" 2>&1

EXIT_CODE=$?
END_TIME=$(date "+%Y/%m/%d %H:%M")

# 提取时间和内存信息
WALL_TIME=$(grep "Elapsed (wall clock)" "${TIME_FILE}" | awk '{print $NF}')
MAX_MEM_KB=$(grep "Maximum resident set size" "${TIME_FILE}" | awk '{print $NF}')

if [ -n "$MAX_MEM_KB" ]; then
    MAX_MEM_GB=$(echo "scale=2; ${MAX_MEM_KB} / 1024 / 1024" | bc)
else
    MAX_MEM_KB="N/A"
    MAX_MEM_GB="N/A"
fi

if [ -n "$WALL_TIME" ]; then
    if echo "$WALL_TIME" | grep -q ":"; then
        WALL_SEC=$(echo "$WALL_TIME" | awk -F: '{if(NF==3){print $1*3600+$2*60+$3}else{print $1*60+$2}}')
    else
        WALL_SEC=$WALL_TIME
    fi
else
    WALL_SEC="N/A"
fi

if [ $EXIT_CODE -eq 0 ] && [ -f "${OUTPUT_DIR}/eccDNA_final.txt" ]; then
    STATUS="完成"
    echo "[完成] ${PREFIX} - 用时: ${WALL_TIME}, 内存: ${MAX_MEM_GB} GB"
else
    STATUS="失败"
    echo "[失败] ${PREFIX} - 查看日志: ${LOG_FILE}"
fi

# 写入统计文件（加锁避免并发写入冲突）
(
    flock -x 200
    echo -e "${PREFIX}\t${STATUS}\t${WALL_SEC}\t${MAX_MEM_KB}\t${MAX_MEM_GB}\t${START_TIME}\t${END_TIME}" >> ${STATS_FILE}
) 200>${LOG_DIR}/.stats.lock

echo "=== 任务结束 ==="
