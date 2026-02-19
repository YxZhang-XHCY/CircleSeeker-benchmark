#!/bin/bash
#SBATCH --job-name=eccDNA_RCA
#SBATCH --output=logs_eccDNA_RCA/%x_%A_%a.out
#SBATCH --error=logs_eccDNA_RCA/%x_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --array=0-26

# eccDNA_RCA_nanopore 批量运行脚本 (SLURM array job)
# 用于 ONT 长读长数据的 eccDNA 检测

eval "$(conda shell.bash hook)"
conda activate circleseeker

THREADS=16
# === 用户配置 ===
# 请修改为实际的数据根目录
BASE_DIR="/path/to/eccDNA_simulation"
REF_DIR="${BASE_DIR}/Ref"
# 请修改为 eccDNA_RCA_nanopore.py 的实际路径
TOOL_PATH="/path/to/eccDNA_RCA_nanopore/eccDNA_RCA_nanopore.py"

# 日志目录
LOG_DIR="${BASE_DIR}/logs_eccDNA_RCA"
mkdir -p ${LOG_DIR}

# 汇总统计文件
STATS_FILE="${LOG_DIR}/eccDNA_RCA_benchmark_stats.tsv"
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

# 构建路径
INPUT_DIR="${BASE_DIR}/${DATASET}/${REP}/${DEPTH}"
INPUT_FILE="${INPUT_DIR}/sim_ecc_${REP}.ONT.fastq"
OUTPUT_DIR="${INPUT_DIR}/eccDNA_RCA_output"
PREFIX="${DATASET}_${REP}_${DEPTH}"
LOG_FILE="${LOG_DIR}/${PREFIX}.log"
TIME_FILE="${LOG_DIR}/${PREFIX}.time"

echo "=== eccDNA_RCA_nanopore Task ${SLURM_ARRAY_TASK_ID} ==="
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
if [ -f "${OUTPUT_DIR}/${PREFIX}_info.tsv" ]; then
    echo "[已完成] ${PREFIX}"
    exit 0
fi

START_TIME=$(date "+%Y/%m/%d %H:%M")

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 运行 eccDNA_RCA_nanopore
/usr/bin/time -v -o "${TIME_FILE}" bash -c "
    cd ${OUTPUT_DIR}

    # Step 1: Minimap2 比对生成 PAF 文件
    echo '[Step 1/2] Running minimap2...'
    minimap2 -t ${THREADS} -c ${REF_FA} ${INPUT_FILE} > ${PREFIX}.paf

    # Step 2: 运行 eccDNA_RCA_nanopore 检测
    echo '[Step 2/2] Running eccDNA_RCA_nanopore...'
    python ${TOOL_PATH} \
        --fastq ${INPUT_FILE} \
        --paf ${OUTPUT_DIR}/${PREFIX}.paf \
        --info ${OUTPUT_DIR}/${PREFIX}_info.tsv \
        --seq ${OUTPUT_DIR}/${PREFIX}_eccDNA.fa \
        --var ${OUTPUT_DIR}/${PREFIX}_var.tsv \
        --reference ${REF_FA} \
        --verbose
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

if [ $EXIT_CODE -eq 0 ] && [ -f "${OUTPUT_DIR}/${PREFIX}_info.tsv" ]; then
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
