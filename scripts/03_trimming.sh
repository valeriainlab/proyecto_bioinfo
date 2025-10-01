#!/usr/bin/env bash
set -euo pipefail
IN="data/raw_data"
OUTD="data/trimmed_data"
QCD="qc/fastp"
mkdir -p qc/fastp


# Illumina paired-end

fastp \
  -i  "${IN}/anc_R1.fastq" -I "${IN}/anc_R2.fastq" \
  -o "${OUTD}/anc_R1.trimmed.fastq" \
  -O "${OUTD}/anc_R2.trimmed.fastq" \
  --trim_front2 0 \
  --trim_tail2 0 \
  --detect_adapter_for_pe \
  --qualified_quality_phred 28 \
  --length_required 80 \
  --trim_front1 0 --trim_tail1 0\
  --cut_front --cut_tail --cut_mean_quality 28 \
  --thread 4 \
  --html "${QCD}/fastp_anc.html" \


fastp \
  -i  "${IN}/evol1_R1.fastq" -I "${IN}/evol1_R2.fastq" \
  -o "${OUTD}/evol1_R1.trimmed.fastq" \
  -O "${OUTD}/evol1_R2.trimmed.fastq" \
  --trim_front2 0 \
  --trim_tail2 0 \
  --detect_adapter_for_pe \
  --qualified_quality_phred 28 \
  --length_required 80 \
  --trim_front1 0 --trim_tail1 0\
  --cut_front --cut_tail --cut_mean_quality 28 \
  --thread 4 \
  --html "${QCD}/fastp_evol1.html" \
  
fastp \
  -i  "${IN}/evol2_R1.fastq" -I "${IN}/evol2_R2.fastq" \
  -o "${OUTD}/evol2_R1.trimmed.fastq" \
  -O "${OUTD}/evol2_R2.trimmed.fastq" \
  --trim_front2 0 \
  --trim_tail2 0 \
  --detect_adapter_for_pe \
  --qualified_quality_phred 28 \
  --length_required 80 \
  --trim_front1 0 --trim_tail1 0\
  --cut_front --cut_tail --cut_mean_quality 28 \
  --thread 4 \
  --html "${QCD}/fastp_evol2.html" \
  
echo "[OK] trimming realizado con éxito en ${OUTD}"

