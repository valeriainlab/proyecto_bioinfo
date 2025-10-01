#!/usr/bin/env bash
set -euo pipefail
IN="data/raw_data"
OUT="qc/pre_qc"
mkdir -p "${OUT}"

fastqc -o "${OUT}" -t 4 \
  "${IN}/anc_R1.fastq" \
  "${IN}/anc_R2.fastq" \
  "${IN}/evol1_R1.fastq" \
  "${IN}/evol1_R2.fastq" \
  "${IN}/evol2_R1.fastq" \
  "${IN}/evol2_R2.fastq"

multiqc -o "${OUT}" "${OUT}"
echo "[OK] QC inicial en ${OUT}"