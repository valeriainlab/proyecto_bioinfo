#!/usr/bin/env bash
set -euo pipefail
IN="data/trimmed_data"
OUT="qc/post_qc"
mkdir -p "${OUT}"

fastqc -o "${OUT}" -t 4 ${IN}/*.trimmed.fastq


multiqc -o "${OUT}" "${OUT}"
echo "[OK] QC post trimming en ${OUT}"