#!/usr/bin/env bash
set -euo pipefail
# Crear la estructura de carpetas
mkdir -p entregables/{quality_results/{qc_pre,qc_post,fastp},assembly_results,index,informe}

# Copiar los archivos
#qc
cp qc/pre_qc/*.html entregables/quality_results/qc_pre/
cp qc/fastp/*.html entregables/quality_results/fastp/
cp qc/post_qc/*.html entregables/quality_results/qc_post/

#Scaffold y Quast 
cp assemble/*.fasta entregables/assembly_results/
cp assemble/quast_results/results_2025_09_30_22_52_13/report.html entregables/assembly_results/

#BAM
cp mapping/filtered_aligned_trimmed_evol1.bam entregables/index/
cp mapping/filtered_aligned_trimmed_evol2.bam entregables/index/
