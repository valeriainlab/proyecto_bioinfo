#!/usr/bin/env bash
set -euo pipefail

mkdir -p proyecto_grupo4
mv scripts/ proyecto_grupo4
cd proyecto_grupo4

wget -O data.tar.gz https://osf.io/2jc4a/download
tar -xvzf data.tar.gz
cd data
gunzip *.gz 
mkdir -p raw_data trimmed_data
mv *.fastq raw_data/
cd ..

mkdir -p qc/{pre_qc,post_qc,fastp}  assemble mapping scripts 

echo "Todo listo!! Ahora si, empecemos."