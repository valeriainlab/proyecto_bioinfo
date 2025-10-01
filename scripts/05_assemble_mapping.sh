#!/usr/bin/env bash
set -euo pipefail


# Crear carpeta para resultados de ensamblaje y archivos extra
mkdir -p  assemble/extra_info/

# Corre SPAdes con datos recortados (trimmed), indicando archivos paired-end -1 y -2
# --isolate recomendado para aislamientos bacterianos, optimiza parámetros para Illumina
spades.py -o assemble/extra_info -1 data/trimmed_data/anc_R1.trimmed.fastq \
-2 data/trimmed_data/anc_R2.trimmed.fastq --isolate

# Mover archivo scaffolds a directorio principal y renombrar para evitar confusión
cd assemble/extra_info 
mv scaffolds.fasta ../
cd ..
mv scaffolds.fasta scaffolds_anc_trimmed.fasta 

# Quita scaffolds de menos de 100 pb de tamaño para mejorar calidad y evitar secuencias poco útiles
seqkit seq -m 100 scaffolds_anc_trimmed.fasta > filter_scaffold_anc_trimmed.fasta

cd ..

# Corre SPAdes con datos sin recortar (raw), mismos parámetros
spades.py -o assemble/extra_info -1 data/raw_data/anc_R1.fastq \
-2 data/raw_data/anc_R2.fastq --isolate

# Mover scaffolds y renombrar para separar ensamblajes raw y trimmed
cd assemble/extra_info 
mv scaffolds.fasta ../
cd ..
mv scaffolds.fasta scaffolds_anc.fasta 

# Filtra scaffolds cortos para ensamblaje raw también
seqkit seq -m 100 scaffolds_anc.fasta > filter_scaffold_anc_raw.fasta

echo "[OK] genoma ensamblado en scaffolds, raw_data y trimmed_data en assemble/"

# Comparar calidad entre ensamblajes filtered raw y trimmed usando QUAST
quast.py filter_scaffold_anc_raw.fasta filter_scaffold_anc_trimmed.fasta

echo "[OK] Quast instalado"

cd extra_info

bwa index -p filter_scaffold_anc_raw ../filter_scaffold_anc_raw.fasta

bwa index -p filter_scaffold_anc_trimmed ../filter_scaffold_anc_trimmed.fasta

cd ../../
# ensamble del genoma en archivo sam, paired end entregando scafolds y reads 
bwa mem -t 8 assemble/extra_info/filter_scaffold_anc_trimmed \
data/trimmed_data/evol1_R1.trimmed.fastq \
data/trimmed_data/evol1_R2.trimmed.fastq > mapping/aligned_trimmed_evol1.sam


bwa mem -t 8 assemble/extra_info/filter_scaffold_anc_raw \
data/trimmed_data/evol1_R1.trimmed.fastq \
data/trimmed_data/evol1_R2.trimmed.fastq > mapping/aligned_raw_evol1.sam


bwa mem -t 8 assemble/extra_info/filter_scaffold_anc_trimmed \
data/trimmed_data/evol2_R1.trimmed.fastq \
data/trimmed_data/evol2_R2.trimmed.fastq > mapping/aligned_trimmed_evol2.sam

bwa mem -t 8 assemble/extra_info/filter_scaffold_anc_raw \
data/trimmed_data/evol2_R1.trimmed.fastq \
data/trimmed_data/evol2_R2.trimmed.fastq > mapping/aligned_raw_evol2.sam

#depuracion por calidad -q 30 quitando reads no mapeados -F 4 
# con fixmate y sort y markdup es para quitar repetidos 
samtools view -b -F 4 -q 30 mapping/aligned_raw_evol1.sam | \
samtools sort -n -o mapping/trash.bam -

samtools fixmate -m mapping/trash.bam \
mapping/trash.fixmate.bam

samtools sort -o - mapping/trash.fixmate.bam | \
samtools markdup -r - mapping/filtered_aligned_raw_evol1.bam
rm -rf mapping/trash* 

samtools view -b -F 4 -q 30 mapping/aligned_raw_evol2.sam | \
samtools sort -n -o mapping/trash.bam -

samtools fixmate -m mapping/trash.bam \
mapping/trash.fixmate.bam

samtools sort -o - mapping/trash.fixmate.bam | \
samtools markdup -r - mapping/filtered_aligned_raw_evol2.bam
rm -rf mapping/trash* 

 samtools view -b -F 4 -q 30 \
  mapping/aligned_trimmed_evol1.sam\
| samtools sort -o mapping/filtered_aligned_trimmed_evol1.bam

 samtools view -b -F 4 -q 30 \
  mapping/aligned_trimmed_evol2.sam\
| samtools sort -o mapping/filtered_aligned_trimmed_evol2.bam



cd mapping 

# archivo que requiere el qualimap para colocarme los dos en un mismo reporte 

qualimap bamqc -bam filtered_aligned_trimmed_evol1.bam -outdir Qualimap_evo1

qualimap bamqc -bam filtered_aligned_trimmed_evol2.bam -outdir Qualimap_evo2


mkdir -p unmapped_reads 

samtools view -b -f 4 aligned_trimmed_evol1.sam > unmapped_reads/unmapped_evol1.bam
samtools view -b -f 4 aligned_trimmed_evol2.sam > unmapped_reads/unmapped_evol2.bam

cd ..

echo "[OK] assemble and mapping"