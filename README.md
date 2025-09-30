# Proyecto Genoma Bioinformática
En este proyecto analizaremos datos de secuenciación de _Escherichia coli_ obtenidos de un ensayo de evolución experimental. A partir de una población ancestral, las bacterias fueron sometidas a distintas condiciones de cultivo durante varias generaciones, lo que favorece la aparición de mutaciones que pueden conferir ventajas adaptativas y la comparación entre la población ancestral y las derivadas permite identificar los cambios genómicos responsables de esas adaptaciones. Para este análisis se simulará un flujo de trabajo, mostrado a continuación:


<img src="https://github.com/user-attachments/assets/7449d4cb-1cd1-41a3-949a-95f2c5e8944f" width="600" />


# **Pregunta Biológica**
¿Qué cambios genómicos surgen en _Escherichia coli_ bajo condiciones de cultivo controladas a lo largo de un experimento de evolución, y cómo estos cambios pueden estar relacionados con ventajas adaptativas en el fenotipo?

# 🗂 Estructura del repositorio
```
├── scripts/            # Scripts para la ejecución del análisis
├── quality_results/    # Reportes de calidad generados (FASTQC/MultiQC HTML)
├── assembly_results/   # scaffold.fasta + reporte QUAST
├── index/              # Archivos BAM procesados e indexados
├── informe/            # Informe PDF con análisis biológico
└── README.md           # Este archivo
```

# 🚀 Ejecución del Flujo de trabajo 
Este repositorio contiene un conjunto de scripts bash que implementan un pipeline reproducible para analizar los datos de secuenciación en un contexto de evolución experimental, desde la estructura en los archivos inicial hasta el análisis de mapeo.
| Script           | Descripción                                                                        |
| ---------------- | ---------------------------------------------------------------------------------- |
| `setup_00.sh`    | Configura la estructura de carpetas, para mantener los reportes organizados y ejecuta el entorno `qc-reads` que contiene todos los programas necesarios para este análisis      |
| `qc_01.sh`       | Ejecuta **FastQC** y **MultiQC** sobre las lecturas crudas.                        |
| `trimming_02.sh` | Aplica **fastp** para trimming y filtrado por calidad.                             |
| `assembly_03.sh` | Ensamblaje de novo con **SPAdes** y evaluación con **Quast**.                      |
| `mapping_04.sh`  | Alineamiento de lecturas evolucionadas con **BWA** + procesamiento con `samtools`. |
