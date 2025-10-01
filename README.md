# Proyecto Genoma Bioinformática
En este proyecto analizaremos datos de secuenciación de _Escherichia coli_ obtenidos de un ensayo de evolución experimental. A partir de una población ancestral, las bacterias fueron sometidas a distintas condiciones de cultivo durante varias generaciones, lo que favorece la aparición de mutaciones que pueden conferir ventajas adaptativas y la comparación entre la población ancestral y las derivadas permite identificar los cambios genómicos responsables de esas adaptaciones. Para este análisis se simulará un flujo de trabajo, mostrado a continuación:


<img src="https://github.com/user-attachments/assets/7449d4cb-1cd1-41a3-949a-95f2c5e8944f" width="600" />


# **Pregunta Biológica**
¿Qué cambios genómicos surgen en _Escherichia coli_ bajo condiciones de cultivo controladas a lo largo de un experimento de evolución, y cómo estos cambios pueden estar relacionados con ventajas adaptativas en el fenotipo?

# 🗂 Estructura del repositorio
```
├── scripts/          # Scripts para la ejecución del análisis
├── quality_results/  # Reportes de calidad generados (FASTQC/MultiQC HTML)
│ ├── pre_qc/
│ ├── post_qc/
│ └── fastp/
├── assembly_results/ # scaffold.fasta + reporte QUAST
├── index/            # Archivos BAM procesados e indexados
├── informe/          # Informe PDF con análisis biológico
└── README.md         # Este archivo
```
# Prerequisitos para correr Scripts 
1. Crear un environment con todos los paquetes descargados, para eso coloque en una carpeta los scripts en conjunto con el archivo environment.yml y corra el siguiente codigo
   conda env create -f environment.yml
2. Active en entorno recien creado
   conda activate qc-reads
3. Dele premiso  de ejecucion a los scripts
  chmod +x *.sh
4. corra cada srcipt en orden numerico, ejemplo de codigo
 ./ejemplo.sh 

# 🚀 Ejecución del Flujo de trabajo 
Este repositorio contiene un conjunto de scripts bash que implementan un pipeline reproducible para analizar los datos de secuenciación en un contexto de evolución experimental, desde la estructura en los archivos inicial hasta el análisis de mapeo.
| Script           | Descripción                                                                        |
| ---------------- | ---------------------------------------------------------------------------------- |
| `setup.sh`    | Configura la estructura de carpetas, para mantener los reportes organizados        |
| `02_pre_qc.sh`       | Ejecuta **FastQC** y **MultiQC** sobre las lecturas crudas.                        |
| `03_trimming.sh` | Aplica **fastp** para trimming y filtrado por calidad.                             |
| `05_assemble_mapping.sh` | Ensamblaje de novo con **SPAdes** y evaluación con **Quast** y alineamiento de lecturas evolucionadas con **BWA** + procesamiento con `samtools`. |

# Explicación de códigos y flags importantes

## trimming_02

**fastp** comando que realiza el control de calidad para datos de secuenciación.

- `--detect_adapter_for_pe`  
  Detecta adaptadores con algoritmos específicos para *paired-end* (lecturas pareadas). Por defecto, esta detección automática está desactivada para datos *paired-end* y activarla puede mejorar la limpieza de las secuencias.

- `-qualified_quality_phred 28`  
  Elimina reads que tengan mas del 40 % de las bases por debajo de un valor phred de 28 

- `--length_required 80`  
  Retira cualquier lectura que tenga una longitud menor a 80 pares de bases (pb).

- `--cut_front --cut_tail --cut_mean_quality 28`  
  Corta bases al principio y al final de la secuencia ve la calidad promedio de las 4 primeras y ultimas bases si esta esta por debajo de un valor phred de 28 las corta y analisa las siguientes 4 bases si el valor phred esta por encima continua con el siguiente

## Assembly

**SPAdes** es un paquete que permite ensamblar genomas *de novo* o utilizando un genoma de referencia.

- `-1` indica el archivo con las lecturas de la hebra forward (forward strand).  
- `-2` indica el archivo con las lecturas de la hebra reverse (reverse strand).  
- `--isolate` es un flag recomendado para datos de Illumina; ajusta parámetros internos optimizados para lecturas Illumina y para ensamblajes de aislados bacterianos.

seqkit seq -m 100 scaffolds_anc_trimmed.fasta > filter_scaffold_anc_trimmed.fasta

Este comando elimina los scaffolds de menos de 100 pb, ya que:

- Scaffolds con menos de 100 bases no permiten identificar genes correctamente.
- Además, suelen ser más cortos que las lecturas que se mapearán sobre ellos.
- Los scaffolds filtrados (de longitud mayor o igual a 100) se guardan en un nuevo archivo FASTA (`filter_scaffold_anc_trimmed.fasta`).
  
## Mapping 

bwa index -p filter_scaffold_anc_raw ../filter_scaffold_anc_raw.fasta
Indexa o ordena los scaffolds para poder luego generar el archivo bam  
-p les pone al todos los archivos creados el nombre filter_scaffold_anc_raw 

 `bwa mem`  
  Ejecuta el algoritmo BWA-MEM que crea el archivo SAM 

- `assemble/extra_info/filter_scaffold_anc_trimmed`  
  Prefijo de los archivos índice del genoma de referencia (los archivos generados previamente con `bwa index`).

- `data/trimmed_data/evol1_R1.trimmed.fastq`  
  Archivo FASTQ con las lecturas de la hebra forward 

- `data/trimmed_data/evol1_R2.trimmed.fastq`  
  Archivo FASTQ con las lecturas de la hebra reverse 

**Samtools** 

samtools view -b -F 4 -q 30 mapping/aligned_raw_evol1.sam | \
samtools sort -n -o mapping/trash.bam -

samtools fixmate -m mapping/trash.bam \
mapping/trash.fixmate.bam

samtools sort -o - mapping/trash.fixmate.bam | \
samtools markdup -r - mapping/filtered_aligned_raw_evol1.bam
rm -rf mapping/trash* 


- `-b`  
  Genera la salida en formato BAM (binario), más eficiente para procesamiento.

- `-F 4`  
  Excluye las lecturas con el flag 4, que indica lecturas **no mapeadas**.

- `-q 30`  
  Excluye lecturas con calidad de mapeo menor a 30.

- Combinación de `sort -n` y `fixmate`:  
  Ordena las lecturas por nombre para asegurar que pares estén juntos y corrige información del par, preparando los datos para la eliminación de duplicados.

- `markdup -r`  
  Marca y remueve lecturas duplicadas.

- `rm -rf mapping/trash*`  
  Elimina archivos temporales generados en el proceso.
