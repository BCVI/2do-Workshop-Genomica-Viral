# 2do Workshop: Bioinformática en Genómica Viral
## Sección 3: Metagenómica viral 
Febrero 20, 2025
Luis Chica

### Parte 0: Datos crudos
Usaremos una de las muestras generadas en:
Rose, G., Shaw, A. G., Sim, K., Wooldridge, D. J., Li, M. S., Gharbia, S., ... & Kroll, J. S. (2017). Antibiotic resistance potential of the healthy preterm infant gut microbiome. PeerJ, 5, e2928.

Utilizando las rutas descritas a continuación, encontrarán dos subdirectorios llamados `reads` y `contigs`. 
- **`/home/luis/Dia3_metagenomicaViral/reads`**: Lecturas crudas en formato FASTQ.
- **`/home/luis/Dia3_metagenomicaViral/contigs`**: Contigs ensamblados y archivos de anotación.

## Parte 1: Preprocesamiento de secuencias y ensamblaje 
#### Paso 1: Análisis de calidad y limpieza
Antes del ensablaje de los genomas virales, es necesario evaluar la calidad y eliminar las lecturas de baja calidad. Para esto, usaremos las herramientas de **FastQC** y **TrimGalore**.

Debemos crear un soft link de las secuencias para no tener problemas de sobreescritura o borrado de datos

```bash
ln -s /home/luis/Dia3_metagenomicaViral/reads/ERR1600426_subsample01_R1.fastq /home/usuarioX/ERR1600426_subsample01_R1.fastq
ln -s /home/luis/Dia3_metagenomicaViral/reads/ERR1600426_subsample01_R2.fastq /home/usuarioX/ERR1600426_subsample01_R2.fastq
``` 
**Importante: cambiar la ruta de destino por la ruta de cada usuario**

#### Verificación del número de lecturas
Podemos verificar cuántas lecturas hay en cada archivo (R1 y R2) usando el comando `grep` con el flag `-c`. Para esto, buscamos un patrón único para cada lectura, como `@ERR`.

```bash
grep -c "@ERR" *.fastq
```

**Descripción del comando**:
- `grep -c`: Cuenta el número de líneas que coinciden con el patrón especificado.
- `@ERR`: Patrón que identifica el inicio de una lectura en archivos FASTQ. Puede cambiar
- `*.fastq`: Busca en todos los archivos con extensión `.fastq`.

**Salida esperada**:
```
Reads-R1.fastq: 698503
Reads-R2.fastq: 698503
```
#### Calidad de las secuencias
Para revisar la calidad de las secuencias usaremosr la herramienta **FastQC**.

```bash
fastqc ERR1600426_subsample01_R1.fastq ERR1600426_subsample01_R2.fastq 
```

**Descripción del comando**:
- `fastqc`: Herramienta que genera un informe de calidad y contenido de adaptadores de las lecturas en formato HTML.
- `ERR1600426_subsample01_R1.fastq` y `ERR1600426_subsample01_R2.fastq1`: Archivos de lecturas crudas.

**Salida esperada**:
- Se generarán dos archivos HTML (uno para cada archivo FASTQ) que contienen gráficos y estadísticas sobre la calidad de las lecturas.

#### Filtrado de lecturas por calidad
Para eliminar lecturas de baja calidad, usamos **TrimGalore**.

```bash
trim_galore --fastqc -q 30 --paired ERR1600426_subsample01_R1.fastq ERR1600426_subsample01_R2.fastq -o Quality_Reads
```

**Descripción del comando**:
- `trim_galore`: Herramienta que combina `Cutadapt` y `FastQC` para filtrar lecturas.
- `--fastqc`: Ejecuta FastQC después del filtrado.
- `-q 30`: Elimina bases con calidad menor a 30.
- `--paired`: Indica que los archivos de entrada provienen de paired-end-sequencing.
- `-o Quality_Reads`: Directorio de salida para las lecturas filtradas.

**Salida esperada**:
- Archivos filtrados: `Reads-R1_val_1.fq` y `Reads-R2_val_2.fq`.
- Informes de calidad: Archivos HTML generados por FastQC.

### Paso 2: Ensamblaje metagenómico
Para ensamblar los genomas virales a partir de las lecturas filtradas, usamos **SPAdes** en modo metagenómico. 

```bash
spades.py -1 Quality_Reads/ERR1600426_subsample01_R1_val_1.fq -2 Quality_Reads/ERR1600426_subsample01_R2_val_2.fq --meta -o spades -t 12
```
**Descripción del comando**:
- `spades.py`: Herramienta de ensamblado de genomas.
- `-1` y `-2`: Especifican los archivos de lecturas pareadas.
- `--meta`: Modo metagenómico para ensamblar genomas de múltiples organismos.
- `-o spades`: Directorio de salida para los resultados.
- `-t 12`: Usa 12 threads para acelerar el proceso.

**Salida esperada**:
- Archivo de contigs: `contigs.fasta`.
- Estadísticas de ensamblado: Archivos de texto y gráficos en el directorio de salida.

### Paso 3: Evaluación del ensamblaje
Para evaluar la calidad del ensamblado, usamos **QUAST**.

```bash
quast.py spades/contigs.fasta -o quast_results
```
**Descripción del comando**:
- `quast.py`: Herramienta que evalúa la calidad de los ensamblados.
- `contigs.fasta`: Archivo de contigs generado por SPAdes.
- `-o quast_results`: Directorio de salida para los resultados.

**Salida esperada**:
- Informe con métricas de calidad (número de contigs, longitud, N50, etc.).

Gracias a las estadísticas de calidad, se observa una gran cantidad de contigs cortos (< 1kb), normalmente estos contigs son filtrados ya que su información genómica es mínima. 

## Parte 2: Identificación de elementos virales

Para esta parte del tutorial usaremos contigs previamente ensamblados de un dataset completo. Estos contigs tuvieron el mismo procesamiento descrito anteriomente, mas un filtro de longitud de 2kb. Ya que comunmente se trabaja con más de una muestra, los ensamblajes por muestra se concatenan y se clusterizan los contigs para remover redundancia.

  ```bash
ln -s /home/luis/Dia3_metagenomicaViral/contigs/contigs_fullSet.fasta /home/usuarioX/contigs_fullSet.fasta
   ```
### Paso 1: Jaeger
**Jaeger** es una herramienta que utiliza modelos de machine learning (redes neuronales convolucionales) para identificar secuencias virales en contigs. Es especialmente útil para detectar virus en metagenomas complejos, donde las herramientas tradicionales basadas en homología pueden fallar.

1. **Activar el entorno de Jaeger**:
   ```bash
   conda activate jaeger
   ```
2. **Ejecutar Jaeger**:
   Usa el siguiente comando para analizar los contigs:
   ```bash
   jaeger run -i contigs_fullSet.fasta -o jaeger -s 2.5 --fsize 1000 
   ```
   **Descripción de los parámetros**:
   - `-i contigs_fullSet.fasta`: Archivo de entrada con los contigs.
   - `-o jaeger`: Directorio de salida donde se guardarán los resultados.
   - `-s 2.5`: Sensitividad del algoritmo de de extracción (valores entre 0 y 4).
   - `--fsize 1000`: Tamaño de la ventana de análisis (1000 pb).

2. **Salida esperada**:

En el directorio de salida (`jaeger`), encontraremos el archivo: **<input_file>_default.jaeger.tsv**. ste archivo contiene las predicciones para cada contig, con varias columnas que proporcionan información detallada sobre las secuencias. De toda la información contenida en el archivo, las 8 primeras columnas proporcionan la mayor información para análisis posteriores.

- contig_id: Nombre o identificador del contig.
- length: Longitud del contig en pares de bases (pb).
- prediction: Predicción principal de Jaeger: Phage, bacteria....
- entropy: Medida de la entropía de la secuencia, que indica su complejidad. Valores más altos sugieren mayor diversidad de nucleótidos.
- reliability_score: valor de confianza de la predicción.
- host_contam: Indica si hay contaminación por secuencias del hospedero.
- prophage_contam: Indica si hayseñales de profagos (virus integrados en genomas bacterianos).
- G+C: Contenido de guanina y citosina (G+C) en el contig.

### Paso 2: PhiSpy
PhiSpy es una herramienta diseñada específicamente para la detección de profagos en genomas bacterianos anotados. Utiliza una combinación de heurísticas y machine learning para identificar regiones que presentan características típicas de profagos, como genes asociados a la integración y excisión, desviaciones en el contenido de G+C y la presencia de sitios de recombinación específicos.

Uno de los inputs requeridos por PhiSpy es la anotación de regiones codificantes en formato GenBank. Para esto es necesario anotar nuestro set de contigs. 

Para esta parte usaremos un archivo .gbk previamente generado. 
  ```bash
ln -s /home/luis/Dia3_metagenomicaViral/contigs/Streptococcus_pyogenes_M1_GAS.gb /home/usuarioX/Streptococcus_pyogenes_M1_GAS.gb
   ```
1. **Activar el entorno de phispy**:
   ```bash
   conda activate phispy
   ```
2. **Correr comando principal**:
```bash
conda activate phispy
phispy -o phispy --threads 12 -u 1000 Streptococcus_pyogenes_M1_GAS.gb
```
**Descripción del comando**:
- `-o phispy`: Directorio de salida.
- `--threads 12`: Usa 12 hilos.
- `Streptococcus_pyogenes_M1_GAS.gbk`: Archivo de entrada en formato GenBank. Argumento posicional
- `-u 1000`: Tamaño mínimo de región para identificar como profago.

**Salida esperada**:
- prophage_coordinates.tsv: Archivos con las coordenadas de las regiones identificadas como profagos.

Las columnas de dicho archivo son: 
- Número del profago (ID)
- Contig del cual se identifica el profago
- Coordenada de inicio del profago. 
- Coordenada de fin del profago. 
Si se detcectan sitios de recombinación, tambien se describe: 
- Incio de attL
- Fin de attL
- Incio de attR
- Fin de att
- Secuencia de attL;
- Secuencia de attR;
- Explicación de por que se selecionó el sitio de integración.

### Paso 3: geNomad

Recientemente se han desarrollado herramientas que predicen secuencias virales en metagenomas y, a su vez, predicen profagos y elementos extracromosomales como plásmidos. Una de las más recientes es *geNomad*, que usa una combinación de información genómica y de regiones codificantes, así como bases de datos de "marker proteins"

1. **Activar el entorno de phispy**:
   ```bash
   conda activate genomad
   ```

2. **Correr comando principal**:
```bash
genomad end-to-end --min-score 0.6 --cleanup --threads 12 contigs_fullSet.fasta geNomad /6TB/DBs/genomad_db
```
**Descripción del comando**:
- `--min-score 0.6`: Filtra elementos con un score mínimo de 0.6.
- `--cleanup`: Elimina archivos temporales.
- `--threads 12`: Usa 12 threads.
- `geNomad`: Directorio de salida.
- `/6TB/DBs/genomad_db`: Ruta a la base de datos.

La salida esperada se discutirá al final del la sesión.  

**Con esto concluimos el análisis de metagenómica viral en este taller. Hemos recorrido desde la preprocesamiento de secuencias hasta la identificación de elementos virales, aplicando herramientas clave como FastQC, TrimGalore, SPAdes, QUAST, Jaeger y Phispy. Esperamos que esta guía  sirva como referencia para futuros análisis. ¡Gracias por participar!**
