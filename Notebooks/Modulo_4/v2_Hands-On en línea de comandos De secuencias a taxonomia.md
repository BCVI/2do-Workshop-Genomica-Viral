# Hands-On en línea de comandos: De secuencias a taxonomia

A partir de secuencias de fagos, se va a realizar la clasificación, clustering y anotación taxonómica.  Vamos a usar dos archivos con grupos de fagos diferentes.  

- `phages_arthrobacter_v2.fasta`: fagos que infectan *Arthrobacter*
- `phages_gordonia.fasta`: fagos que infectan *Gordonia*

## 1. Realizar clasificación y clustering usando el genoma (VIRIDIC) 

la herramienta **VIRIDIC** no solo se puede cirrer usando la interfaz web, en caso de tener una gran cantidad de contigs/genomas es mejor correrlo usando un servidor de HPC (High-Performance Computing).  

Ahora en su `home`cree un diretorio para la el modulo 4 clusterización y taxonomia

```bash
#Entrar al directorio de home, $HOME es una variable predefinida que contiene el path hacia su home. 
cd $HOME 
#crear directorio para modulo 4
mkdir dia4_clustering_taxonomia
#copiar los contigs en el nuevo directorio
cp /home/came/dia4_clustering_taxonomia/phages_*.fasta $HOME/dia4_clustering_taxonomia
```

Para correr VIRIDIC en el servidor, se necesita como entrada las secuencias fasta en un solo archivo (igual que en el servidor web).

```bash
#Atención: este comando no deja autocompletar así que deben estar muy seguros de como estan escritos los directorios. 
#projdir = nombre del directorio donde quieren guardar la salida de VIRIDIC
#in = archivo con las secuencias fasta.
/opt/viridic_v1.1/viridic.bash projdir=/home/came/dia4_clustering_taxonomia/resultados_viridic_arthrobacter.out in=/home/came/dia4_clustering_taxonomia/phages_arthrobacter_v2.fasta
```

 Ahora tiene un directorio llamado `resultados_viridic_arthrobacter.out/`donde se encuentran todos los archivos temporales y de salida de VIRIDIC. Dentro de este directorio el resultado final se encuentra en `04_VIRIDIC_out`. 

```bash
#Para ver la matriz de similitud 
cat $HOME/dia4_clustering_taxonomia/resultados_viridic_arthrobacter.out/04_VIRIDIC_out/sim_MA_genCol.csv
#Para ver los clusters generados
cat $HOME/dia4_clustering_taxonomia/resultados_viridic_arthrobacter.out/04_VIRIDIC_out/clusters.csv
```

¿Cuántos clusters se generaron? ¿Coincide con la matriz de similitud?

Ahora para visualizar las figuras con la agrupación de genomas y los valores de identidad generados por VIRIDIC, se descargan los siguientes archivos: 

```bash
#RECUERDE REMPLAZAR MI USUARIO (came) POR EL SUYO 
scp came@132.248.32.20:$HOME/dia4_clustering_taxonomia/resultados_viridic_arthrobacter.out/04_VIRIDIC_out/Heatmap.PDF heatmap_viridic_arthrobacter.PDF
```

Las figuras concuerdan con los clusters y matriz de similitud generados ? 

**Actividad (opcional). Realice el mismo análisis con los fagos de Gordonia.** 

## 2. Generación de un dendrograma/árbol filogénetico a partir de la matriz de similitud. 

**Actividad de preparación. ¿Cómo generar un árbol/dendrograma?  ** Descargar documento **Ejericicio_ClusteringJerarquico.docx** del github. En este ejercicio se usa el método de agrupación "complete", el usado en árboles filogéneticos es [UPGMA](http://teacheng.illinois.edu/PhylogeneticTree/#/upgma) o "average"

Otro algoritmo de agrupación en árboles es Neighbour joining: [Puede revisar más información aquí](https://www.tenderisthebyte.com/blog/2022/08/31/neighbor-joining-trees/)

Para completar este paso vamos a usar el software R. Para abir R: 

```bash
cd $HOME/dia4_clustering_taxonomia
mkdir trees_nucleotideLevel
cd trees_nucleotideLevel
R 
```

En R debe leer el archivo de matriz de similitud, pasarlo a matriz de distancias y generar un dendrograma. Tenga en cuenta que en R la variable $HOME no funciona, asi que debe copiar la dirección al archivo con su nombre de usuario. 

```R
#Importar la matriz de similitud y volverla una matriz numerica
#RECUERDE REMPLAZAR MI USUARIO (came) POR EL SUYO 
mat_similitud=read.table("/home/came/dia4_clustering_taxonomia/resultados_arthrobacter.out/04_VIRIDIC_out/sim_MA_genCol.csv",sep="\t",h=T)
rownames(mat_similitud)=mat_similitud$genome
mat_similitud=mat_similitud[,-1]
#Pasarla a distancias de 0 a 1
mat_distancias=1-(mat_similitud/100)
#Tansformarla a un objeto de distancias en R
distancias <- as.dist(mat_distancias)

#Crear un arbol con UPGMA 
tree_upgma <- hclust(distancias, method = "average")

#Generar dendrograma
png("upgma_tree_arthrobacter.png")
plot(tree_upgma)
dev.off()

#Cerrar sesion
q() #presionar Y seguido de Enter. 
```

Ahora, puede descargar los archivos para visualizarlos en su computador

```bash
#RECUERDE REMPLAZAR MI USUARIO (came) POR EL SUYO 
scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/trees_nucleotideLevel/upgma_tree_arthrobacter.png .
```

Ahora puede ver el dendrograma generado, coincide con la matriz de similitud ? muestra resultados similares  a los obtenidos con VIRIDIC ? 

**Actividad (opcional). Realice el mismo análisis con los fagos de Gordonia.** 

## 3. Virclust

Esta herramienta permite hacer clustering a nivel de proteínas en diferentes pasos.

1. **Agrupamiento de proteínas**

​	Step1A. predecir y traducir genes para obtener proteínas

​	Step2A. agrupar protepinas en clusters de proteínas (PCs) usando BLASTp

2. **Generar agrupamiento de genomas basado en distancias intergenóminas basadas en el número de proteínas compartidas**

​	Step3A. Clustering jerarquico de genomas – usa método "complete" pero se puede cambiar a "average". En este paso también se puede hacer Bootstrap resampling para calcular la significancia de cada rama del árbol. 

​	Step4A. Dividir los genomas en clusters de genomas virales (VGCs) y generar estadísticas. 

3. Identificación de proteínas core.

   Step5A. Identifica las proteínas core. 

4. Anotación funcional de proteínas usando distintas bases de datos (este paso no lo cubriremos en el taller por restricciones de tiempo)

```bash
conda activate VirClust
#Correr VirClust
Rscript /6TB/apps/VirClust/VirClust/vir_clust_standalone/VirClust_MASTER.R sing=conda condaenvpath=/opt/miniforge3/envs/VirClust projdir=$HOME/dia4_clustering_taxonomia/resultados_virclust_arthrobacter.out infile=$HOME/dia4_clustering_taxonomia/phages_arthrobacter_v2.fasta step1A=T eval_PC=0.0000000001 cov_PC=69 step2A=T step3A=T boot_pv_a=yes step3A_Plot=T step4A=T step4A_Plot=T clust_dist_a=0.3 step5A=T

#Descargar los resultados 
#RECUERDE REMPLAZAR MI USUARIO (came) POR EL SUYO 

scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/resultados_virclust_arthrobacter.out/04a-06a_genome_clustering_PC/06-Heatmap_PC.PDF . 

scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/resultados_virclust_arthrobacter.out/04a-06a_genome_clustering_PC/04/Dist_heatmap_PC_all_genomes.PDF .

scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/resultados_virclust_arthrobacter.out/04a-06a_genome_clustering_PC/04/04/pv_trees/pv_tree_au.newick
```

 Para visualizar dirijase al sitio web iTOL (interactive Tree Of Life https://itol.embl.de/ ), cree un usuario o inicie sesión. Ingrese a "Upload a Tree", seleccion el archivo pv_tree_au.newick. Ahora puede revisar el árbol. 

**Actividad (opcional). Realice el mismo análisis con los fagos de Gordonia.** 

## 4. Visualización sintenia

Ahora para analizar el proteoma a partir del genoma, se debe realizar primero una anotación, es decir, se deben predecir los genes en la secuencia. Para esto vamos a usar el programa pharokka (diseñado para fagos - https://github.com/gbouras13/pharokka). 

```bash
#Cree un directorio para este nuevo ejercicio
mkdir $HOME/dia4_clustering_taxonomia/proteoma

#Para este ejercicio necesitamos cada genoma del archivo phages_arthrobacter_v2.fasta en un archivo diferente, cópielos en su carpeta

cp -r /home/came/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split $HOME/dia4_clustering_taxonomia/proteoma

#luego correr pharokka, -f solo sirve para que corra rapido, si desean hacer un análisis mas exahustivo deben quitarselo
#ATENCION. Este programa demora en correr al menos 5 minutos por genoma.
conda activate pharokka
cd $HOME/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/

#Vamos a generar el proteoma de solo 3 de los 6 fagos. 
pharokka.py -i $HOME/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/phages_arthrobacter_v2.part_002.fa -o Noely_annot -d /6TB/gama/Workshop_2025/Final_Viral_Contigs/pharokka/pharokka_db/ -t 3 -f

pharokka.py -i $HOME/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/phages_arthrobacter_v2.part_004.fa -o CabbageMan_annot -d /6TB/gama/Workshop_2025/Final_Viral_Contigs/pharokka/pharokka_db/ -t 3 -f

pharokka.py -i $HOME/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split//phages_arthrobacter_v2.part_005.fa -o Corgi_annot -d /6TB/gama/Workshop_2025/Final_Viral_Contigs/pharokka/pharokka_db/ -t 3 -f

#Descargar la salida. 
#RECUERDE REMPLAZAR MI USUARIO (came) POR EL SUYO
scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/Noely_annot/pharokka.gbk Noely_pharokka.gbk

scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/CabbageMan_annot/pharokka.gbk CabbageMan_pharokka.gbk

scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/Corgi_annot/pharokka.gbk Corgi_pharokka.gbk

```

Puede ver la sintenia de las proteinas usando  la herramienta clinker https://cagecat.bioinformatics.nl/tools/clinker. Para Clinker se necesita un archivo .gbk por cada genoma, pharokka puede anotar usando el archivo `phages_arthrobacter_v2.fasta` directamente. 

En clinker diríjase a "Genome files", seleccioné todos los archivos .gbk generados >> Submit. Luego de unos minutos podrá ver como las proteínas de los 3 genomas están organizadas de manera similar (sintenia), a pesar de que no sean idénticas. 

Pharokka también permite visualizar la anotación de cada genoma 

```
pharokka_plotter.py -i phages_arthrobacter_v2.part_002.fa -n pharokka_plot -o Noely_annot/

scp came@132.248.32.20:/home/came/dia4_clustering_taxonomia/proteoma/phages_arthrobacter_v2.fa.split/Noely_annot/pharokka_plot.png Noely_plot.png
```

**Actividad (opcional). Realice el mismo análisis con el resto de los fagos de Arthrobacter y luego con los fagos de Gordonia.** 

```bash
#Para los fagos de gordonia:

cp -r /home/came/dia4_clustering_taxonomia/proteoma/phages_gordonia.fasta.split $HOME/dia4_clustering_taxonomia/proteoma
```

