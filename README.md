# Code for nickel adaptations manuscirpt

## Pangenomics

### 1 - Preparing prokka-annotated genomes using gff and fasta files
Install anvio v7
download gff_parser.py as instructed [here](https://merenlab.org/2017/05/18/working-with-prokka/)
```
for i in `ls *fna | awk 'BEGIN{FS=".fna"}{print $1}'`
do
    python ./../gff_parser.py $i.gff --gene-calls $i.gene_calls.txt --annotation $i.gene_annot.txt
    anvi-gen-contigs-database -f $i.fna  -o $i.contigs.db --external-gene-calls $i.gene_calls.txt --ignore-internal-stop-codons
    anvi-import-functions -c $i.contigs.db \
                         -i $i.gene_annot.txt
done
```
### 2 - Generate genomes database
Prepare list of external genomes as described [here](https://anvio.org/help/7/artifacts/external-genomes/) as list.txt

```
anvi-gen-genomes-storage -e list.txt --output-file all-GENOMES.db --gene-caller Prodigal
```

### 3 - Run anvi'o pangenome
```
anvi-pan-genome -g all-GENOMES.db -n "all" \
--num-threads 16 --output-dir all-GENOMES --mcl-inflation 6 \
--minbit 0.6  --skip-homogeneity  --skip-hierarchical-clustering  \
--skip-alignments

```
### 4 - obtain pangenome presence/absence table
```
anvi-export-table all-GENOMES.db --table gene_cluster_presence_absence

```
## GWAS
Download pyseer version 1.3.6
We used a linear mixed model (FaST-LMM)
Clade 1 and Clade 2 run seperately
```
python pyseer/scripts/phylogeny_distance.py --calc-C RAxML_bestTree.protein_fas2.tre > phylogeny_PDDIST.tsv
pyseer --lmm  --phenotypes cl1.mean.od.mic.txt --phenotype-column MIC --similarity phylogeny_PDDIST.tsv --pres cl1_pantable.txt   --cpu 8 > MIC_cl1_lmm.txt
```


