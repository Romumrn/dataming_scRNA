# Projet : Genes network with scRNA data

## Why ?

A gene regulatory network represents the interactions between biological chemicals (messenger RNAs (mRNAs) and proteins. These mRNA and proteins interact with each other with various degrees of specificity. When portrayed as a network, nodes will represent genes, and each arrow represents the links between those genes. The aim of this script is to easily compute a representation of a gene network .


## Fonctionnement

This script use Apriori algorithm to create assocation rules between each genes. It proceeds by identifying the frequent individual gene in the gene count matrix (previously transformed into boolean all in order to remove non expressed genes and those with count too small ) and extending them to larger and larger group of genes. 

Apriori algorihm uses 3 values : 
- Support : Fraction of cells who expressed genes (or a group of genes)
    -  Support(gene1) = Number of cells who express gene1 / total cells
- Confidence :  Fraction of gene expressed if another gene (or a group of genes) is expressed
    - Confidence( Gene1 -> gene2 ) = Support( gene1 , gene2 ) / Support (gene1)
- Lift : Measure of correlation, between a gene and another gene (or groups of genes)
    - lift(Gene1 -> gene2 ) = Support( gene1 , gene2 ) / Support (gene1 )* Support( gene2)


## Script

This script allows you to analyze single cell RNA data in order to find relations between genes and create a network graph of gene expression. Based on the Apriori algorithm, commonly used for both frequent itemset mining and association rule-learning over relational databases. It proceeds by identifying the frequent individual items in the database and extending them to larger and larger item sets as long as those item sets appear sufficiently often in the database (in our case an item is a gene and an itemsets is a groupe of genes).

Exemple of count matrix : 

|   | gene 1| gene 2   | gene 3   | gene 4  |
|---|---|---|---|---|
| Cell 1  |  0 |  1000 | 10000  | 0   |
| Cell 2  |  0 |  0 | 9000  |   0|
| Cell 3  |  0 |  1300 |  11000 |   0|


Support threshold allows the user to control gene frequency, maximum support allows to remove genes who are expressed in all cells (default=0.8 ) and minimum support allow for better discovery of low frequency network (default=0.3).

The graph are created with association between genes, genes are node and link are represented with a line.

Different results files are created : 
- log.txt : contain all informations of script execution 
- association_rules.json : contain link between genes in json format 
- results.txt : contain groups of genes in txt format 
- graf file : representation of association between genes

Exemple of assocation rules in json format : 

```json
"GENE 1": {
        "support": 0.5698989412897016,
        "GENE 2": {
            "support": 0.5170837343599615,
            "confidence": 0.9073253113785096,
            "lift": 1.1536925176959112
        },
        "GENE 3": {
            "support": 0.5167228103946102,
            "confidence": 0.9066919991555837,
            "lift": 1.154653576985018
        } 
},
"GENE 2": {
        "support": 0.7864533205004812,
        "GENE 4": {
            "support": 0.4684793070259865,
            "confidence": 0.9553483807654564,
            "lift": 1.2147553527493458
        },
```

## Pre required : 

These packages are required to generate the network graph : 
- package networkx (https://networkx.github.io/)
- package pyvis (https://pyvis.readthedocs.io/en/latest/)

## Option 

This script has differents options to perform custom or more accurate analyses.


Select action to realize : 
- 'datamining' : Only computes association rules between gene 
- 'graph' : only draws a graph with the results of association rules
- 'full' : association rules + graph 

Select input file : 
- -i or --input : Name of input file. It should be a count matrix with cell in line and gene in column, in .tsv or .csv format. A JSON file of association rules can also be used in case of graph drawing.

example : 
```bash
python script_datamining.py full -i matrixcount.csv
```

optional arguments:

- --normalize            : Need to normalize matrix before analysis (default : false)

- --transpose           : Need to transpose matrix to have gene in column
                        (default : false)
- -o str, --output str  : Output directory (default name: results_analyse).
- -l int, --max-length int :   Max length of relations (default: 4).
- -s float, --min-support float :  Minimum support ratio (must be > 0, default: 0.3).
- -t float, --max-support float :   Maximum support ratio (must be < 1 , default: 0.8). In order to remove too express gene (who express in all     cells)
- -c float, --min-confidence float :  Minimum confidence (default: 0.9).
- -r str, --rowremove str : Remove row of matrix (allow to remove useless gene), for example : -r MT-gene1, MT-gene2
- -p int, --processor int : Number of processor available (improve compute time)


Exemple of command line : 
```bash
python3 script_datamining.py datamining -i countmatrix.tsv --transpose -p 8 -o resultat_supp0.1_0.5 -c 0.95 -s 0.1 -t 0.5 -l 4
```
