# Projet : Genes network with scRNA data

## Why ?

A gene regulatory networks represent interaction between biological chemicals (messenger RNAs (mRNAs) and proteins). These mRNA and proteins interact with each other with various degrees of specificity. In this case a nodes of this network represent genes and each arrow represent link between genes. The aim of this script is to compute easily a representation of gene network .


## Fonctionnement

This script use Apriori algorithm to create assocation rules between each genes. It proceeds by identifying the frequent individual gene in the gene count matrix (prealablement transformed into boolean all in removing non expressed gene and too small count ) and extending them to larger and larger group of genes. 

Apriori algorihm use 3 values : 
- Support : Fraction of cells who expressed genes (or a group of genes)

- Confidence :  Fraction of gene expressed if an other gene (or a group of genes) is expressed
- Lift : Measure of correlation, beetwen a gene and other gene (or groups of genes)

# Add explication therorique assocation rules et thressold json, graph 

## Script

This script allowed you to analyze single cell RNA data in order to find realtion between genes and create a network graph of gene expression. Based on the Apriori algorithm, basically used for frequent item set mining and association rule learning over relational databases. It proceeds by identifying the frequent individual items in the database and extending them to larger and larger item sets as long as those item sets appear sufficiently often in the database (in our case an item is a gene and a itemsets i a groupe of genes).

exemple of count matrix : 

|   | gene 1| gene 2   | gene 3   | gene 4  |
|---|---|---|---|---|
| Cell 1  |  0 |  1000 | 10000  | 0   |
| Cell 2  |  0 |  0 | 9000  |   0|
| Cell 3  |  0 |  1300 |  11000 |   0|


# Add explication sur le deroulement json, graph 




## Pre required : 

These packages are required to print nteworkgraph : 
- package networkx (https://networkx.github.io/)
- package pyvis (https://pyvis.readthedocs.io/en/latest/)

## Option 

This script have differents options to do more accurate analysis 


Select action to realize : 
- 'datamining' : Only doing accoatioation rules bewteen gene 
- 'graph' : only draw a graph with result of association rules
- 'full' : assocation rules + graph 

Select input file : 
- -i or --input : Name of input file. It should be a count matrix with cell in line and gene in collumn in .tsv or .csv or a JSON file of association rules in case of graph drawing

exemple : 
```bash
pyhton script_datamining.py full -i matrixcount.csv
```

optional arguments:

- --normalize            : Need to normalize matrix before analysis (default : false)

- --transpose           : Need to transpose matrix to have gene in collumn
                        (default : false)
- -o str, --output str  : Output directory (default name: results_analyse).
- -l int, --max-length int :   Max length of relations (default: 4).
- -s float, --min-support float :  Minimum support ratio (must be > 0, default: 0.3).
- -t float, --max-support float :   Maximum support ratio (must be < 1 , default: 0.8). In order to remove too express gene (who express in all     cells)
- -c float, --min-confidence float :  Minimum confidence (default: 0.9).
- -r str, --rowremove str : Remove row of matrix (allow to remove useless gene), for exemple : -r MT-gene1, MT-gene2
- -p int, --processor int : Number of processor avaible (improve compute time)




