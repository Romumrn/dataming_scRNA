#! usr/bin/python3


"""
A simple implementation of datamaining algorithm to scRNA seq usage.

"""


import pandas as pd
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Meta informations.
__version__ = '1'
__author__ = 'Romuald marin'
__author_email__ = 'romuald.mrn@outlook.fr'


def parse_args(argv):
    """
    Parse commandline arguments.

    """
    parser = argparse.ArgumentParser(description='Blabla help...' )
    parser.add_argument(
        'input', metavar='inpath',
        help='Name of input file in tsv. It should be a count matrix with cell in line and gene in collumn',
        type=str)
    parser.add_argument(
        '-n', '--normalize', metavar='bool',
        help='Matrice no normalize. (default : true)',
        type=bool, default=True)
    parser.add_argument(
        '-o', '--output', metavar='str',
        help='Output file (default name: resultat.txt).',
        type=str, default='resultat.txt')
    parser.add_argument(
        '-l', '--max-length', metavar='int',
        help='Max length of relations (default: infinite).',
        type=int, default=None)
    parser.add_argument(
        '-s', '--min-support', metavar='float',
        help='Minimum support ratio (must be > 0, default: 0.1).',
        type=float, default=0.1)
    parser.add_argument(
        '-c', '--min-confidence', metavar='float',
        help='Minimum confidence (default: 0.5).',
        type=float, default=0.5)
    parser.add_argument(
        '-r', '--rowremove', metavar='str',
        help="remove row of matrix for exemple : -r N_unmapped,N_multimapping,ASAT1 ",
        type=str, default='')
    args = parser.parse_args(argv)
    return args

def bool_matrix( df ):
    """
    Take a normalize count matrix and tranform into boolean matrix
    """
    #Remove collumn full of 0, that is mean gene doesnt expresded 
    df = df.loc[:, (df != 0).any(axis=0)]
    # if value 0 -> false
    df = ( df != 0 )
    return df

def bool_and_normalize_matrix( df ):
    """
    Take a NON-normalize count matrix and tranform into boolean matrix
    """
    df = df.loc[:, (df != 0).any(axis=0)]
    print('Normalise')
    df = df.applymap( np.log2 )
    df = df / df.max()
    df = (df > df.quantile(0.1) )
    return df

def apriori( data, minSupport, nb_transaction, output_file , max_len):
    print('Lauch new apriori')
    out = open("RESULTAT.txt", "w")
    current_lenght = 1
    print( "Generate C"+str(current_lenght))
    list_resultat_all =[]
    # 1st step : generate list of item who pass the threshold 
    resultat_C1 = new_generate_C1(data, minSupport,nb_transaction )
    
    list_candidat1 = []
    for item in resultat_C1:
        #add resultat to global results list
        list_resultat_all.append( item )
        #wirte in file 
        out.write( str(item)+'\n')
        #create a list of item who pass the threshold
        list_candidat1.append( item[0] )

    list_candidat1 =  list_candidat1[:1000] #LINE TEST to compute with x genes only !!!!!!!!!!!!!

    # 2nd step : Generate Candidat of lenght 2
    current_lenght = 2
    print( "Generate C"+str(current_lenght))
    res_candidat2 = []
    position_in_list = 0
    for i in range( position_in_list, len( list_candidat1 )):
        for j in range( position_in_list+1 , len( list_candidat1 )):
            support = len ( df[(df[ list_candidat1[i]]== True) | (df[list_candidat1[j]] == True)].index) / nb_transaction 
            if support > minSupport:
                new_result = [ (list_candidat1[i],list_candidat1[j]) , support ]
                res_candidat2.append( new_result  )
                list_resultat_all.append( new_result )
                # Write in out file
                out.write( str(new_result)+'\n')
        position_in_list += 1

    print( '    number of itemsets find : ' ,len(res_candidat2))
    assoc_rules = calc_assoc_rules(list_resultat_all)
    print_graph( assoc_rules )

    itemsets_prec = res_candidat2
    while current_lenght < max_len:
        current_lenght += 1
        print( "Generate C"+str(current_lenght))
        result = test_auto_apriori( data, list_candidat1, itemsets_prec, minSupport )
        for itemset in result:
            list_resultat_all.append( itemset )
            # Write in out file
            out.write( str(itemset)+'\n')
        itemsets_prec = result
    
    return list_resultat_all
    #fichier.close()

def new_generate_C1(data, minSupport, len_transaction):
    """
    Take count matrix and min support threshold to return Candidat K = 1 
    """
    #number of transaction to calc support 
    c1 = []
    for i in data.columns:
        support_values = data[i].sum() / len_transaction
        if support_values > minSupport :
            c1.append(  [i , support_values ] )
    return c1

def test_auto_apriori( data, C1, res_itemsets_prec, min_support):
    resultat = []
    #extrait la liste d'itemset de la liste itemset et resultat en tuple
    itemsets_prec = list(map(lambda x: list(x[0]), res_itemsets_prec) )
    print( '    number of itemsets find : ' ,len(itemsets_prec))
    for add_new_item in C1:
        for itemset in itemsets_prec:
            if add_new_item not in itemset:
                new_item = itemset+ [add_new_item]
                support = df[new_item].all(axis='columns').sum() / nb_transaction 
                if support > min_support:
                    new_res = [(new_item) , support ]
                    resultat.append( new_res )
            else:
                pass
    return resultat

def calc_assoc_rules( list_itemsets):
    """
    Take results variable and genereate asssociation rules between different genes
    confidence( X -> Y)= support(X , Y ) / support(X)
    """
    """
    exemple : 
    ['TRNE', 0.8428571428571429]
    ['TRNT', 0.8857142857142857]
    [('VWA1', 'CASZ1'), 0.37142857142857144]
    [('MMP23B', 'CDK11A'), 0.35714285714285715]
    [('MMP23B', 'CASZ1'), 0.37142857142857144]
    """
    dico_item = {}
    for item in list_itemsets:
        if type(item[0]) == str:
            dico_item[ item[0] ] = { 'support' : item[1] }
        else: 
            dico_item[ item[0][0] ][ item[0][1] ] = {  'support' : item[1] , 'confidence' : (item[1] / dico_item[ item[0][0] ]['support']) }  # : support(X , Y ) / support(X) } }
    return dico_item

def print_graph( asssociation_rules ):
    G = nx.Graph()
    G.support = {}

    node_list = []
    for gene in asssociation_rules.keys():
        for gene_link in asssociation_rules[gene].keys():
            if gene_link != 'support' :
                G.add_edge(gene, gene_link, weight=asssociation_rules[gene][gene_link]['confidence'] )

                node_list.extend( [ gene,gene_link ])

    #uniq element
    node_list = list(set(node_list)) 
    
    


    for gene in node_list:
        G.add_node( gene)
        G.support[gene] = asssociation_rules[gene]['support']*20
    #cmap information : https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    nx.draw( G , node_color=  [len( G.edges(n)) for n in G.nodes] ,edge_color= [d['weight']*20 for (u, v, d) in G.edges(data=True)], edge_cmap=plt.cm.Greys, node_cmap=plt.cm.YlGnBu, width = [d['weight']/2 for (u, v, d) in G.edges(data=True)], node_size=[G.support[n]*4 for n in G.nodes], with_labels=True , font_size=5)
    plt.show()

"""
Executes Apriori algorithm and print its result.
"""
args = parse_args(sys.argv[1:])
print(args)
#Import dataset in tsv
df = pd.read_csv(args.input, sep='\t', index_col=0)

#remove row 
for_removing = args.rowremove.split(',')
df = df.drop( for_removing ,axis = 1)

#transform into boolean matrix
if args.normalize : 
    matrix_bool = bool_and_normalize_matrix( df )
else :
    matrix_bool = bool_matrix( df )

print(matrix_bool.head())

#define max len 
if args.max_length:
    max_len = args.max_length
else:
    max_len = len(matrix_bool.columns )

print( "Maximun lenght of itemset is : ", max_len)

nb_transaction = len(matrix_bool.index)
res = apriori( matrix_bool, args.min_support, nb_transaction , args.output, max_len)