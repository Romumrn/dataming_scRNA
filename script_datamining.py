#! usr/bin/python3

import pandas as pd
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from multiprocessing import Process, Queue, Pool

# Meta informations.
__version__ = '1'
__author__ = 'Romuald marin'
__author_email__ = 'romuald.mrn@outlook.fr'


def parse_args(argv):
    """
    Parse commandline arguments.
    For exemple : python3 script_datamining.py ../count_matrix -r N_unmapped,N_multimapping,N_noFeature,N_ambiguous -s 0.60 -l 5
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
    parser.add_argument(
        '-p', '--processor', metavar='int',
        help="Number of processor avaible (improve speed) ",
        type=int, default=1)
    parser.add_argument(
        '-g', '--graph',  action='store_true',
        help="Create  graph " )
    args = parser.parse_args(argv)
    return args

def bool_matrix( dataframe ):
    """
    Take a normalize count matrix and tranform into boolean matrix
    """
    #Remove collumn full of 0, that is mean gene doesnt expresded 
    dataframe = dataframe.loc[:, (dataframe != 0).any(axis=0)]
    # if value 0 -> false
    dataframe = ( dataframe != 0 )
    return dataframe

def bool_and_normalize_matrix( dataframe ):
    """
    Take a NON-normalize count matrix and tranform into boolean matrix
    """
    dataframe = dataframe.loc[:, (dataframe != 0).any(axis=0)]
    print('Normalise')
    dataframe = dataframe.applymap( np.log2 )
    dataframe = dataframe / dataframe.max()
    dataframe = (dataframe > dataframe.quantile(0.1) )
    return dataframe

def splitlist(li , n):
    """
    Take a list and return a new list compost of sub list with equal number of element when it's possible
    """
    newli = []
    div = int ( len(li) / n )
    i = 0 
    while i < (n -1 ) :
        subli = li[i*div : (i+1)*div]
        newli.append( subli)
        i += 1 
    #add rest of list
    subli = li[ (i*div) : ]
    newli.append( subli)
    return newli


def main_apriori( data, minSupport, nb_transaction, output_file , max_len, processor , graph):
    print('Lauch new apriori')
    out = open("RESULTAT.txt", "w")
    current_lenght = 1
    print( "Generate C"+str(current_lenght))
    list_resultat_all =[]

    # 1st step : generate list of item who pass the threshold 
    resultat_C1 = generate_C1(data, minSupport,nb_transaction , processor)
    
    list_candidat1 = []
    for item in resultat_C1:
        #add resultat to global results list
        list_resultat_all.append( item )
        #wirte in file 
        out.write( str(item)+'\n')
        #create a list of item who pass the threshold
        list_candidat1.append( item[0] )

    print( '    number of itemsets find : ' ,len(list_candidat1))


    list_candidat1 =  list_candidat1[:2000] #LINE TEST to compute with x genes only !!!!!!!!!!!!!
    print( '    new number of itemsets find : ' ,len(list_candidat1))

    itemsets_prec = list_candidat1
    while current_lenght < max_len:
        if itemsets_prec == []:
            print( 'No frequent itemsets avaible')
            quit()
        current_lenght += 1
        print( "Generate C"+str(current_lenght))
        result = do_apriori_multip( data, list_candidat1, itemsets_prec, minSupport , processor )
        print( '    number of itemsets find : ' ,len(result))
        for itemset in result:
            list_resultat_all.append( itemset )
            # Write in out file
            out.write( str(itemset)+'\n')
        if current_lenght == 2 and graph == True:
            print( "Draw ntework graph")
            assoc_rules = calc_assoc_rules(list_resultat_all)
            print_graph( assoc_rules )
            
        itemsets_prec = list(map(lambda x: list(x[0]), result) )
    
    return list_resultat_all


def calcsupport( data, item, totlen):
    support =  ( data[ item ].all(axis='columns').sum() ) / nb_transaction
    return support

def generate_C1(data, minSupport, len_transaction, proc):
    """
    Take count matrix and min support threshold to return Candidat lenght K = 1 
    """
    #number of transaction to calc support 
    c1 = []

    #create pool of process
    pool = Pool(processes=proc)

    list_item_split = splitlist ( data.columns , proc)

    var_pool = []
    for i in range(proc):
        print( 'proc. ',i, len(list_item_split[i]) )
        var_pool.append( ( list_item_split[i] , data, minSupport , len_transaction ) )
    
    #calcul 
    res = pool.map( processC1, var_pool ) 

    #concatenate lists of results into one list
    results = []
    for subres in res:
        print( len(subres))
        results = results+subres

    pool.close()
    
    return results

def processC1( args ):
    ( liitem, data, minSupport, len_transaction  ) = args
    res = []
    for i in liitem:
        support_values = data[i].sum() / len_transaction
        if support_values >= minSupport :
            res.append(  [ (i) , support_values ] )
    return res


def auto_apriori_process( args_list ):
    """
    Function process 
    """
    #unpack list of arg 
    (data, C1, itemsets_prec , min_support ) = args_list

    res = []

    for add_new_item in C1:
        for itemset in itemsets_prec:
            if add_new_item not in itemset:
                # if str that mean is the list of C1 for compute lengh 2
                if type(itemset) == str:
                    newitem = [itemset , add_new_item ]
                else:
                    newitem = itemset+ [add_new_item]
                support = calcsupport( data, newitem , nb_transaction )
                if support >= min_support:
                    newres = [(newitem) , support ]
                    res.append( newres )
            else:
                pass
    return res

def do_apriori_multip(data, C1, res_itemsets_prec, min_support, nb_process):
    #create pool of process
    pool = Pool(processes=nb_process)
    
    # create list of variable.
    list_var_pool = []

    #split list to give into process 
    itemsetsplit = splitlist( res_itemsets_prec , nb_process )

    for i in range(nb_process):
        list_var_pool.append( (data, C1, itemsetsplit[i] , min_support ) )
    
    # Use map - blocks until all processes are done.

    res = pool.map(auto_apriori_process, list_var_pool ) 
    #concatenate 4 list of results into one list
    results = []
    for i in res:
        results = results+i

    pool.close()

    return results


def calc_assoc_rules( list_itemsets):
    """
    Take results variable and genereate asssociation rules between differents genes
    confidence( X -> Y)= support(X , Y ) / support(X)
    """
    dico_item = {}
    for item in list_itemsets:
        if type(item[0]) == str:
            dico_item[ item[0] ] = { 'support' : item[1] }
        else: 
            dico_item[ item[0][0] ][ item[0][1] ] = {  'support' : item[1], 'confidence' : (item[1] / dico_item[ item[0][0] ]['support']), 
            'lift' : (item[1]/(dico_item[ item[0][0]]['support']* dico_item[ item[0][1] ]['support']) ) }  
            # confidence : s(X , Y ) / s(X) 
            # lift : S(X , Y ) / S(X)*S(Y)
    return dico_item

def print_graph( asssociation_rules ):
    G = nx.Graph()
    G.support = {}
    print( asssociation_rules )
    node_list = []
    count_threshold = 0
    for gene in asssociation_rules.keys():
        for gene_link in asssociation_rules[gene].keys():
            if gene_link != 'support' : #le support est enregrister comme une clé dans le dico
                if asssociation_rules[gene][gene_link]['lift'] >= 1 and asssociation_rules[gene][gene_link]['confidence'] > 0.95:
                    print( [ gene,gene_link ] )
                    G.add_edge(gene, gene_link, weight=asssociation_rules[gene][gene_link]['confidence'] )
                    node_list.extend( [ gene,gene_link ])
                else: 
                    count_threshold += 1
    print( "refused link : " , count_threshold)

    #uniq element
    node_list = list(set(node_list)) 
    
    for gene in node_list:
        G.add_node( gene)
        G.support[gene] = asssociation_rules[gene]['support']*20
    node_color=  [len( G.edges(n)) for n in G.nodes]
    edge_color= [d['weight']*20 for (u, v, d) in G.edges(data=True)]
    width =  [d['weight']/2 for (u, v, d) in G.edges(data=True)]
    nx.draw( G , node_color= node_color ,edge_color=edge_color, edge_cmap=plt.cm.gray, node_cmap=plt.cm.Spectral, width = width , node_size=[G.support[n]*4 for n in G.nodes], with_labels=True , font_size=5)
    plt.show()


# --------------------------------------------------------------------------------

if __name__ == "__main__":

    #parse arguments given
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


main_apriori( matrix_bool, args.min_support, nb_transaction , args.output, max_len, args.processor , args.graph)