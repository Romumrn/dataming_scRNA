#! usr/bin/python3

import pandas as pd
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool
import json
import os
import networkx as nx


# Meta informations.
__version__ = '1'
__author__ = 'Romuald marin'
__author_email__ = 'romuald.mrn@outlook.fr'


"""
 ArgumentParser.add_argument_group(title=None, description=None)
 ajouter des groupes pour une meilleur comprehension
"""

def parse_args(argv):
    """
    Parse commandline arguments.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''\
This scrit allow to realise a datamining analysis with apriori algortyhm. You can choose the number of process to improve the speed. 

You should give a count matrix of this shape : 
            Gene1   Gene2   Gene3   ... Gene
        Cell1   0      13   0   ... 0
        Cell2   0      0   0   ... 0 
        Cell3   2      13   0   ... 0    
        
For exemple : 
    python3 script_datamining.py full -i ../new_matrix/countmatrix.tsv  -s 0.90 -l 3 -p 4
    ''')

    parser.add_argument(
        'do', choices=['datamining', 'graph', 'full'],
        help='Select action to realize : datamanining analysis, only graph compute or full analysis',
        type=str)

    parser.add_argument(
        '-i','--input',
        help='Name of input file. It should be a count matrix with cell in line and gene in collumn in tsf or a JSON file of asso rules',
        type=str)

    parser.add_argument(
        '--normalize', action='store_true',
        help='Need to normalize matrix before analysis (default : false)', default=False)

    parser.add_argument(
        '--transpose', action='store_true',
        help='Need to transpose matrix to have gene in collumn (default : false)', default=False)

    parser.add_argument(
        '-o', '--output', metavar='str',
        help='Output directory (default name: results_analyse).',
        type=str, default='results_analyse')

    parser.add_argument(
        '-l', '--max-length', metavar='int',
        help='Max length of relations (default: infinite).',
        type=int, default=None)

    parser.add_argument(
        '-s', '--min-support', metavar='float',
        help='Minimum support ratio (must be > 0, default: 0.5).',
        type=float, default=0.5)

    parser.add_argument(
        '-t', '--max-support', metavar='float',
        help='Maximum support ratio (must be < 1 , default: 0.95). In order to remove gene who express in all cells',
        type=float, default=0.95)

    parser.add_argument(
        '-c', '--min-confidence', metavar='float',
        help='Minimum confidence (default: 0.9).',
        type=float, default=0.9)

    parser.add_argument(
        '-r', '--rowremove', metavar='str',
        help="remove row of matrix for exemple : -r N_unmapped,N_multimapping,ASAT1 ",
        type=str, default='')

    parser.add_argument(
        '-p', '--processor', metavar='int',
        help="Number of processor avaible (improve speed) ",
        type=int, default=1)

    args = parser.parse_args(argv)
    return args

def printandlog( openfile , *a):
    print( *a)
    for e in a:
        log.write( str(e)+'\n')


def bool_matrix( dataframe ):
    """
    Take a normalize count matrix and tranform into boolean matrix
    """
    #Remove collumn full of 0, that is mean gene doesnt expresded 
    dataframe = dataframe.loc[:, (dataframe != 0).any(axis=0)]
    # if value 0 -> false
    dataframe = ( dataframe != 0 )
    return dataframe

def readmtx( ):
    #https://www.biostars.org/p/312933/
    return 0

def bool_and_normalize_matrix( dataframe ):
    """
    Take a NON-normalize count matrix and tranform into boolean matrix
    """
    dataframe = dataframe.loc[:, (dataframe != 0).any(axis=0)]
    printandlog( log, 'Normalise Matrix')
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

def calcsupport( data, item, totlen):
    support =  ( data[ item ].all(axis='columns').sum() ) / nb_transaction
    return support

def generate_C1(data, minSupport, maxSupport, len_transaction, proc):
    """
    Take count matrix and min support threshold to return Candidat lenght K = 1 
    """

    #create pool of process
    pool = Pool(processes=proc)

    list_item_split = splitlist ( data.columns , proc)

    var_pool = []
    for i in range(proc): 
        var_pool.append( ( list_item_split[i] , data, minSupport ,maxSupport, len_transaction ) )
    
    #calcul 
    res = pool.map( processC1, var_pool ) 

    #concatenate lists of results into one list
    results = []
    for subres in res:
        results = results+subres

    pool.close()
    
    return results

def processC1( args ):
    ( liitem, data, minSupport,maxSupport, len_transaction  ) = args
    res = []
    for i in liitem:
        support_values = data[i].sum() / len_transaction
        if support_values >= minSupport :
            if support_values > maxSupport : 
                printandlog( log , "    refused :",i , support_values)
            else: 
                res.append(  [ (i) , support_values ] )

    return res


def auto_apriori_process( args_list ):
    """
    Function process 
    """
    #unpack list of arg 
    (data, C1, itemsets_prec , min_support,min_confidence ) = args_list

    res = []
    for add_new_item in C1:
        for itemset in itemsets_prec:
            itemsetname = itemset[0]

            if add_new_item not in itemsetname:
                # if str that mean is the list of C1 for compute lengh 2
                if type(itemsetname) == str:
                    newitem = [itemsetname , add_new_item ]
                else:
                    newitem = itemsetname+ [add_new_item]

                support = calcsupport( data, newitem , nb_transaction )
                confidence = support/ itemset[1]

                if support >= min_support and confidence >= min_confidence:
                    newres = [(newitem) , support, confidence]
                    res.append( newres )
            else:
                pass
    return res


def do_apriori_multip(data, C1, res_itemsets_prec, min_support, min_confidence, nb_process):
    #create pool of process
    pool = Pool(processes=nb_process)
    
    # create list of variable.
    list_var_pool = []

    #split list to give into process 
    itemsetsplit = splitlist( res_itemsets_prec , nb_process )

    for i in range(nb_process):
        list_var_pool.append( (data, C1, itemsetsplit[i] , min_support , min_confidence) )
    
    # Use map - blocks until all processes are done.

    res = pool.map(auto_apriori_process, list_var_pool ) 
    #concatenate 4 list of results into one list
    results = []
    for i in res:
        results = results+i

    pool.close()

    printandlog( log , "    Remove clone. Old number : ", len(results))
    results = uniqlist(results)
    printandlog( log , "        New number : ", len(results))
    return results

    """
    take a list of results and remove clone
    """

def uniqlist(li):
    li2 = list(map(lambda x: [sorted(x[0]), x[1], x[2]] , li))
    litest = []
    lireturn = []
    for i in li2:
        if str(i[0]) not in litest:
            litest.append( str(i[0]))
            lireturn.append( i )
    return lireturn

def calcAssociationRules( list_itemsets ):
    """
    Take results variable and genereate asssociation rules between differents genes
    confidence( X -> Y)= support(X , Y ) / support(X)
    """
    dico_item = {}
    for item in list_itemsets:
        if type(item[0]) == str:
            dico_item[ item[0] ] = { 'support' : item[1] }
        else: 
            dico_item[ item[0][0] ][ item[0][1] ] = {  
            'support' : item[1], 
            'confidence' : item[2] , 
            'lift' : item[1]/(dico_item[ item[0][0]]['support']* dico_item[ item[0][1] ]['support'])  }  
            # confidence : s(X , Y ) / s(X) 
            # lift : S(X , Y ) / S(X)*S(Y)
    
    return dico_item

from functools import reduce  # forward compatibility for Python 3
import operator

def getFromDict(dataDict, mapList):  
    for k in mapList: dataDict = dataDict[k]
    return dataDict

def setInDict(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value

#['THOC2','THOC3','THOC4','THOC9'], 0.9857142857142858
def add_association_rule( dicorules, listtoadd):
    for iteminlist in listtoadd:
        list_item = iteminlist[0]
        support = iteminlist[1]
        confidence = iteminlist[2]
        try:
            lift = support / getFromDict(dicorules, list_item[:-1]+['support']) * getFromDict(dicorules, [list_item[-1]]+['support']) 
            setInDict(dicorules, list_item + ['lift'] , lift )
        except:
            print( list_item, list_item[:-1]+['support'], [list_item[-1]]+['support'] )

        setInDict(dicorules, list_item + ['support'] , support )
        setInDict(dicorules, list_item + ['confidence'], confidence )

    return dicorules

def color_generator( numberfloat , limite ):
    # 0.8714285714285714 

    R = 0
    G = 500
    B = 100

    numberfloat = int ( (numberfloat - limite) * 6000 )
    for i in range(numberfloat ):
        if G < 220:
            G += 1
        else:
            B += 1

    return '#%02x%02x%02x' % (R,G,B)

def print_graph( asssociation_rules , path , minconf ):
    from pyvis.network import Network
    import networkx as nx

    G = Network(height="900px", width="100%", bgcolor="#222222", font_color="white")
    G.support = {}
    node_list = []
    count_threshold = 0
    for gene in asssociation_rules.keys():
        for gene_link in asssociation_rules[gene].keys():
            if gene_link != 'support' : #le support est enregrister comme une clé dans le dico
                #remove too hight support 
                """ if asssociation_rules[gene]['support'] > 0.9 or asssociation_rules[gene_link]['support'] > 0.9 :
                    #printandlog( log ,gene, asssociation_rules[gene]['support'], gene_link, asssociation_rules[gene_link]['support'])
                    continue"""
                if 'lift' in asssociation_rules[gene][gene_link].keys():
                    if asssociation_rules[gene][gene_link]['lift'] >= 1 and asssociation_rules[gene][gene_link]['confidence'] >= minconf :
                        if gene not in node_list and gene_link not in node_list:
                            G.add_node(gene, title=gene)
                            G.add_node(gene_link, title=gene_link)
                            G.support[gene] = asssociation_rules[gene]['support']
                            node_list.append( gene)
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                        elif gene in node_list and gene_link not in node_list:
                            G.add_node(gene_link, title=gene_link)
                            G.support[gene_link] = asssociation_rules[gene_link]['support']
                            node_list.append( gene_link)
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                        elif gene not in node_list and gene_link in node_list:
                            G.add_node(gene, title=gene)
                            G.support[gene] = asssociation_rules[gene_link]['support']
                            node_list.append( gene)
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                        else:
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                    else: 
                        count_threshold += 1
    printandlog( log , "refused link : " , count_threshold)
    
    neighbor_map = G.get_adj_list()

    # add neighbor data to node hover data
    for node in G.nodes:
        node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
        node["value"] = len(neighbor_map[node["id"]])

    #node_color= [node["value"] for n in G.nodes]
    #edge_color= [d['weight']*20 for (u, v, d) in G.edges(data=True)]
    #width =  [d['weight']/2 for (u, v, d) in G.edges(data=True)]
    #nx.draw( G , node_color= node_color , node_cmap=plt.cm.Spectral)
    #,edge_color=edge_color, edge_cmap=plt.cm.gray,, width = width , node_size=[G.support[n]*4 for n in G.nodes], with_labels=True , font_size=5 
    #G.show_buttons(filter_=['nodes', 'edges', 'physics'])

    G.set_options("""
    var options = {
  "edges": {
    "color": {
      "inherit": true
    },
    "smooth": false
  },
  "physics": {
    "forceAtlas2Based": {
      "gravitationalConstant": -70,
      "springLength": 100
    },
    "minVelocity": 0.75,
    "solver": "forceAtlas2Based"
  }
}
    """)
    G.show(path+"/graphSCRNA.html")


def print_graph2( asssociation_rules , path ):
    G = nx.Graph()
    G.support = {}
    node_list = []
    count_threshold = 0

    for gene in asssociation_rules.keys():
        for gene_link in asssociation_rules[gene].keys():
            if gene_link != 'support' : #le support est enregrister comme une clé dans le dico
                #remove too hight support 
                """ if asssociation_rules[gene]['support'] > 0.9 or asssociation_rules[gene_link]['support'] > 0.9 :
                    #printandlog( log ,gene, asssociation_rules[gene]['support'], gene_link, asssociation_rules[gene_link]['support'])
                    continue"""
                if 'lift' in asssociation_rules[gene][gene_link].keys():
                    if asssociation_rules[gene][gene_link]['lift'] >= 1 :
                        if gene not in node_list and gene_link not in node_list:
                            G.add_node(gene, title=gene)
                            G.support[gene] = asssociation_rules[gene]['support']*20
                            G.add_node(gene_link, title=gene_link)
                            G.support[gene_link] = asssociation_rules[gene_link]['support']*20
                            node_list.append( gene_link)
                            node_list.append( gene)
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                        elif gene in node_list and gene_link not in node_list:
                            G.add_node(gene_link, title=gene_link)
                            G.support[gene_link] = asssociation_rules[gene_link]['support']*20
                            node_list.append( gene_link)
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                        elif gene not in node_list and gene_link in node_list:
                            G.add_node(gene, title=gene)
                            G.support[gene] = asssociation_rules[gene_link]['support']*20
                            node_list.append( gene)
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                        else:
                            G.add_edge(gene, gene_link, weight = asssociation_rules[gene][gene_link]['lift'])
                    else: 
                        count_threshold += 1

    print( "refused link : " , count_threshold)

    #uniq element
    node_list = list(set(node_list)) 
    
    node_color=  [len( G.edges(n)) for n in G.nodes]
    edge_color= [d['weight']*20 for (u, v, d) in G.edges(data=True)]
    width =  [d['weight']/2 for (u, v, d) in G.edges(data=True)]
    nx.draw( G , node_color= node_color ,edge_color=edge_color, edge_cmap=plt.cm.gray, node_cmap=plt.cm.Spectral, width = width , node_size=[G.support[n]*4 for n in G.nodes], with_labels=True , font_size=5)
    plotpath = path+"/graphSCRNA.pdf"
    plt.savefig(plotpath)
    #plt.show()


def print_graph_profondeur( asssociation_rules , path ):
    G = nx.Graph()
    G.support = {}
    resultatentuple =[]
    count = 0 
    dicoedges = {}
    removecount =0
    for gene in asssociation_rules.keys():
        resultatentuple =[]
        #print( count,' / ', len(asssociation_rules.keys()))
        parcourir( asssociation_rules, [gene], resultatentuple)
        for i in resultatentuple:
            for node in i[0]:
                #print(i)
                if i[3] > 1:
                    G.add_node(node, title=node)
                    G.support[node] = asssociation_rules[node]['support']*20
                    for node2 in i[0]:
                        if node == node2:
                            pass
                        else:
                            keyyy = str(sorted([node, node2]))
                            if keyyy in dicoedges.keys():
                                dicoedges[ keyyy]['poid'] = dicoedges[ keyyy]['poid'] +1
                            else:
                                dicoedges[ keyyy ] = {'node1': sorted([node, node2])[0], 'node2':sorted([node, node2])[1] , 'poid': 1}
                else:
                    removecount += 1
        count+=1
    print( 'remove :',removecount)
    maxrelation = 0
    for edge in dicoedges:
        G.add_edge( dicoedges[edge]['node1'],dicoedges[edge]['node2'], weight = dicoedges[edge]['poid'])
        if maxrelation < dicoedges[edge]['poid']:
            maxrelation = dicoedges[edge]['poid']


    node_color= [G.support[n]*4 for n in G.nodes]
    edge_color= [d['weight']*120/maxrelation for (u, v, d) in G.edges(data=True)]
    width =  [d['weight']/maxrelation for (u, v, d) in G.edges(data=True)]
    #nx.draw( G , node_color= node_color ,edge_color=edge_color, edge_cmap=plt.cm.binary, node_cmap=plt.cm.Spectral, width = width , node_size=[len( G.edges(n)) for n in G.nodes], with_labels=True , font_size=3)

    #nx.draw( G , node_color= node_color ,node_cmap=plt.cm.Spectral,with_labels=True , font_size=5, width = width )

    plotpath = path+"/graphSCRNAdeep.pdf"
    plt.savefig(plotpath)
    plt.show()


def parcourir( dico , list_gene, resultatentuple):
    listkgu = list_gene
    #forme list [gene1, gene2 ..., support, confidence]
    #print(listkgu)
    list_link = [] #A chaque fois la fonction crée une nouvelle liste !!!!!!!!!!!!!!
    dicornot = True
    while dicornot:
        clef_gene = list( getFromDict( dico, listkgu).keys())
        #print( clef_gene)
        while len(clef_gene) > 0:
            new_gene = clef_gene[0]
            #print(clef_gene)
            if new_gene in ['lift','confidence','support']:
                if clef_gene == ['support','confidence','lift']:
                    #print( 'en cours',listkgu)
                    conf = getFromDict(dico, listkgu+['confidence'])
                    supp = getFromDict(dico, listkgu+['support'])
                    lift = getFromDict(dico, listkgu+['lift'])
                    return (listkgu , supp, conf,lift)
            else:
                #print('in ', listkgu+[new_gene])
                res = parcourir( dico, listkgu+[new_gene] , resultatentuple)
                if type(res) is tuple :
                    resultatentuple.append(res)
            clef_gene = clef_gene[1:]
        dicornot = False
    return list_link


# --------------------------------------------------------------------------------

if __name__ == "__main__":

    
    #parse arguments given
    args = parse_args(sys.argv[1:])

    path = args.output

    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

    log = open( path+'/log.txt', "w")
    printandlog( log ,args)

    if args.do == 'full' or args.do == 'datamining':
    
        matrixfile = args.input
        printandlog( log , 'Loading data from file '+matrixfile)
        if matrixfile.endswith('tsv'):
            #Import dataset in tsv
            df = pd.read_csv( matrixfile, sep='\t', index_col=0)
        elif matrixfile.endswith('csv'):
            #Import dataset in tsv
            df = pd.read_csv( matrixfile, sep=',', index_col=0)
        else:
            df = pd.read_csv( matrixfile, sep='\t', index_col=0)

        if args.transpose :
            printandlog( log , 'Transpose matrix')
            df = df.transpose()
        #remove row 
        if args.rowremove:
            printandlog( log , 'Remove row', args.rowremove)
            log.write( args)
            for_removing = args.rowremove.split(',')
            df = df.drop( for_removing ,axis = 1)

        #transform into boolean matrix
        if args.normalize : 
            matrix_bool = bool_and_normalize_matrix( df )
        else :
            matrix_bool = bool_matrix( df )

        printandlog( log ,matrix_bool.head())

        #define max len 
        if args.max_length:
            max_len = args.max_length
        else:
            max_len = len(matrix_bool.columns )

        printandlog( log , "Maximun lenght of itemset is : ", max_len)
        nb_transaction = len(matrix_bool.index)

        minSupport = args.min_support
        maxSupport = args.max_support
        minConf = args.min_confidence
        processor = args.processor

        printandlog( log ,'Lauch new apriori')
        out = open( path+'/results.txt', "w")
        current_lenght = 1
        printandlog( log , "Generate C"+str(current_lenght))
        list_resultat_all =[]

        # 1st step : generate list of item who pass the threshold 
        resultat_C1 = generate_C1(matrix_bool, minSupport, maxSupport, nb_transaction , processor)
        
        C1 = []
        list_candidat1 = []
        for item in resultat_C1:
            #add resultat to global results list
            list_resultat_all.append( item )
            #wirte in file 
            out.write( str(item)+'\n')
            #create a list of item who pass the threshold
            list_candidat1.append( item )
            C1.append( item[0])

        printandlog( log , '    number of itemsets find : ' ,len(list_candidat1))

        #list_candidat1 =  list_candidat1[:1000] #LINE TEST to compute with x genes only !!!!!!!!!!!!!
        printandlog( log , '    new number of itemsets find : ' ,len(list_candidat1))

        itemsets_prec = list_candidat1

        while current_lenght < max_len:
            if itemsets_prec == []:
                printandlog( log , 'No frequent itemsets avaible')
                break
            current_lenght += 1
            printandlog( log , "Generate C"+str(current_lenght))
            result = do_apriori_multip( matrix_bool, C1, itemsets_prec, minSupport , minConf, processor )
            printandlog( log , '    number of itemsets find : ' ,len(result))
            for itemset in result:
                list_resultat_all.append( itemset )
                # Write in out file
                out.write( str(itemset)+'\n')

            if current_lenght == 2 :
                assoc_rules = calcAssociationRules(list_resultat_all)
                #Stock into pickles obj dictionnary (can be use later to build a graph)
                with open(path+'/association_rules.json', 'w') as filejson:
                    json.dump(assoc_rules, filejson)
                if args.do == 'full':
                    printandlog( log , "Draw network graph")  
                    #print_graph( assoc_rules , path, args.min_confidence)
                    print_graph2( assoc_rules , path)

            if current_lenght > 2 :
                assoc_rules = add_association_rule(assoc_rules, result)
                with open(path+'/association_rules.json', 'w') as filejson:
                    json.dump(assoc_rules, filejson)

            itemsets_prec = result
        print_graph_profondeur( assoc_rules , path)
    
    elif args.do == 'graph':
        printandlog( log , "Draw network graph")
        with open( args.input, 'r') as f:
            association_rules = json.load(f)
        #print_graph( association_rules , path, args.min_confidence )
        #print_graph2( association_rules , path)
        print_graph_profondeur( association_rules , path)