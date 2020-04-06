#! usr/bin/python3


"""
A simple implementation of datamaining algorithm to scRNA seq usage.

"""


import pandas as pd
import argparse
import sys
import numpy as np

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

def apriori(data, minSupport, nb_transaction):
    print( 'go apriori')
    list_resultat_all =[]
    # 1st step : generate list of item 2 who pass the threshold 
    dico_candidat1 = generate_C1( data, minSupport,nb_transaction )
    list_candidat1 =  list(dico_candidat1.keys())[:10] #TEST
    res_candidat2 = []
    position_in_list = 0
    for i in range( position_in_list, len( list_candidat1 )):
        for j in range( position_in_list+1 , len( list_candidat1 )):
            res = len ( df[(df[ list_candidat1[i]]== True) | (df[list_candidat1[j]] == True)].index) / nb_transaction 
            if res > minSupport:
                res_candidat2.append( (list_candidat1[i],list_candidat1[j], res ) )
        position_in_list += 1
    
    res_candidat3 = generate_3(data, res_candidat2 , minSupport, nb_transaction)
    list_resultat_all = list(dico_candidat1)+res_candidat2+res_candidat3
    return list_resultat_all


def new_apriori( data, minSupport, nb_transaction, output_file , max_len):
    print( 'Lauch new apriori')
    out = open("RESULTAT.txt", "a")
    
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

    list_candidat1 =  list_candidat1[:100] #LINE TEST !!!!!!!!!!!!!

    # 2nd step : generate Candidat of lenght 2
    print( 'Generate C2')
    res_candidat2 = []
    position_in_list = 0
    for i in range( position_in_list, len( list_candidat1 )):
        for j in range( position_in_list+1 , len( list_candidat1 )):
            res = len ( df[(df[ list_candidat1[i]]== True) | (df[list_candidat1[j]] == True)].index) / nb_transaction 
            if res > minSupport:
                new_result = [ (list_candidat1[i],list_candidat1[j]) , res ]
                res_candidat2.append( new_result  )
                list_resultat_all.append( new_result )
        position_in_list += 1

    test_auto_apriori( data, list_candidat1, [], minSupport  )
    return list_resultat_all
    #fichier.close()

def new_generate_C1(data, minSupport, len_transaction):
    """
    Take count matrix and min support threshold to return Candidat K = 1 
    """
    print('generate C1')
    #number of transqction to calc support 
    c1 = []
    for i in data.columns:
        support_values = data[i].sum() / len_transaction
        if support_values > minSupport :
            c1.append(  [i , support_values ] )
    return c1

def test_auto_apriori( data, C1, itemsets_prec, min_support ):
    resultat = []
    for itemsets in itemsets_prec:
        for new_item in C1:
            if new_item in itemsets:
                pass
            else : 
                itemsets_to_test = itemsets.append(new_item)
                res = len ( df[(df[ itemsets_to_test == True ] )].index) / nb_transaction 
                if res > min_support:
                    pass
                    #blabla


def generate_3(data, res_precedent , minSupport, len_transaction):
    """
    But : a partir des candidat de taille precedente faire des nouveau candidat avec item(s) en commun 
    exemple : 
        DD2 DD3
        DD1 DD3
        --> DD2 DD1 DD3 
    faire selon la taille de la liste de candidat generé
    essayer de ranger les resultats dans des liste ou dico ou objet (?) 
        afin de rendre la lecture plus faciles ou generé un tableau a la fin (en créeant un fichier au fur et a mesure par exemple)
    """ 
    print( 'Generate C3')
    #generation de la liste de candidat a tester avec les resultats precedents
    #Quick and dirty need an optimization 
    list_candidat = []
    for i in range(0 , len(res_precedent)):
        for y in res_precedent[i+1:]:
            if res_precedent[i][0] in y:
                new_candidat = list(y[:-1]) 
                new_candidat.append( res_precedent[i][1] )
                list_candidat.append( new_candidat)

    #Calc item set de taille 3
    print('Calc C3')
    res_candidat3 = []
    position_in_list = 0
    for i in list_candidat:
        res = len ( df[(df[ i[0]]== True) | (df[ i[1]] == True) | (df[ i[2]] == True) ].index) / nb_transaction 
        if res > minSupport:
            res_candidat3.append( (i, res ) )
        position_in_list += 1
    return res_candidat3



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

#trasnform into boolean matrix
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
res = new_apriori( matrix_bool, args.min_support, nb_transaction , args.output, max_len)

for i in res:
    print(i)