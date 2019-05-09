import re

import pandas as pd
import numpy as np
import xlmhg as hg
import scipy.stats as ss
import time
import math

FLOAT_PRECISION = 0.001
  
    



def combination_product(discrete_exp,c_list,coi,xlmhg):
    count = 0
    trips_list=[]
    for index,row in xlmhg.iterrows():
        if count == 50:
            break
        else:
            trips_list.append(row['gene_1'])
            count = count + 1


    def quads_matrix_gen(matrix):
        Cth = time.time()
        quads_matrix_list = []
        for row1 in np.transpose(matrix):
            for row2 in np.transpose(matrix):
                quads_matrix_list.append(row1&row2)
        quads_matrix = np.asarray(quads_matrix_list)
        print('')
        print('quads matrix times itself')
        print(quads_matrix.shape)
        print('')
        Ath = time.time()
        quads_product = np.matmul(quads_matrix,np.transpose(quads_matrix))
        #print('N^2 x N^2 full count matrix hopefully')
        #print(quads_product)
        #print(quads_product.shape)
        #print('')
        Bth = time.time()
        print(str(Bth-Ath) + ' mult seconds')
        return quads_product

    ###############
    #Experimental Section
    if trips_list == None:
        pass
    else:
        for column in discrete_exp:
            #print(column)
            if str(column) in trips_list:
                continue
            else:
                discrete_exp.drop(column, axis=1,inplace=True)
    ################
    
    
    start_time = time.time()
    #print('discrete exp matrix')
    #print(discrete_exp)
    in_cls_matrix = discrete_exp[c_list == coi].values
    #print('')
    #print('in cls matrix')
    #print(in_cls_matrix)
    #print('')
    total_matrix = discrete_exp.values
    #print('total matrix')
    #print(total_matrix)
    #print('')
    gene_count = len(discrete_exp.columns)
    first = time.time()
    quads_in_cls_product = quads_matrix_gen(in_cls_matrix)
    quads_total_product = quads_matrix_gen(total_matrix)

    gene_map = discrete_exp.columns.values
    #print(gene_map)
    odd_gene_map = []
    count = 0
    for gene in gene_map:
        for x in range(gene_count):
            odd_gene_map.append(gene)
    #print(odd_gene_map)
    #print('')
    odd_gene_map = pd.Index(odd_gene_map)
    even_gene_map = []
    count = 0
    for gene in gene_map:
        for x in gene_map:
            even_gene_map.append(x)
    #print(even_gene_map)
    even_gene_map = pd.Index(even_gene_map)

        
    quads_indices = np.triu_indices((gene_map.size*gene_map.size),1)
    return (
        quads_in_cls_product,quads_total_product,quads_indices,odd_gene_map,even_gene_map
        )





def quads_hg(gene_map,in_cls_count,pop_count,quads_in_cls,quads_total,quads_indices,odd_gene_map,even_gene_map):

    def tp(taken_in_cls):
        return taken_in_cls / in_cls_count

    def tn(taken_in_cls, taken_in_pop):
        return (
            ((pop_count - in_cls_count) - (taken_in_pop - taken_in_cls))
            / (pop_count - in_cls_count)
        )
    st = time.time()
    tp_result = np.vectorize(tp)(quads_in_cls[quads_indices])
    tn_result = np.vectorize(tn)(
        quads_in_cls[quads_indices], quads_total[quads_indices]
    )
    
    vhg = np.vectorize(ss.hypergeom.sf, excluded=[1, 2, 4], otypes=[np.float])

    
    hg_result = vhg(
        quads_in_cls[quads_indices],
        pop_count,
        in_cls_count,
        quads_total[quads_indices],
        loc=1
    )
    #0th index in quads_indices refers to 0th pair of genes
    
    #print(quads_indices[0])
    #print(quads_indices[1])
    #print(odd_gene_map)
    #print(even_gene_map)
    #print(gene_map)
    #print('HG + TP/TN done')
    output = pd.DataFrame({
        'gene_1': odd_gene_map[quads_indices[0]],
        'gene_2': even_gene_map[quads_indices[0]],
        'gene_3': odd_gene_map[quads_indices[1]],
        'gene_4': even_gene_map[quads_indices[1]],
        'HG_stat': hg_result,
        'TP' : tp_result,
        'TN': tn_result
    }, columns=['gene_1', 'gene_2', 'gene_3', 'gene_4', 'HG_stat','TP','TN'])

    en = time.time()
    print(str(en-st) + ' seconds')
    print('end HG/TP/TN')
    print('')
    print('begin filter')
    filt = time.time()
    output = output.sort_values(by='HG_stat', ascending=True)

    used_genes = []
    counter=0
    prev_stat = 0
    dropped=0
    for index, row in output.iterrows():
        #row[0] = gene1
        #row[1] = gene2
        #row[2] = gene3
        if counter == 1000:
            break
        #if row[0] != 'LY6D' or row[1] != 'CD3G_c':
        #    output.drop([index],inplace=True)
        #    dropped=dropped + 1
        #    print(dropped)
        #    continue
        #if row[-1] < .9:
        #    output.drop([index],inplace=True)
        #    continue
        if row[0]==row[1] or row[1]==row[2] or row[0]==row[2] or row[0]==row[3] or row[1]==row[3] or row[2]==row[3]:
            output.drop([index],inplace=True)
            continue
        if row[3] == prev_stat:
            output.drop([index],inplace=True)
            continue
        else:
            prev_stat = row[3]
        counter = counter+1

    endfilt = time.time()
    print(str(endfilt-filt) + ' seconds')
    print('end filter')
    return output
