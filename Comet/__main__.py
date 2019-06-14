import os
import argparse
import datetime
import re

import pandas as pd
import numpy as np

from . import hgmd
from . import visualize as vis
from . import quads
import sys
import multiprocessing
import time
import math
import matplotlib.pyplot as plt
import random
import scipy.io
import csv
#from docs.source import conf


def init_parser(parser):
    """Initialize parser args."""
    parser.add_argument(
        'marker', type=str,
        help=("Marker file input")
    )
    parser.add_argument(
        'tsne', type=str,
        help=("tsne file input")
    )
    parser.add_argument(
        'cluster', type=str,
        help=("Cluster file input")
    )
    parser.add_argument(
        '-g', nargs='?', default=None,
        help="Optional Gene list"
    )
    parser.add_argument(
        'output_path', type=str,
        help="the output directory where output files should go"
    )
    parser.add_argument(
        '-C', nargs='?', default=None,
        help="Num of cores avail for parallelization"
    )
    parser.add_argument(
        '-X', nargs='?', default=None,
        help="X argument for XL-mHG"
    )
    parser.add_argument(
        '-L', nargs='?', default=None,
        help="L argument for XL-mHG"
    )
    parser.add_argument(
        '-Abbrev', nargs='?',default=[],
        help="Choose between abbreviated or full 3-gene computation"
    )
    parser.add_argument(
        '-K', nargs='?',default=None,
        help="K-gene combinations to include"
    )
    parser.add_argument(
        '-Down', nargs='?',default=False,
        help="Downsample"
    )
    parser.add_argument(
        '-Trim', nargs='?',default=2000,
        help="Trim output files"
    )
    parser.add_argument(
        '-Count', nargs='?',default=False,
        help="Set to True when count data is being used, for visualizations."
    )
    parser.add_argument(
        '-tenx', nargs='?',default=False,
        help="Set to True when count data is being used, for visualizations."
    )
    parser.add_argument(
        '-online', nargs='?',default=False,
        help="Set to True for online version."
    )
    return parser


def read_data(cls_path, tsne_path, marker_path, gene_path, D, tenx, online):
    """
    Reads in cluster series, tsne data, marker expression without complements
    at given paths.
    """
    
    cls_ser = pd.read_csv(
        cls_path, sep='\t', index_col=0, names=['cell', 'cluster'], squeeze=True
    )
    if np.isnan(cls_ser[0]):
        cls_ser = pd.read_csv(
        cls_path, sep=',', index_col=0, names=['cell', 'cluster'], squeeze=True )

    tsne = pd.read_csv(
        tsne_path, sep='\t', index_col=0, names=['cell', 'tSNE_1', 'tSNE_2']
    )
    if np.isnan(tsne['tSNE_1'][0]):
        tsne = pd.read_csv(
        tsne_path, sep=',', index_col=0, names=['cell', 'tSNE_1', 'tSNE_2'] )

    
    start_= time.time()
    tenx = int(tenx)
    if tenx == 1:
        print('Loading 10X matrix')
        mat = scipy.io.mmread(marker_path+"matrix.mtx")
        features_path = marker_path + "genes.tsv"
        gene_names = [row[1] for row in csv.reader(open(features_path), delimiter="\t")]
        barcodes_path = marker_path + "barcodes.tsv"
        barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
        #construct pandas dataframe w/ the pieces (gene names, barcodes, counts in sparse form)
        matrix = pd.DataFrame(index = gene_names, columns = barcodes )
        now = time.time()
        print('assembling expression matrix')
        for i,j,v in zip(mat.row, mat.col, mat.data):
            matrix.iat[i,j] = v
        matrix.fillna(0,inplace=True)
        noww = time.time()
        print(str(noww-now) + ' seconds')
        print('size: ' + str(matrix.shape))
        no_complement_marker_exp = matrix
        no_complement_marker_exp.rename_axis('cell',axis=1,inplace=True)
    else:

        #Should allow either tab OR comma delimited formats
        try:
            no_complement_marker_exp = pd.read_csv(
                marker_path,sep='\t', index_col=0
                ).rename_axis('cell',axis=1)
            if len(no_complement_marker_exp.columns) == 0:
                raise Exception
            elif len(no_complement_marker_exp.index) == 0:
                raise Exception
            else:
                pass
        except:
            no_complement_marker_exp = pd.read_csv(
                marker_path,sep=',', index_col=0
                ).rename_axis('cell',axis=1)
           
    if no_complement_marker_exp.shape[1] == cls_ser.shape[0]:
        pass
    else:
        for index,row in cls_ser.iteritems():
            if str(index) in list(no_complement_marker_exp):
                continue
            else:
                cls_ser.drop(labels=index,inplace=True)
    
    #gene list filtering
    no_complement_marker_exp = no_complement_marker_exp.loc[~no_complement_marker_exp.index.duplicated(keep='first')]
    no_complement_marker_exp = np.transpose(no_complement_marker_exp)
    no_complement_marker_exp.columns = [x.upper() for x in no_complement_marker_exp.columns]
    no_complement_marker_exp = no_complement_marker_exp.loc[~no_complement_marker_exp.index.duplicated(keep='first')]
    #gene filtering
    #-------------#
    if gene_path is None:
        pass
    else:

        #read the genes
        #Compatible with single line comma list OR one per line no commas OR mix of both
        master_gene_list = []
        with open(gene_path, "r") as genes:
            lines = genes.readlines()
            if len(lines) == 1:
                with open(gene_path, "r") as genes:
                    init_read = genes.read().splitlines()
                    master_str = str.upper(init_read[0])
                    master_gene_list = master_str.split(",")
            else:
                for i, line in enumerate(lines):
                    if '\n' in line:
                        master_gene_list.append(line[:-1])
                    else:
                        master_gene_list.append(line)
                for item in master_gene_list[:]:
                    if ',' in item:
                        new_split = item.split(",")
                        master_gene_list.remove(item)
                        for ele in new_split:
                            master_gene_list.append(str.upper(ele))

    
        new_no_comp_mark_exp = pd.DataFrame()
        for gene in master_gene_list:
            try:
                new_no_comp_mark_exp[gene] = no_complement_marker_exp[gene]
            except:
                pass

        no_complement_marker_exp = new_no_comp_mark_exp
        '''
        for column_name in no_complement_marker_exp.columns:
            if str.upper(column_name) in master_gene_list:
                pass
            else:
                try:
                    no_complement_marker_exp.drop(column_name, axis=1, inplace=True)
                except:
                    pass
        '''
    #-------------#

    #downsampling
    #-------------#
    
    if D is False:
        pass
    else:
        # get number of genes to set downsample threshold
        gene_numb = len(no_complement_marker_exp.columns)
        #print(gene_numb)

        if gene_numb > 3000:
            if int(D) < int(2500):
                pass
            else:
                D = int(2500)
        #total number of cells input
        N = len(cls_ser)
        #print(N)
        #downsample target
        M = int(D)
        if N <= M:
            return (cls_ser, tsne, no_complement_marker_exp, gene_path)
        clusters = sorted(cls_ser.unique())
        counts = { x : 0 for x in clusters}
        for clus in cls_ser:
            counts[clus] = counts[clus]+1
        #at this point counts has values for # cells in cls
        #dict goes like ->{ cluster:#cells }
        take_nums = {x : 0 for x in clusters}
        for clstr in take_nums:
            take_nums[clstr] = math.ceil(counts[clstr]*(M/N))
        summ = 0
        for key in take_nums:
            summ = summ + take_nums[key]
        #print('Downsampled cell num ' + str(summ))
        counts= { x : 0 for x in clusters}
        new_cls_ser = cls_ser.copy(deep=True)
        keep_first = 0
        for index,value in new_cls_ser.iteritems():
            keep_first = keep_first + 1
            if keep_first ==1:
                counts[value] = counts[value]+1
                continue
            new_cls_ser.drop(index,inplace=True)
        cls_ser.drop(cls_ser.index[0],inplace=True)
        #Now new_cls_ser has all removed except first item, which we can keep
        for num in range(N-1):
            init_rand_num = random.randint(0,N-num-1-1)
            if counts[cls_ser[init_rand_num]] >= take_nums[cls_ser[init_rand_num]]:
                cls_ser.drop(cls_ser.index[init_rand_num], inplace=True)
                continue
            new_cls_ser = new_cls_ser.append(pd.Series([cls_ser[init_rand_num]], index=[cls_ser.index[init_rand_num]]))
            counts[cls_ser[init_rand_num]] = counts[cls_ser[init_rand_num]]+1
            cls_ser.drop(cls_ser.index[init_rand_num], inplace=True)
            
        new_cls_ser.rename_axis('cell',inplace=True)
        new_cls_ser.rename('cluster', inplace=True)
        return(new_cls_ser,tsne,no_complement_marker_exp,gene_path)


    #-------------#

            
    return (cls_ser, tsne, no_complement_marker_exp, gene_path)

def process(cls,X,L,plot_pages,cls_ser,tsne,marker_exp,gene_file,csv_path,vis_path,pickle_path,cluster_number,K,abbrev,cluster_overall,Trim,count_data):
    #for cls in clusters:
    # To understand the flow of this section, read the print statements.
    heur_limit = min(50,len(marker_exp.columns))
    start_cls_time = time.time()
    print('########\n# Processing cluster ' + str(cls) + '...\n########')
    print(str(K) + ' gene combinations')
    if len(abbrev) == 0:
        pass
    else:
        print('Heuristic limit set to: ' + str(heur_limit))
    print('Running t test on singletons...')
    try:
        t_test = hgmd.batch_stats(marker_exp, cls_ser, cls)
    except Exception as err:
        print('t-test error')
        print(err)
    print('Calculating fold change')
    try:
        fc_test = hgmd.batch_fold_change(marker_exp, cls_ser, cls)
    except Exception as err:
        print('fold-change error')
        print(err)
    print('Running XL-mHG on singletons...')
    try:
        xlmhg = hgmd.batch_xlmhg(marker_exp, cls_ser, cls, X=X, L=L)
    except Exception as err:
        print('XLMHG error')
        print(err)
    try:
        q_val = hgmd.batch_q(xlmhg)
    except Exception as err:
        print('q-val error')
        print(err)
    # We need to slide the cutoff indices before using them,
    # to be sure they can be used in the real world. See hgmd.mhg_slide()
    try:
        cutoff_value = hgmd.mhg_cutoff_value(
            marker_exp, xlmhg[['gene_1', 'mHG_cutoff']]
            )
        xlmhg = xlmhg[['gene_1', 'mHG_stat', 'mHG_pval']].merge(
            hgmd.mhg_slide(marker_exp, cutoff_value), on='gene_1'
            )
        # Update cutoff_value after sliding
        cutoff_value = pd.Series(
            xlmhg['cutoff_val'].values, index=xlmhg['gene_1']
            )
        xlmhg = xlmhg\
          .sort_values(by='mHG_stat', ascending=True)
    except Exception as err:
        print('error in sliding values')
        print(err)
    print('Creating discrete expression matrix...')
    try:
        discrete_exp = hgmd.discrete_exp(marker_exp, cutoff_value, abbrev, xlmhg)
    except Exception as err:
        print('discrete matrix construction failed')
        print(err)
    '''
    #For checking the sliding issue
    count = 0
    print(discrete_exp['Reg4'].sort_values(ascending=False).head(667))
    #time.sleep(100000)
    print(cls_ser)
    for index in discrete_exp['Reg4'].sort_values(ascending=False).head(667).iteritems():
        for index2 in cls_ser.iteritems():
            if index[0] == index2[0]:
                if index2[1] == 0:
                    count = count +1
    print(count)
    print(discrete_exp['Reg4'].sort_values(ascending=False).head(70))
    print(cls_ser.to_string())
    #print(marker_exp['1600029D21Rik'].sort_values(ascending=False).head(160))
    #time.sleep(100000)
    '''
    discrete_exp_full = discrete_exp.copy()
    print('Finding simple true positives/negatives for singletons...')
    #Gives us the singleton TP/TNs for COI and for rest of clusters
    #COI is just a DF, rest of clusters are a dict of DFs
    (sing_tp_tn, other_sing_tp_tn) = hgmd.tp_tn(discrete_exp, cls_ser, cls, cluster_overall)

    ### Take out any genes with a true positive less than 15% from the expression matrix ###
    for gene in discrete_exp:
        if sing_tp_tn.set_index('gene_1').at[gene,'TP'] <= .15:
            discrete_exp.drop(labels=gene,axis=1,inplace=True)
    ########################################################################################
    ###########
    #OLD HEURISTICS
    #abb = '3'
    #if abb in abbrev:
    #    print('Heuristic Abbreviation initiated for ' + str(abbrev) )
    #    count = 0
    #    trips_list=[]
    #    for index,row in xlmhg.iterrows():
    #        if sing_tp_tn.set_index('gene_1').at[row['gene_1'],'TP'] <= .15:
    #            continue
    #        ##THREE GENE VARIABLE HERE##
    #        if count == trips_heur:
    #            break
    #        else:
    #            trips_list.append(row['gene_1'])
    #            count = count + 1
    #else:
    #    trips_list = None
    ############
    print('Finding pair expression matrix...')
    (
        gene_map, in_cls_count, pop_count,
        in_cls_product, total_product, upper_tri_indices,
        cluster_exp_matrices, cls_counts
    ) = hgmd.pair_product(discrete_exp, cls_ser, cls, cluster_number,cluster_overall)
    if K >= 4:
        print('')
        print('Starting quads')
        print('')
        quads_in_cls, quads_total, quads_indices, odd_gene_mapped, even_gene_mapped = quads.combination_product(discrete_exp,cls_ser,cls,xlmhg)
        print('')
        print('')
        print('')
        print('')
        print('HG TEST ON QUADS')
        quads_fin = quads.quads_hg(gene_map,in_cls_count,pop_count,quads_in_cls,quads_total,quads_indices,odd_gene_mapped,even_gene_mapped)

    if K == 3:
        start_trips = time.time()
        print('Finding Trips expression matrix...')
        try:
            trips_in_cls,trips_total,trips_indices,gene_1_mapped,gene_2_mapped,gene_3_mapped = hgmd.combination_product(discrete_exp,cls_ser,cls,abbrev,heur_limit)
        except Exception as err:
            print('error in 3-gene matrix construction')
            print(err)
        end_trips = time.time()
        print(str(end_trips-start_trips) + ' seconds')

    HG_start = time.time()
    print('Running hypergeometric test on pairs...')
    pair, revised_indices = hgmd.pair_hg(
        gene_map, in_cls_count, pop_count,
        in_cls_product, total_product, upper_tri_indices, abbrev, heur_limit
    )
    pair_q = hgmd.pairs_q(pair)
    
    HG_end = time.time()
    print(str(HG_end-HG_start) + ' seconds')
    pair_out_initial = pair\
    .sort_values(by='HG_pval', ascending=True)
    
    pair_out_initial['rank'] = pair_out_initial.reset_index().index + 1
    #spair_out_print = pair_out_initial.head(Trim)
    pair_out_initial.to_csv(
    csv_path + '/cluster_' + str(cls) + '_pair_full_unranked.csv'
    )
    if K == 3:
        HG_start = time.time()
        print('Running hypergeometric test & TP/TN on trips...')
        trips = hgmd.trips_hg(
            gene_map,in_cls_count,pop_count,
            trips_in_cls,trips_total,trips_indices,
            gene_1_mapped,gene_2_mapped,gene_3_mapped,abbrev,heur_limit
            )
        #print(trips)
        HG_end = time.time()
        print(str(HG_end-HG_start) + ' seconds')


    # Pair TP/TN FOR THIS CLUSTER
    print('Finding simple true positives/negatives for pairs...')
    pair_tp_tn = hgmd.pair_tp_tn(
        gene_map, in_cls_count, pop_count,
        in_cls_product, total_product, upper_tri_indices, abbrev, revised_indices
    )
    #accumulates pair TP/TN vals for all other clusters
    ##NEW
    other_pair_tp_tn = {}
    for key in cluster_exp_matrices:
        new_pair_tp_tn = hgmd.pair_tp_tn(
            gene_map, cls_counts[key], pop_count,
            cluster_exp_matrices[key], total_product, upper_tri_indices,
            abbrev, revised_indices
        )
        other_pair_tp_tn[key] = new_pair_tp_tn
        other_pair_tp_tn[key].set_index(['gene_1','gene_2'],inplace=True)
    pair = pair\
        .merge(pair_tp_tn, on=['gene_1', 'gene_2'], how='left')\
        .merge(pair_q, on=['gene_1','gene_2'], how='left')
        
    #Supplementary figure 2 code#
    '''
    for cl in other_pair_tp_tn:
        print(cl)
        print(other_pair_tp_tn[cl].loc[('LY6D', 'CD3G_negation')])
    for cl in other_sing_tp_tn:
        print(cl)
        print(other_sing_tp_tn[cl].loc[('LY6D')])
    time.sleep(10000)
    '''
    pair_tp_tn.set_index(['gene_1','gene_2'],inplace=True)
    sing_tp_tn.set_index(['gene_1'], inplace=True)
    rank_start = time.time()
    print('Finding NEW Rank')
    ranked_pair,histogram = hgmd.ranker(pair,xlmhg,sing_tp_tn,other_sing_tp_tn,other_pair_tp_tn,cls_counts,in_cls_count,pop_count)
    rank_end = time.time()
    print(str(rank_end - rank_start) + ' seconds')
    # Save TP/TN values to be used for non-cluster-specific things
    print('Pickling data for later...')
    sing_tp_tn.to_pickle(pickle_path + 'sing_tp_tn_' + str(cls))
    pair_tp_tn.to_pickle(pickle_path + 'pair_tp_tn_' + str(cls))
    #trips_tp_tn.to_pickle(pickle_path + 'trips_tp_tn' + str(cls))
    print('Exporting cluster ' + str(cls) + ' output to CSV...')
    try:
        sing_output = xlmhg\
          .merge(t_test, on='gene_1')\
          .merge(fc_test, on='gene_1')\
          .merge(sing_tp_tn, on='gene_1')\
          .merge(q_val, on='gene_1')\
          .set_index('gene_1')\
          .sort_values(by='mHG_stat', ascending=True)
    except Exception as err:
        print(err)
        sing_output = xlmhg.sort_values(by='mHG_stat',ascending=True)
    sing_output.sort_values(by='gene_1',ascending=True).to_csv(
        csv_path + '/cluster_' + str(cls) + '_singleton_full_unranked.csv'
    )
        
    sing_output = sing_output.loc[sing_output['TP'] >= .15]

    for index, row in sing_output.iterrows():
        if re.compile(".*_negation$").match(index):
            if row['TN'] >= .5:
                pass
            else:
                sing_output.drop(labels=index,axis=0,inplace=True)


    # make a copy of sing_output
    # run the stuff below on both the original (fine) and the copy (no negations)
    sing_output_positives = sing_output.copy(deep=True)
    for index, row in sing_output_positives.iterrows():
        if re.compile(".*_negation$").match(index):
            sing_output_positives.drop(labels=index,axis=0,inplace=True)
    sing_output_positives['hgrank'] = sing_output_positives.reset_index().index + 1
    sing_output_positives.sort_values(by='Log2FoldChangeAbs', ascending=False, inplace=True)
    sing_output_positives['fcrank'] = sing_output_positives.reset_index().index + 1
    sing_output_positives['finrank'] = sing_output_positives[['hgrank', 'fcrank']].mean(axis=1)
    sing_output_positives.sort_values(by='finrank',ascending=True,inplace=True)
    sing_output_positives['rank'] = sing_output_positives.reset_index().index + 1
    sing_output_positives.drop('finrank',axis=1, inplace=True)
    count = 1
    for index,row in sing_output_positives.iterrows():
        if count == 100:
            break
        sing_output_positives.at[index,'Plot'] = 1
        count = count + 1
    sing_output_positives.to_csv(
        csv_path + '/cluster_' + str(cls) + '_singleton_positive_markers_ranked.csv'
    )

                
    sing_output['hgrank'] = sing_output.reset_index().index + 1
    sing_output.sort_values(by='Log2FoldChangeAbs', ascending=False, inplace=True)
    sing_output['fcrank'] = sing_output.reset_index().index + 1
    sing_output['finrank'] = sing_output[['hgrank', 'fcrank']].mean(axis=1)
    sing_output.sort_values(by='finrank',ascending=True,inplace=True)
    sing_output['rank'] = sing_output.reset_index().index + 1
    sing_output.drop('finrank',axis=1, inplace=True)
    count = 1
    for index,row in sing_output.iterrows():
        if count == 100:
            break
        sing_output.at[index,'Plot'] = 1
        count = count + 1
        
    sing_output.to_csv(
        csv_path + '/cluster_' + str(cls) + '_singleton_all_ranked.csv'
    )
    sing_stripped = sing_output[
        ['mHG_stat', 'TP', 'TN']
    ].reset_index().rename(index=str, columns={'gene_1': 'gene_1'})
    

    ranked_print = ranked_pair.head(Trim)
    ranked_print.to_csv(
        csv_path + '/cluster_' + str(cls) + '_pair_final_ranking.csv'
    )
    #Add trips data pages
    #does not currently do new rank scheme
    if K == 3:
        trips_output = trips
        #  .sort_values(by='HG_stat', ascending=True)
          #print(trips_output)
        trips_output['rank'] = trips_output.reset_index().index + 1
        trips_print = trips_output.head(Trim)
        trips_print.to_csv(
            csv_path + '/cluster_' + str(cls) + '_trips.csv'
            )
    else:
        trips_output = int(1)
    if K >= 4:
        quads_final = quads_fin\
          .sort_values(by='HG_stat', ascending=True)
        quads_final['rank'] = quads_final.reset_index().index + 1
        quads_print = quads_final.head(Trim)
        quads_print.to_csv(
            csv_path + '/cluster_' + str(cls) + '_quads.csv'
            )
    else:
        quads_final = int(1)
    print('Drawing plots...')
    #plt.bar(list(histogram.keys()), histogram.values(), color='b')
    #plt.savefig(vis_path + '/cluster_' + str(cls) + '_pair_histogram')

    #if cls == fincls:
    # cls = 0
    vis.make_plots(
        pair=ranked_pair,
        sing=sing_output,
        sing_tp_tn=sing_tp_tn,
        xlmhg=xlmhg,
        trips=trips_output,
        quads_fin=quads_final,
        tsne=tsne,
        discrete_exp=discrete_exp_full,
        marker_exp=marker_exp,
        plot_pages=plot_pages,
        combined_path=vis_path + '/cluster_' + str(cls) + '_pairs_as_singletons',
        sing_combined_path=vis_path + '/cluster_' +
        str(cls) + '_singleton',
        discrete_path=vis_path + '/cluster_' + str(cls) + '_discrete_pairs',
        tptn_path=vis_path + 'cluster_' + str(cls) + 'pair_TP_TN',
        trips_path=vis_path + 'cluster_' + str(cls) + '_discrete_trios',
        quads_path=vis_path + 'cluster_' + str(cls) + '_discrete_quads',
        sing_tptn_path=vis_path + 'cluster_' + str(cls) + '_singleton_TP_TN',
        count_data = count_data
        )
    end_cls_time=time.time()
    print(str(end_cls_time - start_cls_time) + ' seconds')
    #time.sleep(10000)


def main():
    """Hypergeometric marker detection. Finds markers identifying a cluster.

    Reads in data from single-cell RNA sequencing. Data is in the form of 3
    CSVs: gene expression data by gene by cell, 2-D tSNE data by cell, and the
    clusters of interest by cell. Creates a list of genes and a list of gene
    pairs (including complements), ranked by hypergeometric and t-test
    significance. The highest ranked marker genes generally best identify the
    cluster of interest. Saves these lists to CSV and creates gene expression
    visualizations.
    """
    # TODO: more precise description

    #ten_x = input("10X data? [y/n]")
    
    start_dt = datetime.datetime.now()
    start_time = time.time()
    print("Started on " + str(start_dt.isoformat()))
    
    args = init_parser(argparse.ArgumentParser(
        description=("Hypergeometric marker detection. Finds markers identifying a cluster. Documentation available at https://hgmd.readthedocs.io/en/latest/index.html")
    )).parse_args()
    
    output_path = args.output_path
    C = args.C
    K = args.K
    Abbrev = args.Abbrev
    Down = args.Down
    X = args.X
    L = args.L
    marker_file = args.marker
    tsne_file = args.tsne
    cluster_file = args.cluster
    gene_file = args.g
    Trim = args.Trim
    count_data = args.Count
    tenx = args.tenx
    online = args.online
    
    plot_pages = 30  # number of genes to plot (starting with highest ranked)

    # TODO: gene pairs with expression ratio within the cluster of interest
    # under [min_exp_ratio] were ignored in hypergeometric testing. This
    # functionality is currently unimplemented.
    # min_exp_ratio = 0.4

    csv_path = output_path + 'data/'
    vis_path = output_path + 'vis/'
    pickle_path = output_path + '_pickles/'
    try:
        os.makedirs(csv_path)
    except:
        os.system('rm -r ' + csv_path)
        os.makedirs(csv_path)

    try:
        os.makedirs(vis_path)
    except:
        os.system('rm -r ' + vis_path)
        os.makedirs(vis_path)

    try:
        os.makedirs(pickle_path)
    except:
        os.system('rm -r ' + pickle_path)
        os.makedirs(pickle_path)

    if Trim is not None:
        Trim = int(Trim)
    else:
        Trim = int(2000)
    if C is not None:
        C = abs(int(C))
    else:
        C = 1
    if X is not None:
        try:
            X = float(X)
        except:
            raise Exception('X param must be a number between 0 and 1')
        if X > 1:
            X = int(1)
        elif X <= 0:
            X = int(0)
        else:
            X = float(X)
        print("Set X to " + str(X) + ".")
    if L is not None:
        L = int(L)
        print("Set L to " + str(L) + ".")
    if K is not None:
        K = int(K)
    else:
        K = 2
    if K > 4:
        K = 4
        print('Only supports up to 4-gene combinations currently, setting K to 4')
    if count_data is not None:
        if count_data == str(True):
            count_data = 1
            print('Count Data')
        elif count_data == 'yes':
            count_data = 1
            print('Count Data')
        else:
            count_data = int(0)
    else:
        count_data = int(0)
    if tenx is not None:
        if tenx == str(True):
            tenx = int(1)
        elif tenx == 'yes':
            tenx = int(1)
        else:
            tenx = int(0)
    else:
        tenx = int(0)
    if online is not None:
        if online == str(True):
            online = int(1)
        elif online == 'yes':
            online = int(1)
        else:
            online = int(0)
    else:
        online = int(0)
    print("Reading data...")
    if gene_file is None:
        (cls_ser, tsne, no_complement_marker_exp, gene_path) = read_data(
            cls_path=cluster_file,
            tsne_path=tsne_file,
            marker_path=marker_file,
            gene_path=None,
            D=Down,
            tenx=tenx,
            online=online
        )
    else:
        (cls_ser, tsne, no_complement_marker_exp, gene_path) = read_data(
            cls_path=cluster_file,
            tsne_path=tsne_file,
            marker_path=marker_file,
            gene_path=gene_file,
            D=Down,
            tenx=tenx,
            online=online
        )
    print("Generating complement data...")
    marker_exp = hgmd.add_complements(no_complement_marker_exp)
    #throw out vals that show up in expression matrix but not in cluster assignments
    for ind,row in marker_exp.iterrows():
        if ind in cls_ser.index.values.tolist():
            continue
        else:
            marker_exp.drop(ind, inplace=True)
        #print(marker_exp.index.values.tolist().count(str(ind)))
        #print(marker_exp[index])
    #throw out gene rows that are duplicates and print out a message to user
    
            
    '''
    #throw out cls_ser vals not in marker_exp
    for index in cls_ser.index.values.tolist():
        if index in marker_exp.columns:
            continue
        else:
            cls_ser.drop(index,inplace=True)
    '''
    marker_exp.sort_values(by='cell',inplace=True)
    cls_ser.sort_index(inplace=True)

    # Process clusters sequentially
    clusters = cls_ser.unique()
    clusters.sort()
    cluster_overall=clusters.copy()

    #Only takes a certain number of clusters (cuts out smallest ones)
    if online == 1:
        max_clus_size = 15
        if len(clusters) <= max_clus_size:
            pass
        else:
            cls_helper = list(clusters.copy())
            cls_size_count = {}
            for item in cls_ser:
                if item in cls_size_count:
                    cls_size_count[item] = cls_size_count[item] + 1
                else:
                    cls_size_count[item] = 1
            for counted in cls_size_count:
                cls_size_count[counted] = cls_size_count[counted] / len(cls_ser)
            while len(cls_helper) > max_clus_size:
                lowest = 1
                place = 0
                for key in cls_size_count:
                    if cls_size_count[key] < lowest:
                        place = key
                        lowest = cls_size_count[key]
                cls_helper.remove(place)
                del cls_size_count[place]
            clusters = np.array(cls_helper)

            
    #Below could probably be optimized a little (new_clust not necessary),
    #cores is number of simultaneous threads you want to run, can be set at will
    cores = C
    cluster_number = len(clusters)
    # if core number is bigger than number of clusters, set it equal to number of clusters
    if cores > len(clusters):
        cores = len(clusters)
    #below loops allow for splitting the job based on core choice
    group_num  = math.ceil((len(clusters) / cores ))
    for element in range(group_num):
        new_clusters = clusters[:cores]
        print(new_clusters)
        jobs = []
        #this loop spawns the workers and runs the code for each assigned.
        #workers assigned based on the new_clusters list which is the old clusters
        #split up based on core number e.g.
        #clusters = [1 2 3 4 5 6] & cores = 4 --> new_clusters = [1 2 3 4], new_clusters = [5 6]
        for cls in new_clusters:
            p = multiprocessing.Process(target=process,
                args=(cls,X,L,plot_pages,cls_ser,tsne,marker_exp,gene_file,csv_path,vis_path,pickle_path,cluster_number,K,Abbrev,cluster_overall,Trim,count_data))
            jobs.append(p)
            p.start()
        p.join()
        new_clusters = []
        clusters = clusters[cores:len(clusters)]

    end_time = time.time()

    

    # Add text file to keep track of everything
    end_dt = datetime.datetime.now()
    print("Ended on " + end_dt.isoformat())
    metadata = open(output_path + 'metadata.txt', 'w')
    metadata.write("Started: " + start_dt.isoformat())
    metadata.write("\nEnded: " + end_dt.isoformat())
    metadata.write("\nElapsed: " + str(end_dt - start_dt))
    #metadata.write("\nGenerated by COMET version " + conf.version)


    print('Took ' + str(end_time-start_time) + ' seconds')
    print('Which is ' + str( (end_time-start_time)/60 ) + ' minutes')

if __name__ == '__main__':
    main()
