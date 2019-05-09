"""
Set of modularized components of COMET's HGMD testing.

For marker expression, float comparisions are fuzzy to 1e-3.  Marker expression
must therefore be normalized to a point where a difference of 0.001 is
insignificant.  I.e. 15.001 and 15.000 are treated as equivalent expression
values.
"""

import re

import pandas as pd
import numpy as np
import xlmhg as hg
import scipy.stats as ss
import time
import math
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from . import qvalue

# TODO: idea, replace 'gene_1', 'gene_2', and 'gene' columns with
# indices/multi-indices

# Used for comparision of marker expression values.
FLOAT_PRECISION = 0.001


def add_complements(marker_exp):
    """Adds columns representing gene complement to a gene expression matrix.

    Gene complements are represented simplistically: gene expression values for
    a given gene X are multiplied by -1 and become a new column, labeled X_negation.
    "High" expression values of X become "low" values of X_negation, and vice versa,
    where discrete expression corresponds to a "high" value, and discrete
    non-expression to a "low" value.

    marker_exp should have cell row labels, gene column labels, gene expression
    float values.

    :param marker_exp: gene expression DataFrame whose rows are cell
        identifiers, columns are gene identifiers, and values are float values
        representing gene expression.

    :returns: A DataFrame of same format as marker_exp, but with a new column
              added for each existing column label, representing the column
              label gene's complement.

    :rtype: pandas.DataFrame
    """
    for gene in marker_exp.columns:
        #marker_exp[gene + '_negation'] = 1/(1+marker_exp[gene])
        marker_exp[gene + '_negation'] = -marker_exp[gene]
    '''
    for gene in marker_exp.columns:
        marker_exp[gene + '_c'] = max(marker_exp[gene]) - marker_exp[gene]
    '''
    return marker_exp


def batch_xlmhg(marker_exp, c_list, coi, X=None, L=None):
    """Applies XL-mHG test to a gene expression matrix, gene by gene.

    Outputs a 3-column DataFrame representing statistical results of XL-mHG.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.
    :param X: An integer to be used as argument to the XL-mHG test.
    :param L: An integer to be used as argument to the XL-mHG test.

    :returns: A matrix with arbitrary row indices, whose columns are the gene
              name, stat, cutoff, and pval outputs of the XL-mHG test; of
              float, int, and float type respectively.  Their names are 'gene',
              'mHG_stat', 'mHG_cutoff', and 'mHG_pval'.

    :rtype: pandas.DataFrame
    """
    # * 1 converts to integer
    mem_list = (c_list == coi) * 1
    #count_n = 0
    #count_2n = 0
    #Count the number of cells in the cluster, store into count_n
    count_n = np.sum(mem_list)
    '''
    for cell in mem_list:
        if cell == 1:
            count_n = count_n + 1
        else:
            continue
    '''
    #Twice the number of cells in the cluster
    #count_2n = count_n * 2
    #Set X and L params
    if X is None:
        X = np.int(.15*count_n)
    else:
        X = np.int(X*count_n)
    if L is None:
        if 2*count_n >= marker_exp.shape[0]:
            L = np.int(marker_exp.shape[0])
        else:
            L = np.int(2*count_n)
       #L = marker_exp.shape[0]
    print('X = ' + str(X))
    print('L = ' + str(L))
    print('Cluster size ' + str(count_n))
    xlmhg = marker_exp.apply(
        lambda col:
        hg.xlmhg_test(
            mem_list.reindex(
                col.sort_values(ascending=False).index
            ).values,
            X=X,
            L=L
        )
    )
    xlmhg_1 = marker_exp.apply(
        lambda col:
        hg.get_xlmhg_test_result(N=len(mem_list),
                                     indices=np.array(np.where(np.array(mem_list.reindex(col.sort_values(ascending=False).index).values==1))[0],dtype='int64').astype('uint16')
                                     ,X=X,L=L,pval_thresh=1e-12,escore_pval_thresh=1e-1,tol=1e-12).escore
        )
    output = pd.DataFrame()
    output['gene_1'] = xlmhg.index
    output[['mHG_stat', 'mHG_cutoff', 'mHG_pval']] = pd.DataFrame(
        xlmhg.values.tolist(),
        columns=['mHG_stat', 'mHG_cutoff', 'mHG_pval']
    )
    '''
    print(xlmhg_1)
    output['escore'] = np.array(xlmhg_1,dtype=float)
    output.fillna(0,inplace=True)
    print(output.sort_values(by='escore', ascending=False))
    time.sleep(1000)
    '''
    return output



def batch_q(xlmhg):
    #Takes in the xlmhg dataframe which includes gene names and pvalue assignments
    #Returns a list of q-value, corresponding to each gene
    try:
        q_vals = qvalue.estimate(np.array(xlmhg['mHG_pval']))
    except:
        print('error in qval calculation')
    q_val_out = pd.DataFrame(data=xlmhg['gene_1'].tolist(),columns=['gene_1'])
    q_val_out['q_value'] = q_vals

    return q_val_out

def pairs_q(pair):
    try:
        q_vals = qvalue.estimate(np.array(pair['HG_pval']))
        q_val_out_pair = pd.DataFrame(data=pair['gene_1'].tolist(),columns=['gene_1'])
        q_val_out_pair['gene_2'] = pair['gene_2']
        q_val_out_pair['q_value'] = q_vals
        return q_val_out_pair
    except:
        print('error in qval calculation')
        q_val_out_pair = pd.DataFrame(data=pair['gene_1'].tolist(),columns=['gene_1'])
        q_val_out_pair['gene_2'] = pair['gene_2']
        q_val_out_pair['q_value'] = 0
        return q_val_out_pair

def batch_stats_extended(marker_exp, c_list, coi):
    """Applies t test , wilcoxon test, and likelihood ratio test (Based on logistic regression)
    to a gene expression matrix, gene by gene. Also gives simple up versus down regulation test (difference between means).

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix with arbitary row indices whose columns are the gene, t
              statistic, then t p-value; the last two being of float type.
              Their names are 'gene', 't_stat' , 't_pval' , w_stat, w_pval , LRT_pval, up/down regulated

    :rtype: pandas.DataFrame
    """
    def LRT_LogReg(df):
        # Define model matrix and response
        X = np.matrix(df.drop('cluster', axis=1))
        y = df['cluster']
        # Train logistic regression with full model
        logreg1 = LogisticRegression(solver='lbfgs').fit(X,y)
        ll1 = -log_loss(y,logreg1.predict_proba(X),normalize=False)
        # Train logistic regression with null model (only intercept)
        logreg0 = LogisticRegression(solver='lbfgs').fit([[0]]*len(X) ,y)
        ll0 = -log_loss(y,logreg0.predict_proba(X),normalize=False)
        # Likelihood ratio test
        stat = 2*(ll1-ll0)
        pval = ss.chi2.sf(stat, 1)
        return(pval)
    LRT_pvals = []
    up_v_down_vals = []
    for column in marker_exp:
        log_reg_in = pd.DataFrame(data=[marker_exp[column]])
        log_reg_in = np.transpose(log_reg_in)
        c_list_2 = np.array(c_list)
        c_list_2 = np.array(c_list_2 == coi,dtype=int)
        c_list_2 = np.transpose(c_list_2)
        log_reg_in['cluster'] = c_list_2
        in_cls = marker_exp[column][c_list == coi].values
        out_cls = marker_exp[column][c_list != coi].values
        out_cls_mean = np.sum(out_cls) / len(out_cls)
        in_cls_mean = np.sum(in_cls) / len(in_cls)
        test = in_cls_mean - out_cls_mean
        if test <= 0:
            up_v_down_vals.append('down')
        else:
            up_v_down_vals.append('up')
        
        LRT_pval = LRT_LogReg(log_reg_in)
        LRT_pvals.append(LRT_pval)
        
    t = marker_exp.apply(
        lambda col:
        ss.ttest_ind(
            col[c_list == coi],
            col[c_list != coi],
            equal_var=False
        )
    )
    ws = marker_exp.apply(
        lambda col:
        ss.ranksums(
            col[c_list == coi],
            col[c_list != coi]
        )
    )
    output = pd.DataFrame()
    output['gene_1'] = t.index
    #output['gene_1'] = ws.index
    output[['t_stat', 't_pval']] = pd.DataFrame(
        t.values.tolist(),
        columns=['t_stat', 't_pval']
    )
    output[['w_stat', 'w_pval']] = pd.DataFrame(
        ws.values.tolist(),
        columns=['w_stat', 'w_pval']
    )

    output['up_down'] = up_v_down_vals
    output['LRT_pval'] = LRT_pvals
    return output


def batch_stats(marker_exp, c_list, coi):
    """Applies t test & wilcoxon rank sum test to a gene expression matrix, gene by gene.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix with arbitary row indices whose columns are the gene, t
              statistic, then t p-value; the last two being of float type.
              Their names are 'gene', 't_stat' , 't_pval' , 'w_stat' , 'w_pval'

    :rtype: pandas.DataFrame
    """
        
    t = marker_exp.apply(
        lambda col:
        ss.ttest_ind(
            col[c_list == coi],
            col[c_list != coi],
            equal_var=False
        )
    )
    ws = marker_exp.apply(
        lambda col:
        ss.ranksums(
            col[c_list == coi],
            col[c_list != coi]
        )
    )
    output = pd.DataFrame()
    output['gene_1'] = t.index
    #output['gene_1'] = ws.index
    output[['t_stat', 't_pval']] = pd.DataFrame(
        t.values.tolist(),
        columns=['t_stat', 't_pval']
    )
    output[['w_stat', 'w_pval']] = pd.DataFrame(
        ws.values.tolist(),
        columns=['w_stat', 'w_pval']
    )

    return output


def batch_fold_change(marker_exp, c_list, coi):
    """Applies log2 fold change to a gene expression matrix, gene by gene.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.
    :rtype: pandas.DataFrame
    """
    
    def fold_change(col,c_list,coi):
        mean0 = np.mean(col[c_list == coi]) + .000001
        mean1 = np.mean(col[c_list != coi]) + .000001
        if mean0 == 0:
            val = math.nan
            return val
        if mean1 == 0:
            val = math.nan
            return val
        val = (mean0/mean1)
        if val < 0:
            return abs(val)
        return val

    
    fc = marker_exp.apply(
        lambda col:
        (math.log(fold_change(col,c_list,coi),2))
        )
    fca = marker_exp.apply(
        lambda col:
        abs(math.log(fold_change(col,c_list,coi),2))
        )
    output = pd.DataFrame()
    output['gene_1'] = fc.index
    output[['Log2FoldChange']] = pd.DataFrame(
        fc.values.tolist(),
        columns=['Log2FoldChange']
    )
    output['gene_1'] = fca.index
    output[['Log2FoldChangeAbs']] = pd.DataFrame(
        fca.values.tolist(),
        columns=['Log2FoldChangeAbs']
    )
    return output


def mhg_cutoff_value(marker_exp, cutoff_ind):
    """Finds discrete expression cutoff value, from given cutoff index.

    The XL-mHG test outputs the index of the cutoff of highest significance
    between a sample and population.  This functions finds the expression value
    which corresponds to this index.  Cells above this value we define as
    expressing, and cells below this value we define as non-expressing.  We
    therefore choose this value to be between the expression at the index, and
    the expression of the "next-highest" cell.  I.e. for expression [3.0 3.0
    1.5 1.0 1.0] and index 4, we should choose a cutoff between 1 and 1.5. This
    implementation will add epsilon to the lower bound (i.e. the value of
    FLOAT_PRECISION).  In our example, the output will be 1.0 +
    FLOAT_PRECISION.  For FLOAT_PRECISION = 0.001, this is 1.001.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_ind: A DataFrame whose 'gene' column are gene identifiers,
        and whose 'mHG_cutoff' column are cutoff indices

    :returns: A DataFrame whose 'gene' column are gene identifiers, and whose
              'cutoff_val' column are cutoff values corresponding to input
              cutoff indices.

    :rtype: pandas.DataFrame
    """

    def find_val(row):
        gene = row['gene_1']
        val = marker_exp[gene].sort_values(
            ascending=False).iloc[row['mHG_cutoff']]
        if re.compile(".*_negation$").match(gene):
            return val + FLOAT_PRECISION
        else:
            return val + FLOAT_PRECISION

    cutoff_ind.index = cutoff_ind['gene_1']
    cutoff_val = cutoff_ind.apply(
        find_val, axis='columns'
    ).rename('cutoff_val')
    output = cutoff_val.to_frame().reset_index()
    return output


def mhg_slide(marker_exp, cutoff_val):
    """Slides cutoff indices in XL-mHG output out of uniform expression groups.

    The XL-mHG test may place a cutoff index that "cuts" across a group of
    uniform expression inside the sorted expression list.  I.e. for a
    population of cells of which many have zero expression, the XL-mHG test may
    demand that we sample some of the zero-expression cells and not others.
    This is impossible because the cells are effectively identical.  This
    function therefore moves the XL-mHG cutoff index so that it falls on a
    measurable gene expression boundary.

    Example: for a sorted gene expression list [5, 4, 1, 0, 0, 0] and XL-mHG
    cutoff index 4, this function will "slide" the index to 3; marking the
    boundary between zero expression and expression value 1.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_val: A DataFrame whose 'gene' column are gene identifiers,
        and whose 'cutoff_val' column are cutoff values corresponding to input
        cutoff indices.

    :returns: A DataFrame with 'gene', 'mHG_cutoff', and 'cutoff_val' columns,
              slid.

    :rtype: pandas.DataFrame
    """
    '''
    def searcher(row):
        if re.compile(".*_negation$").match(row['gene_1']):
            val_ =  np.searchsorted(
                -marker_exp[row['gene_1']].sort_values(ascending=False).values,
                -row['cutoff_val'], side='left')
        else:
            val_ = np.searchsorted(
                -marker_exp[row['gene_1']].sort_values(ascending=False).values,
                -row['cutoff_val'], side='left')
        
        return val_
        
    cutoff_val.index = cutoff_val['gene_1']
    cutoff_ind = cutoff_val.apply(lambda row: searcher(row),axis='columns')
    '''
    cutoff_val.index = cutoff_val['gene_1']
    cutoff_ind = cutoff_val.apply(
        lambda row:
        np.searchsorted(
            -marker_exp[row['gene_1']].sort_values(ascending=False).values,
            -row['cutoff_val'], side='left'
        ),
        axis='columns'
    )
    output = cutoff_val
    output['mHG_cutoff'] = cutoff_ind
    # Reorder and remove redundant row index
    output = output.reset_index(
        drop=True)[['gene_1', 'mHG_cutoff', 'cutoff_val']]
    return output


def discrete_exp(marker_exp, cutoff_val,abbrev,xlmhg):
    """Converts a continuous gene expression matrix to discrete.

    As a note: cutoff values correspond to the "top" of non-expression.  Only
    cells expressing at values greater than the cutoff are marked as
    "expressing"; cells expressing at the cutoff exactly are not.

    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param cutoff_val: A Series whose rows are gene identifiers, and values are
        cutoff values.

    :returns: A gene expression matrix identical to marker_exp, but with
              boolean rather than float expression values.

    :rtype: pandas.DataFrame
    """
    output = pd.DataFrame()
    for gene in marker_exp.columns:
        output[gene] = (marker_exp[gene] > cutoff_val[gene]) * 1
    ab = '2'
    abb = '3'
    if ab in abbrev or abb in abbrev:
        #Ensures the matrix is in the proper order to get the top however many genes for the heuristic approach
        #create 'new_order', a list with the new ordering. trips_list goes first, rest after in any order
        rest = []
        new_order = xlmhg['gene_1'].tolist()
        output = output[new_order]
    return output


def tp_tn(discrete_exp, c_list, coi, cluster_overall):
    """Finds simple true positive/true negative values for the cluster of
    interest.

    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.

    :returns: A matrix with arbitary row indices, and has 3 columns: one for
              gene name, then 2 containing the true positive and true negative
              values respectively.  Their names are 'gene', 'TP', and 'TN'.

    :rtype: pandas.DataFrame
    """

    #does rest of clusters
    discrete_exp.fillna(0, inplace=True)
    sing_cluster_exp_matrices = {}
    for clstrs in cluster_overall:
        if clstrs == coi:
            continue
        mem_list = (c_list == (clstrs)) * 1
        tp_tn =discrete_exp.apply(
            lambda col: (
                np.dot(mem_list, col.values) / np.sum(mem_list),
                np.dot(1 - mem_list, 1 - col.values) / np.sum(1 - mem_list),
            )
        )
        sing_cluster_exp_matrices[clstrs] = pd.DataFrame()
        sing_cluster_exp_matrices[clstrs]['gene_1'] = tp_tn.index
        sing_cluster_exp_matrices[clstrs][['TP', 'TN']] = pd.DataFrame(
            tp_tn.values.tolist(),
            columns=['TP', 'TN']
        )
        sing_cluster_exp_matrices[clstrs].set_index('gene_1',inplace=True)
        

    #does our cluster of interest
    # * 1 converts to integer
    mem_list = (c_list == coi) * 1
    
    tp_tn = discrete_exp.apply(
        lambda col: (
            np.dot(mem_list, col.values) / np.sum(mem_list),
            np.dot(1 - mem_list, 1 - col.values) / np.sum(1 - mem_list)
        )
    )
    output = pd.DataFrame()
    output['gene_1'] = tp_tn.index
    output[['TP', 'TN']] = pd.DataFrame(
        tp_tn.values.tolist(),
        columns=['TP', 'TN']
    )
    #outputs a DF for COI and a dict of DF's for rest
    return output, sing_cluster_exp_matrices


def pair_product(discrete_exp, c_list, coi, cluster_number,cluster_overall):
    """Finds paired expression counts.  Returns in matrix form.

    (Number of cells in which two genes are coexpressed)

    The product of the transpose of the discrete_exp DataFrame is a matrix
    whose rows and columns correspond to individual genes.  Each value is the
    number of cells which express both genes (i.e. the dot product of two lists
    of 1s and 0s encoding expression/nonexpression for their respective genes
    in the population).  The product therefore encodes joint expression counts
    for any possible gene pair (including a single gene paired with itself).

    This function produces two matrices: one considering only cells inside the
    cluster of interest, and one considering all cells in the population.

    This function also produces a list mapping integer indices to gene names,
    and the population cell count.

    Additionally, only the upper triangular part of the output matrices is
    unique.  This function therefore also returns the upper triangular indices
    for use by other functions; this is a lazy workaround for the issue that
    comes with using columns 'gene_1' and 'gene_2' to store gene pairs; the
    gene pair (A, B) is therefore treated differently than (B, A).  Specifying
    the upper triangular part prevents (B, A) from existing.

    TODO: fix this redundancy using multi-indices

    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param c_list: A Series whose indices are cell identifiers, and whose
        values are the cluster which that cell is part of.
    :param coi: The cluster of interest.
    :param clusters_number: dict with clusters and their sizes

    :returns: (gene mapping list, cluster count, total count, cluster paired
              expression count matrix, population paired expression count
              matrix, upper triangular matrix index)

    :rtype: (pandas.Index, int, int, numpy.ndarray, numpy.ndarray,
            numpy.ndarray)
    """

    '''
    ###############
    #Heuristic Section
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
    '''
    
    gene_map = discrete_exp.columns
    in_cls_matrix = discrete_exp[c_list == coi].values
    total_matrix = discrete_exp.values
    #get exp matrices for each non-interest cluster for new rank scheme
    cluster_exp_matrices = {}
    cls_counts = {}
    for clstrs in cluster_overall:
        if clstrs == coi:
            pass
        else:
            cluster_exp_matrices[clstrs] = discrete_exp[c_list == (clstrs)].values
            cls_counts[clstrs] = np.size(cluster_exp_matrices[clstrs],0)
            cluster_exp_matrices[clstrs]=np.matmul(np.transpose(cluster_exp_matrices[clstrs]),cluster_exp_matrices[clstrs] )

    in_cls_count = np.size(in_cls_matrix, 0)
    pop_count = np.size(total_matrix, 0)
    in_cls_product = np.matmul(np.transpose(in_cls_matrix), in_cls_matrix)
    total_product = np.matmul(np.transpose(total_matrix), total_matrix)
    upper_tri_indices = np.triu_indices(gene_map.size,1)
    
    return (
        gene_map,
        in_cls_count, pop_count,
        in_cls_product, total_product,
        upper_tri_indices, cluster_exp_matrices, cls_counts
    )

  
    

def combination_product(discrete_exp,c_list,coi,abbrev,heur_limit):
    '''
    Makes the count matrix for combinations of K genes where K > 2 :)

    X -> in_cls_matrix

    First constructs the extended gene by cell matrix then multiplies by the transpose
    -> (X^T & X^T) * X
    & represents AND bitwise operation to construct the gene combos by cells matrix
    (need to perform AND on each row in first X^T by each row in second X^T)
    This scales to K genes , just need to do more ANDs with more X^Ts

    Once the first step is complete, simply multiply using np.matmul to get the gene by gene
    matrix. Need to also make a new gene_map

    More in-depth Description:
    the trips matrix is a gene/gene pair by gene matrix. TO construct, we first perform an AND 
    operation on each given row with a different row from the same matrix(X&X), we use two 
    identical matrices and loop through most combinations. The code SKIPS any entry that would be
    the same pair (e.g. AA) and also doesn't do inverses of already counted pairs
    (e.g. do AB, dont do BA). this gives us a somewhat unique looking rectangular matrix that isnt
    quit N^2 long. From here, we  multiply by the original matrix to get our full gene expression
    count matrix. To be absolutely clear, the matrices we use to construct the ~N^2 gene pair 
    matrix are actually the original discrete expression matrix transposed so that the genes
    and cells end up in the right place. We do this for our cluster and the total population,
    just as in the pair case. After this, we need to make a smarter gene mapping since the
    terms are not as easy to predict as in a square matrix. To do this, a pattern was determined
    (specifically, the pairs follow the N'th triangular number scheme), so from there we construct
    the full mapping for genes 1 & 2 to be used later in the data table. 'Trips_indices' gives us
    the indices for the entire rectangular matrix, something to feed into the vectorized HG & TPTN
    when we get there. Gene 3 mapping follows these indices so there is no separate mapping.

    '''

    def trips_matrix_gen_heuristic(matrix,heur_limit):
        #transpose gives us genes as rows
        Cth = time.time()
        trips_matrix_list = []
        row1count = 1
        for row1 in np.transpose(matrix)[:heur_limit]:
            row2count = 1
            for row2 in np.transpose(matrix):
                if row2count<=row1count:
                    row2count = row2count+1
                    continue
                trips_matrix_list.append(row1&row2)
                row2count=row2count+1
            row1count = row1count + 1
        #trips_matrix is in cluster combo gene by cell
        Ath = time.time()
        trips_matrix = np.asarray(trips_matrix_list)
        trips_product = np.matmul(trips_matrix,np.transpose(np.transpose(matrix)[:heur_limit]))
        Bth = time.time()
        print(str(Bth-Ath) + ' mult seconds')
        return trips_product
    

    def trips_matrix_gen(matrix):
        Cth = time.time()
        trips_matrix_list = []
        row1count = 1
        for row1 in np.transpose(matrix):
            row2count = 1
            for row2 in np.transpose(matrix):
                if row2count<=row1count:
                    row2count = row2count+1
                    continue
                trips_matrix_list.append(row1&row2)
                row2count=row2count+1
            row1count = row1count + 1

        #trips_matrix is in cluster combo gene by cell
        Ath = time.time()
        trips_matrix = np.asarray(trips_matrix_list)
        trips_product = np.matmul(trips_matrix,matrix)
        Bth = time.time()
        print(str(Bth-Ath) + ' mult seconds')
        return trips_product

    
    
    start_time = time.time()
    in_cls_matrix = discrete_exp[c_list == coi].values
    total_matrix = discrete_exp.values
    gene_count = len(discrete_exp.columns)
    first = time.time()
    if '3' in abbrev:
        total_pairs =  0
        for i in range(1,heur_limit+1):
            total_pairs = total_pairs + (len(discrete_exp.columns.values) - i)
        total_pairs = total_pairs * len(discrete_exp.columns.values)
        trips_in_cls_product = trips_matrix_gen_heuristic(in_cls_matrix,heur_limit)
        trips_total_product = trips_matrix_gen_heuristic(total_matrix,heur_limit)
    else:
        trips_in_cls_product = trips_matrix_gen(in_cls_matrix)
        trips_total_product = trips_matrix_gen(total_matrix)

    #make a row-wise gene_map scheme
    gene_map = discrete_exp.columns.values
    gene_1_mapped = []
    count = 0
    count_2 = 0
    for gene in gene_map:
        if '3' in abbrev and count_2 >= total_pairs+1:
            break
        for x in range(gene_count-1-count):
            for n in range(gene_count):
                gene_1_mapped.append(gene)
                count_2 = count_2 + 1
        count = count + 1
    gene_1_mapped = pd.Index(gene_1_mapped)
    if '3' in abbrev:
        gene_1_mapped = gene_1_mapped[:total_pairs]
    count_2 = 0
    gene_3_mapped = []
    val = int(len(gene_1_mapped)/gene_count)
    for x in range(val):
        if '3' in abbrev and count_2 >= total_pairs+1:
            break
        count = 0
        for gene in gene_map:
            gene_3_mapped.append(gene_map[count])
            count_2 = count_2 + 1
            count = count + 1
    gene_3_mapped = pd.Index(gene_3_mapped)
    if '3' in abbrev:
        gene_3_mapped = gene_3_mapped[:total_pairs]
    count_2 = 0
    gene_2_mapped = []
    for x in range(gene_count-1):
        if '3' in abbrev and count_2 >= total_pairs+1:
            break
        count = 0
        for gene in gene_map:
            if count == 0:
                count = count + 1
                continue
            for n in range(gene_count):
                gene_2_mapped.append(gene_map[count])
                count_2 = count_2 + 1
            count = count + 1
        gene_map = np.delete(gene_map,0)
    gene_2_mapped = pd.Index(gene_2_mapped)
    if '3' in abbrev:
        gene_2_mapped = gene_2_mapped[:total_pairs]

    #indices for computation
    row_count = int((gene_count*(gene_count-1))/2)
    
    #column coordinate
    count_2 = 0
    list_one_ = []
    for numm in range(row_count):
        if '3' in abbrev and count_2 >= total_pairs+1:
            break
        for num in range(gene_count):
            list_one_.append(num)
            count_2 = count_2 + 1
    list_one = np.array(list_one_)
    
    #row coordinate
    count_2 = 0
    list_two_ = []
    for num in range(row_count):
        if '3' in abbrev and count_2 >= total_pairs+1:
            break
        for numm in range(gene_count):
            list_two_.append(num)
            count_2 = count_2 + 1
    list_two = np.array(list_two_)
    
    fourth = time.time()

    if '3' in abbrev:
        list_one = list_one[:total_pairs]
        list_two = list_two[:total_pairs]
        
    trips_indices = ( list_two , list_one )
    
    return (
        trips_in_cls_product,trips_total_product,trips_indices,
        gene_1_mapped,gene_2_mapped,gene_3_mapped
        )


def pair_hg(gene_map, in_cls_count, pop_count, in_cls_product, total_product,
            upper_tri_indices, abbrev,heur_limit):
    """Finds hypergeometric p-value of gene pairs.

    Takes in discrete single-gene expression matrix, and finds the
    hypergeometric p-value of the sample that includes cells which express both
    of a pair of genes. 

    hypergeometric parameters:
    k = in_cls_product (number of cells in cluster with coexpression of both genes in the given pair)
    M = pop_count (total number of cells)
    n = in_cls_count (number of cells in cluster)
    N = total_product  (total number of cells with coexpression of both genes in the given pair)

    :param gene_map: An Index mapping index values to gene names.
    :param in_cls_count: The number of cells in the cluster.
    :param pop_count: The number of cells in the population.
    :param in_cls_product: The cluster paired expression count matrix.
    :param total_product: The population paired expression count matrix.
    :param upper_tri_indices: An array specifying UT indices; from numpy.utri

    :returns: A matrix with columns: the two genes of the pair, hypergeometric
              test statistics for that pair.  Their names are 'gene_1',
              'gene_2', 'HG_stat'.

    :rtype: pandas.DataFrame

    """
    
    vhg = np.vectorize(ss.hypergeom.sf, excluded=[1, 2, 4], otypes=[np.float])
    ab = '2'
    if ab in abbrev:
        count = 0
        row_count = 1
        heuristic_limit = heur_limit
        for num in range(heuristic_limit):
            count = count + (len(gene_map) - row_count)
            row_count = row_count + 1
        revised_row_indices = np.asarray(upper_tri_indices[0][:count])
        revised_col_indices = np.asarray(upper_tri_indices[1][:count])
        revised_indices = (revised_row_indices,revised_col_indices)
        # Only apply to upper triangular    
        hg_result = vhg(
            in_cls_product[revised_indices],
            pop_count,
            in_cls_count,
            total_product[revised_indices],
            loc=1
            )
        output = pd.DataFrame({
            'gene_1': gene_map[revised_indices[0]],
            'gene_2': gene_map[revised_indices[1]],
            'HG_pval': hg_result
            }, columns=['gene_1', 'gene_2', 'HG_pval'])
        return output, revised_indices

    
    else:
        hg_result = vhg(
            in_cls_product[upper_tri_indices],
            pop_count,
            in_cls_count,
            total_product[upper_tri_indices],
            loc=1
        )
        output = pd.DataFrame({
            'gene_1': gene_map[upper_tri_indices[0]],
            'gene_2': gene_map[upper_tri_indices[1]],
            'HG_pval': hg_result
        }, columns=['gene_1', 'gene_2', 'HG_pval'])
        return output, None



def ranker(pair,xlmhg,sing_tp_tn,other_sing_tp_tn,other_pair_tp_tn,cls_counts,in_cls_count,pop_count):
    """
    :param pair: table w/ gene_1, gene_2, HG_stat as columns (DataFrame (DF) )
    :param xlmhg: DF of mHG stats for each gene(for testing lead vs follow gene in pair)
    :param other_sing_tp_tn: TP/TN values for singletons in all other clusters (dict of DFs)
    :param other_pair_tp_tn: TP/TN values for pairs in all other clusters (dict of DFs)
    :param cls_counts: # of cells in a given cluster (dict of DFs)
    **
    All dicts have style: 
    key -> cluster number
    value -> data
    **

    Statistic to calculate is :
    SUM across all clusters of (TN_after - TN_before) / N
    where:
    TN_after = TN of gene combo in cluster
    TN_before = TN of initial gene from pair in cluster
    N = # of cells in the cluster
    
    THRESHOLDING:
    -We are only taking the first 100 appearances for each gene in all the pairs
    -'Lead' gene is one with smallest p val

    returns: New pair table w/ new columns and ranks.
    ranked-pair is a DEEP copy of pair, meaning value changes in it
    are not reflected in pair 
    """

    def ranked_stat(gene_1,gene_2,lead_gene,cls_counts,other_pair_tp_tn,other_sing_tp_tn,in_cls_count):
        stats=[]
        MGDstats=[]
        stats_debug = {}
        for clstrs in cls_counts:
            TN_before = other_sing_tp_tn[clstrs].at[lead_gene,'TN']
            TN_after = other_pair_tp_tn[clstrs].at[(gene_1,gene_2),'TN']
            N = cls_counts[clstrs]
            #value =  ( TN_after - TN_before ) / N
            value =  ( TN_after - TN_before )
            stats.append(value)
            MGDstats.append(((TN_after - TN_before)*N)/pop_count)
        stat1 = sum(stats)
        MGDstat = sum(MGDstats)
        #stats.remove(stat1)
        #stat2 = max(stats)
        #stat= stat1+stat2
        return stat1,MGDstat

    
    xlmhg_cop = xlmhg.copy()
    xlmhg_cop = xlmhg_cop.set_index('gene_1')
    ranked_pair = pair.copy()
    num_tests = len(ranked_pair)
    ranked_pair.sort_values(by='HG_pval',ascending=True,inplace=True)
    ranked_pair['HG_rank'] = ranked_pair.reset_index().index + 1
    
    #below not used because this does ALL pairs (too many)
    #ranked_pair['CCS'] = ranked_pair.apply(ranked_stat,axis=1,args=(cls_counts,other_pair_tp_tn,other_sing_tp_tn))
    omit_pairs = {}
    count = 1
    if len(ranked_pair.index) < 5000:
        thresh = 1000
    else:
        thresh = 1000
    loopstart = time.time()
    
    for index,row in ranked_pair.iterrows():
        if row[0] in omit_pairs:
            if omit_pairs[row[0]] > 200:
                ranked_pair.drop(index, inplace=True)
                continue
        if row[1] in omit_pairs:
            if omit_pairs[row[1]] > 200:
                ranked_pair.drop(index, inplace=True)
                continue
        
        gene_1 = row[0]
        gene_2 = row[1]
        #determine which of the two genes is 'lead'
        stat_1 = xlmhg_cop.at[gene_1,'mHG_pval']
        stat_2 = xlmhg_cop.at[gene_2,'mHG_pval']
        if stat_1 <= stat_2:
            lead_gene=gene_1
            follow_gene = gene_2
        else:
            lead_gene=gene_2
            follow_gene = gene_1


        if row[0] in omit_pairs:
            omit_pairs[row[0]] = omit_pairs[row[0]] + 1
        if row[1] in omit_pairs:
            omit_pairs[row[1]] = omit_pairs[row[1]] + 1
            
        if row[0] not in omit_pairs:
            omit_pairs[row[0]] = 1
        if row[1] not in omit_pairs:
            omit_pairs[row[1]] = 1
        
        if count == thresh:
            break
        
        ranked_pair.at[index,'CCS'],ranked_pair.at[index,'MGD'] = ranked_stat(gene_1,gene_2,lead_gene,cls_counts,other_pair_tp_tn,other_sing_tp_tn,in_cls_count)
        count = count + 1

    loopend = time.time()
    ranked_pair.sort_values(by='CCS',ascending=False,inplace=True)
    ranked_pair['CCS_rank'] = ranked_pair.reset_index().index + 1
    ranked_pair['finalrank'] = ranked_pair[['HG_rank', 'CCS_rank']].mean(axis=1)
    ranked_pair.sort_values(by='finalrank',ascending=True,inplace=True)
    ranked_pair['rank'] = ranked_pair.reset_index().index + 1
    ranked_pair.drop('finalrank',axis=1,inplace=True)
    omit_genes_lead = {}
    omit_genes_follow = {}
    count = 0
    if len(ranked_pair.index) < 5000:
        plot_num = 50
    else:
        plot_num = 100
    for index,row in ranked_pair.iterrows():
        lead_on = 1
        follow_on = 1
        if count == plot_num:
            break
        #If a lead gene has appeared more than 10 times, dont plot anymore
        if row[0] in omit_genes_lead:
            if omit_genes_lead[row[0]] >= 10:
                lead_on = 0
                #ranked_pair.at[index,'Plot'] = 0
            else:
                # test if the non-adjusted pvalue is less than .05
                if row[2]  >= .05:
                    lead_on = 0
                    #ranked_pair.at[index,'Plot'] = 0
                else:
                    #If true positive rate is less than 15%, dont plot
                    if row[3] <= .15:
                        lead_on = 0
                        #ranked_pair.at[index,'Plot'] = 0
                    else:
                        pass
                        #ranked_pair.at[index,'Plot'] = 1
                omit_genes_lead[row[0]] = omit_genes_lead[row[0]]+1
        else:
            omit_genes_lead[row[0]] = 1
            if row[2] >= .05:
                lead_on = 0
                #ranked_pair.at[index,'Plot'] = 0
            else:
                pass
                #ranked_pair.at[index,'Plot'] = 1
        ##LEAD GENE END##
        
        # Follow Gene Start #
        if row[1] in omit_genes_follow:
            if omit_genes_follow[row[1]] >= 20:
                follow_on = 0
                #ranked_pair.at[index,'Plot'] = 0
            else:
                # test if the non-adjusted pvalue is less than .05
                if row[2]  >= .05:
                    follow_on = 0
                    #ranked_pair.at[index,'Plot'] = 0
                else:
                    #If true positive rate is less than 15%, dont plot
                    if row[3] <= .15:
                        follow_on = 0
                        #ranked_pair.at[index,'Plot'] = 0
                    else:
                        pass
                        #ranked_pair.at[index,'Plot'] = 1
                omit_genes_follow[row[1]] = omit_genes_follow[row[1]]+1
        else:
            omit_genes_follow[row[1]] = 1
            if row[2] >= .05:
                lead_on = 0
                #ranked_pair.at[index,'Plot'] = 0
            else:
                pass
                #ranked_pair.at[index,'Plot'] = 1
        
        ##FOLLOW GENE END##
        if follow_on == 1 and lead_on == 1:
            ranked_pair.at[index,'Plot'] = 1
            count = count + 1
        else:
            ranked_pair.at[index,'Plot'] = 0
        
    pop_list = []
    for key, value in omit_pairs.items():
        if value <= 20:
            pop_list.append(key)

    for item in pop_list:
        omit_pairs.pop(item)
    return ranked_pair,omit_pairs


def pair_tp_tn(gene_map, in_cls_count, pop_count, in_cls_product,
               total_product, upper_tri_indices, abbrev, revised_indices):
    """Finds simple true positive/true negative values for the cluster of
    interest, for all possible pairs of genes.

    :param gene_map: An Index mapping index values to gene names.
    :param in_cls_count: The number of cells in the cluster.
    :param pop_count: The number of cells in the population.
    :param in_cls_product: The cluster paired expression count matrix.
    :param total_product: The population paired expression count matrix.
    :param upper_tri_indices: An array specifying UT indices; from numpy.utri
    :param cluster_exp_matrices: dict containing expression matrices of 
         all clusters except the cluster of interest

    :returns: A matrix with arbitary row indices and 4 columns: containing the
              two genes of the pair, then true positive and true negative
              values respectively.  Their names are 'gene_1', 'gene_2', 'TP',
              and 'TN'.

    :rtype: pandas.DataFrame

    """


    def tp(taken_in_cls):
        return taken_in_cls / in_cls_count

    def tn(taken_in_cls, taken_in_pop):
        return (
            ((pop_count - in_cls_count) - (taken_in_pop - taken_in_cls))
            / (pop_count - in_cls_count)
        )

    ab = '2'
    if ab in abbrev:
        tp_result = np.vectorize(tp)(in_cls_product[revised_indices])
        tn_result = np.vectorize(tn)(
            in_cls_product[revised_indices], total_product[revised_indices]
            )

    
        output = pd.DataFrame({
            'gene_1': gene_map[revised_indices[0]],
            'gene_2': gene_map[revised_indices[1]],
            'TP': tp_result,
            'TN': tn_result
            }, columns=['gene_1', 'gene_2', 'TP', 'TN'])

        return output
        
    else:
        tp_result = np.vectorize(tp)(in_cls_product[upper_tri_indices])
        tn_result = np.vectorize(tn)(
            in_cls_product[upper_tri_indices], total_product[upper_tri_indices]
            )
        output = pd.DataFrame({
            'gene_1': gene_map[upper_tri_indices[0]],
            'gene_2': gene_map[upper_tri_indices[1]],
            'TP': tp_result,
            'TN': tn_result
            }, columns=['gene_1', 'gene_2', 'TP', 'TN'])

        return output



def trips_hg(gene_map,in_cls_count,pop_count,trips_in_cls,trips_total,trips_indices,gene_1_mapped,gene_2_mapped,gene_3_mapped,abbrev,heur_limit):
    '''
    Altered version of pair_hg to do trips w/ the trips expression matrices
    See pair_hg for general var descriptions
    See combination product for full trips discussion
    Uses same vectorization scheme as pair, but has new rectangular indices and 
    third gene output. 'pair_indices' just gives us a list from 0 to however many values
    there are, this is then checked against the gene maps which contain all the mapping
    info for genes 1 & 2 (this setup plays nicely with the dataframe construction, the gene map
    values are already in the correct order). We then trim the matrix, disallowing any values 
    that have repeat genes (e.g. ACC -> no good) and also removing repeat triplets
    (e.g. ABC -> good , BAC -> no good if ABC already exists). This is SLOW, but more or less
    unavoidable. We didn't need to do this in the pair case b/c we knew ahead of time that all
    valid and unique values would be in the upper triangle, that doesn't hold in this case.
    '''


    def tp(taken_in_cls):
        return taken_in_cls / in_cls_count

    def tn(taken_in_cls, taken_in_pop):
        return (
            ((pop_count - in_cls_count) - (taken_in_pop - taken_in_cls))
            / (pop_count - in_cls_count)
        )


    abb = '3'
    if abb in abbrev:
        new_row_indices = []
        new_col_indices = []
        new_gene_1_mapped = []
        new_gene_2_mapped = []
        new_gene_3_mapped = []
        num_genes = len(gene_map)
        heuristic_limit = heur_limit


        #reassign all gene_maps and indices of computation for hg test
        #count -> position in the row to cut off the right half of the matrix
        #Need to cut out the mappings from the right half of the matrix

            
        #new col indices
        count = 0
        for num in trips_indices[1]:
            if count <= heuristic_limit-1:
                new_col_indices.append(num)
            elif count > heuristic_limit-1 and count < len(gene_map)-1:
                pass
            else:
                count = 0
                continue
            count = count + 1
        #new row indices
        count = 0
        for num in trips_indices[0]:
            if count <= heuristic_limit-1:
                new_row_indices.append(num)
            elif count > heuristic_limit-1 and count < len(gene_map)-1:
                pass
            else:
                count = 0
                continue
            count = count + 1
        #New gene 1 mapping
        count = 0
        for gene in gene_1_mapped:
            if count <= heuristic_limit-1:
                new_gene_1_mapped.append(gene)
            elif count > heuristic_limit-1 and count < len(gene_map)-1:
                pass
            else:
                count = 0
                continue
            count = count + 1
        #New gene 2 mapping 
        count = 0
        for gene in gene_2_mapped:
            if count <= heuristic_limit-1:
                new_gene_2_mapped.append(gene)
            elif count > heuristic_limit-1 and count < len(gene_map)-1:
                pass
            else:
                count = 0
                continue
            count = count + 1
        #New gene 3 mapping
        count = 0
        for gene in gene_3_mapped:
            if count <= heuristic_limit-1:
                new_gene_3_mapped.append(gene)
            elif count > heuristic_limit-1 and count < len(gene_map)-1:
                pass
            else:
                count = 0
                continue
            count = count + 1

        new_gene_1_mapped = pd.Index(new_gene_1_mapped)
        new_gene_2_mapped = pd.Index(new_gene_2_mapped)
        new_gene_3_mapped = pd.Index(new_gene_3_mapped)
        revised_row_indices = np.asarray(new_row_indices)
        revised_col_indices = np.asarray(new_col_indices)
        revised_indices = (revised_row_indices,revised_col_indices)

        
        tp_result = np.vectorize(tp)(trips_in_cls[revised_indices])
        tn_result = np.vectorize(tn)(
            trips_in_cls[revised_indices], trips_total[revised_indices]
        )
    
        vhg = np.vectorize(ss.hypergeom.sf, excluded=[1, 2, 4], otypes=[np.float])

    
        hg_result = vhg(
            trips_in_cls[revised_indices],
            pop_count,
            in_cls_count,
            trips_total[revised_indices],
            loc=1
        )
        print('HG + TP/TN done')
        pair_indices = []
        gene_count = int(len(gene_map))
        val_count =int( len(new_gene_1_mapped))
        for x in range(val_count):
            pair_indices.append(x)
        pair_indices = np.array(pair_indices)
        
        output = pd.DataFrame({
            'gene_1': new_gene_1_mapped[pair_indices],
            'gene_2': new_gene_2_mapped[pair_indices],
            'gene_3': new_gene_3_mapped[pair_indices],
            'HG_stat': hg_result,
            'TP' : tp_result,
            'TN': tn_result
        }, columns=['gene_1', 'gene_2', 'gene_3', 'HG_stat','TP','TN'])

    else:
        tp_result = np.vectorize(tp)(trips_in_cls[trips_indices])
        tn_result = np.vectorize(tn)(
            trips_in_cls[trips_indices], trips_total[trips_indices]
        )
    
        vhg = np.vectorize(ss.hypergeom.sf, excluded=[1, 2, 4], otypes=[np.float])

    
        hg_result = vhg(
            trips_in_cls[trips_indices],
            pop_count,
            in_cls_count,
            trips_total[trips_indices],
            loc=1
        )

        print('HG + TP/TN done')
        pair_indices = []
        gene_count = int(len(gene_map))
        val_count =int( len(gene_1_mapped))
        for x in range(val_count):
            pair_indices.append(x)
        pair_indices = np.array(pair_indices)
        output = pd.DataFrame({
            'gene_1': gene_1_mapped[pair_indices],
            'gene_2': gene_2_mapped[pair_indices],
            'gene_3': gene_3_mapped[pair_indices],
            'HG_stat': hg_result,
            'TP' : tp_result,
            'TN': tn_result
        }, columns=['gene_1', 'gene_2', 'gene_3', 'HG_stat','TP','TN'])



    
    output = output.sort_values(by='HG_stat', ascending=True)
    output = output.head(10000)
    #trims off values w/ repeating genes (e.g. ACC)
    #store unique trios in used_genes
    #iteratively check against the list
    used_genes = []
    counter=0
    prev_combo = []
    dropped=0
    buff = 0
    for index, row in output.iterrows():
        dropped = 0
        if counter == 2000:
            break
        if row[0]==row[1] or row[1]==row[2] or row[0]==row[2]:
            output.drop([index],inplace=True)
            continue
        for combo in list(reversed(prev_combo)):
            if row[0] in combo and row[1] in combo and row[2] in combo:
                output.drop([index],inplace=True)
                dropped = 1
                break
        if dropped != 1:
            prev_combo.append( [ row[0],row[1],row[2] ] )
            counter = counter + 1
    return output

