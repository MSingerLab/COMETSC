import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from itertools import repeat
import time
import os
import math
from adjustText import adjust_text

"""
Set of modularized visualization functions for COMET; producing graphs in PDFs.
"""

"""
Goals: Separate singleton and pair graphs for:

    - TP/TN with gene names

    - Side-by-side comparison between discrete and continuous expression

For singletons and pairs combined:

    - Side-by-side comparison of single gene discrete expression and resulting
      pair discrete expression.
"""


CMAP_CONTINUOUS = cm.get_cmap('nipy_spectral')
#CMAP_CONTINUOUS = cm.get_cmap('viridis')
CMAP_DISCRETE = cm.get_cmap('bwr')

def make_plots(
    pair, sing, sing_tp_tn, xlmhg, trips, quads_fin, tsne, discrete_exp, marker_exp, plot_pages,
    combined_path, sing_combined_path, discrete_path, tptn_path,trips_path,quads_path,sing_tptn_path,count_data
):
    """
    General function for all visualization generation.  Arguments should be
    self-explanatory.  See __main__.
    """
    # cutoff maps genes to their (absolute) cutoff. I.e. complements are mapped
    # to positive cutoffs.
    #print(discrete_exp)
    #print(marker_exp)
    cutoff = sing['cutoff_val'].abs()
    rank = sing['rank']
    TP_pair = pair['TP']
    TN_pair = pair['TN']
    TP_sing = sing['TP']
    TN_sing = sing['TN']
    xlmhg.set_index('gene_1',inplace=True)
    pair_sing_only = pair[pair['gene_2'].isnull()]
    sing_rank = pd.Series(
        pair_sing_only['rank'].values, index=pair_sing_only['gene_1']
    )

    p_short = pair[pair['Plot']==1]
    #p_short = pair.iloc[:plot_pages]
    s_short = sing[sing['Plot']==1]
    #s_short = sing.iloc[:plot_pages]

    if count_data == 1:
        marker_exp = marker_exp.applymap(lambda x: math.log(abs(x)+1,2))
        #marker_exp = marker_exp.applymap(lambda x: min(x,1500))
    if type(trips) == int:
        pass
    else:
        t_short = trips.iloc[:plot_pages]
    if type(quads_fin) == int:
        pass
    else:
        q_short=quads_fin.iloc[:plot_pages]
    try:
        vmt = np.vectorize(make_title,excluded=[6,7])
        d_plot_genes = zip(
            zip(
                vmt(
                    p_short['gene_1'], p_short['gene_2'],
                    p_short['rank'], p_short['gene_1'].map(cutoff),
                    p_short['TP'], p_short['TN'], sing_tp_tn,xlmhg
                    ), vmt(
                        p_short['gene_1'], np.nan,
                        p_short['gene_1'].map(rank),
                        p_short['gene_1'].map(cutoff),
                        p_short['gene_1'].map(TP_sing), p_short['gene_1'].map(TN_sing), sing_tp_tn,xlmhg
                        ), vmt(
                            p_short['gene_2'], np.nan,
                            p_short['gene_2'].map(rank),
                            p_short['gene_2'].map(cutoff),
                            p_short['gene_2'].map(TP_sing), p_short['gene_2'].map(TN_sing), sing_tp_tn,xlmhg
                            )
                ), p_short['gene_1'].values, p_short['gene_2'].values
            )
    except:
        print('No genes to plot')
    if type(trips) == int:
        pass
    else:
        try:
            vmt_2 = np.vectorize(make_trips_title)
            t_plot_genes = zip(
                zip(
                    vmt_2(
                        t_short['gene_1'], t_short['gene_2'], t_short['gene_3'],
                        t_short['rank'], t_short['gene_1'].map(cutoff)
                        ), vmt_2(
                            t_short['gene_1'], np.nan, np.nan,
                            t_short['gene_1'].map(rank),
                            t_short['gene_1'].map(cutoff)
                            ), vmt_2(
                                t_short['gene_2'], np.nan, np.nan,
                                t_short['gene_2'].map(rank),
                                t_short['gene_2'].map(cutoff)
                                ), vmt_2(
                                    t_short['gene_3'], np.nan, np.nan,
                                    t_short['gene_3'].map(rank),
                                    t_short['gene_3'].map(cutoff)
                                    )
                    ), t_short['gene_1'].values, t_short['gene_2'].values, t_short['gene_3'].values
                )
        except:
            print('No genes to plot')
    

    if type(quads_fin) == int:
        pass
    else:
        vmt_3 = np.vectorize(make_quads_title)
        q_plot_genes = zip(
            zip(
                vmt_3(
                    q_short['gene_1'], q_short['gene_2'], q_short['gene_3'], q_short['gene_4'],
                    q_short['rank'], q_short['gene_1'].map(cutoff)
                    ), vmt_3(
                        q_short['gene_1'], np.nan, np.nan, np.nan,
                        q_short['gene_1'].map(rank),
                        q_short['gene_1'].map(cutoff)
                        ), vmt_3(
                            q_short['gene_2'], np.nan, np.nan, np.nan,
                            q_short['gene_2'].map(rank),
                            q_short['gene_2'].map(cutoff)
                            ), vmt_3(
                                q_short['gene_3'], np.nan, np.nan, np.nan,
                                q_short['gene_3'].map(rank),
                                q_short['gene_3'].map(cutoff)
                                ), vmt_3(
                                q_short['gene_3'], np.nan, np.nan, np.nan,
                                q_short['gene_3'].map(rank),
                                q_short['gene_3'].map(cutoff)
                                )
                ), q_short['gene_1'].values, q_short['gene_2'].values, q_short['gene_3'].values, q_short['gene_4'].values)


    
    print("Drawing discrete plots for pairs...")
    try:
        make_discrete_plots(
            tsne, discrete_exp, d_plot_genes, discrete_path, 2
        )
    except Exception as err:
        print('No plots generated')
        print(err)
    if type(trips) == int:
        pass
    else:
        print("Drawing discrete plots for trips...")
        try:
            make_discrete_plots(
                tsne, discrete_exp, t_plot_genes, trips_path, 3
                )
        except Exception as err:
            print('no plots generated')
            print(err)
    
    if type(quads_fin) == int:
        pass
    else:
        print("Drawing discrete plots for quads...")
        try:
            make_discrete_plots(
                tsne, discrete_exp, q_plot_genes, quads_path, 4
                )
        except Exception as err:
            print('no plots generated')
            print(err)
    try:
        c_plot_genes = zip(
            zip(
                vmt(
                    p_short['gene_1'], np.nan,
                    p_short['gene_1'].map(rank),
                    p_short['gene_1'].map(cutoff),
                    p_short['gene_1'].map(TP_sing),
                    p_short['gene_1'].map(TN_sing),sing_tp_tn,xlmhg
                    ), vmt(
                        p_short['gene_2'], np.nan,
                        p_short['gene_2'].map(rank),
                        p_short['gene_2'].map(cutoff),
                        p_short['gene_2'].map(TP_sing),
                        p_short['gene_2'].map(TN_sing),sing_tp_tn,xlmhg
                        )
                ), p_short['gene_1'].values, p_short['gene_2'].values
            )
    except:
        print('no genes to plot')
    print("Drawing combined plots...")
    try:
        make_combined_plots(
            tsne, discrete_exp, marker_exp, c_plot_genes, combined_path
            )
    except Exception as err:
        print('no plots generated')
        print(err)
    c_s_plot_genes = zip(
        zip(
            vmt(
                s_short.index, np.nan,
                s_short['rank'], s_short['cutoff_val'],
                s_short['TP'], s_short['TN'],sing_tp_tn,xlmhg
            ), repeat(np.nan)
        ), s_short.index, repeat(np.nan)
    )
    print("Drawing singleton combined plots...")
    try:
        make_combined_plots(
            tsne, discrete_exp, marker_exp, c_s_plot_genes, sing_combined_path
            )
    except Exception as err:
        print('no plots generated')
        print(err)
    pair_tp_tn = pair[['gene_1', 'gene_2', 'TP', 'TN']]
    sing_tp_tn = sing[['TP', 'TN']]
    print("Drawing true positive/negative plots...")
    try:
        make_tp_tn_plot(
            zip(p_short['gene_1'], p_short['gene_2']),
            sing_tp_tn, pair_tp_tn, tptn_path, 0
            )
        make_tp_tn_plot(
            zip(s_short.index, repeat(np.nan)),
            sing_tp_tn, pair_tp_tn, sing_tptn_path, 1
            )
    except Exception as err:
        print('no plots generated')
        print(err)

def make_title(gene_1, gene_2, rank, cutoff_val, TP, TN, sing_tp_tn,xlmhg):
    """Makes a plot title for a gene or gene pair.

    Formatting: for pairs, 'rank $rank: $gene_1+$gene_2', and for singletons,
    'rank $rank: $gene_1 $cutoff_val'.  gene_2 should be None for singletons.

    :param genes: A DataFrame with columns 'gene_1', 'gene_2', 'rank',
        'cutoff_val'.

    :returns: A list of strings containing the titles, with indices
              corresponding to the input DataFrame.

    :rtype: string list
    """
    if np.isnan(rank):
        rank = 'Filtered'
        cutoff_val = xlmhg.at[gene_1,'cutoff_val']
        TP = sing_tp_tn.at[gene_1,'TP']
        TN = sing_tp_tn.at[gene_1,'TN']
    if pd.isna(gene_2):
        try:
            return ("rank %.0f: %s %.3f - TP:%.2f TN:%.2f" % (rank, gene_1, cutoff_val, TP, TN))
        except:
            return ("rank %s: %s %.3f - TP:%.2f TN:%.2f" % (rank, gene_1, cutoff_val, TP, TN))
    else:
        return ("rank %.0f: %s+%s - TP:%.2f TN:%.2f" % (rank, gene_1, gene_2, TP, TN))

def make_trips_title(gene_1,gene_2,gene_3,rank,cutoff_val):
    
    if pd.isna(gene_2):
        return ("rank %.0f: %s %.3f" % (rank, gene_1, cutoff_val))
    elif pd.isna(gene_3):
        return ("rank %.0f: %s+%s" % (rank, gene_1, gene_2))
    else:
        return ("rank %.0f: %s+%s+%s" % (rank, gene_1, gene_2, gene_3)) 

def make_quads_title(gene_1,gene_2,gene_3,gene_4,rank,cutoff_val):

    
    if pd.isna(gene_2):
        return ("rank %.0f: %s %.3f" % (rank, gene_1, cutoff_val))
    elif pd.isna(gene_3):
        return ("rank %.0f: %s+%s" % (rank, gene_1, gene_2))
    elif pd.isna(gene_4):
        return ("rank %.0f: %s+%s+%s" % (rank, gene_1, gene_2, gene_3))
    else:
        return ("rank %.0f: %s+%s+%s+%s" % (rank, gene_1, gene_2, gene_3, gene_4))
    
def make_plot(ax, title, coords, cmap, draw_cbar=False):
    """
    Make a single graph on ax with given specs.  Plots only absolute values.
    """
    ax.set_title(title)#,fontdict={'fontsize':6})
    ax.set_xlabel('tSNE_1')
    ax.set_ylabel('tSNE_2')
    sc = ax.scatter(
        x=coords[0],
        y=coords[1],
        c=abs(coords[2]),
        s=4,
        cmap=cmap
    )
    if draw_cbar:
        plt.colorbar(sc, ax=ax)


def make_discrete_plots(tsne, discrete_exp, plot_genes, path,num):
    """Plots discrete gene expression of paired genes to PDF.

    For each gene pair listed in plot_genes, make three scatterplots.  First, a
    plot showing joint expression.  Then, two plots showing singleton
    expression for each genes.  If a single gene is passed, plot only its
    expression, and make two blank plots.  Save each gene/gene pair as a PDF
    page, then save to path.

    :param tsne: A DataFrame with 'cell', 'tSNE_1', and 'tSNE_2' columns.
    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param plot_genes: A list of 3-tuples, where the first element of each
        tuple is another 3-tuple containing the three plot titles to be used.
        The other 2 elements are the gene names to be plotted.
    :param path: The path to which the PDF will be saved.

    :returns: Nothing.
    """
    def make_quads_discrete_page(fig, ax_triple, titles, gene_1, gene_2, gene_3,gene_4):
        """Make page with quads discrete plots given titles and genes."""
        coords_df = tsne.merge(discrete_exp[[gene_1, gene_2, gene_3,gene_4]], on='cell')
        coords_df['pair'] = coords_df[gene_1] * coords_df[gene_2]
        coords_df['trips'] = coords_df[gene_1] * coords_df[gene_2] * coords_df[gene_3]
        coords_df['quads_fin'] = coords_df[gene_1] * coords_df[gene_2] * coords_df[gene_3] * coords_df[gene_4]

        for (graph_index, z_label) in ((0, 'quads_fin'), (1, gene_1), (2, gene_2), (3, gene_3), (4, gene_4)):
            make_plot(
                ax=ax_triple[graph_index],
                title=titles[graph_index],
                coords=(
                    coords_df['tSNE_1'].values,
                    coords_df['tSNE_2'].values,
                    coords_df[z_label].values
                ),
                cmap=CMAP_DISCRETE
            )

    
    def make_trips_discrete_page(fig, ax_triple, titles, gene_1, gene_2, gene_3):
        """Make page with trips discrete plots given titles and genes."""
        coords_df = tsne.merge(discrete_exp[[gene_1, gene_2, gene_3]], on='cell')
        coords_df['pair'] = coords_df[gene_1] * coords_df[gene_2]
        coords_df['trips'] = coords_df[gene_1] * coords_df[gene_2] * coords_df[gene_3]

        for (graph_index, z_label) in ((0, 'trips'), (1, gene_1), (2, gene_2), (3, gene_3)):
            make_plot(
                ax=ax_triple[graph_index],
                title=titles[graph_index],
                coords=(
                    coords_df['tSNE_1'].values,
                    coords_df['tSNE_2'].values,
                    coords_df[z_label].values
                ),
                cmap=CMAP_DISCRETE
            )
        
    def make_pair_discrete_page(fig, ax_triple, titles, gene_1, gene_2):
        """Make page with three discrete plots given titles and genes."""
        coords_df = tsne.merge(discrete_exp[[gene_1, gene_2]], on='cell')
        coords_df['pair'] = coords_df[gene_1] * coords_df[gene_2]

        for (graph_index, z_label) in ((0, 'pair'), (1, gene_1), (2, gene_2)):
            make_plot(
                ax=ax_triple[graph_index],
                title=titles[graph_index],
                coords=(
                    coords_df['tSNE_1'].values,
                    coords_df['tSNE_2'].values,
                    coords_df[z_label].values
                ),
                cmap=CMAP_DISCRETE
            )

    def make_single_discrete_page(fig, ax_triple, title, gene):
        """Make page with one discrete plot given title and gene"""
        print(title)
        coords_df = tsne.merge(discrete_exp[[gene]], on='cell')
        make_plot(
            ax=ax_triple[0],
            title=title[0],
            coords=(
                coords_df['tSNE_1'].values,
                coords_df['tSNE_2'].values,
                coords_df[gene].values
            ),
            cmap=CMAP_DISCRETE
        )

    #with PdfPages(path) as pdf:
    count = 1
    try:
        os.makedirs(path)
    except:
        os.system('rm -r ' + path)
        os.makedirs(path)
    for plot_gene in plot_genes:
        # print(plot_gene)
        if num == 4:
            fig, ax_triple = plt.subplots(ncols=5, figsize=(15, 5))
            if pd.isnull(plot_gene[2]):
                make_single_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    title=plot_gene[0],
                    gene=plot_gene[1]
                    )
            elif pd.isnull(plot_gene[3]):
                make_pair_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    titles=plot_gene[0],
                    gene_1=plot_gene[1],
                    gene_2=plot_gene[2]
                    )
            elif pd.isnull(plot_gene[4]):
                make_trips_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    titles=plot_gene[0],
                    gene_1=plot_gene[1],
                    gene_2=plot_gene[2],
                    gene_3=plot_gene[3]
                    )
            else:
                make_quads_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    titles=plot_gene[0],
                    gene_1=plot_gene[1],
                    gene_2=plot_gene[2],
                    gene_3=plot_gene[3],
                    gene_4=plot_gene[4]
                    )
        elif num == 3:
            fig, ax_triple = plt.subplots(ncols=4, figsize=(15, 5))
            if pd.isnull(plot_gene[2]):
                make_single_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    title=plot_gene[0],
                    gene=plot_gene[1]
                    )
            elif pd.isnull(plot_gene[3]):
                make_pair_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    titles=plot_gene[0],
                    gene_1=plot_gene[1],
                    gene_2=plot_gene[2]
                    )
            else:
                make_trips_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    titles=plot_gene[0],
                    gene_1=plot_gene[1],
                    gene_2=plot_gene[2],
                    gene_3=plot_gene[3]
                    )
        else:
            fig, ax_triple = plt.subplots(ncols=3, figsize=(15, 5))
            if pd.isnull(plot_gene[2]):
                make_single_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    title=plot_gene[0],
                    gene=plot_gene[1]
                    )
            else:
                make_pair_discrete_page(
                    fig=fig, ax_triple=ax_triple,
                    titles=plot_gene[0],
                    gene_1=plot_gene[1],
                    gene_2=plot_gene[2]
                    )
                
            #pdf.savefig(fig)
        if len(plot_gene) >= 4:
            plt.savefig( path + '/' + 'rank_' + str(count) + '.png')
        else:
            try:
                if ':' in str(plot_gene[0][0][:9]):
                    if ':' in str(plot_gene[0][0][:8]):
                        if ':' in str(plot_gene[0][0][:7]):
                            plt.savefig( path + '/' + str(plot_gene[0][0][:6]) + '.png')
                        else:
                            plt.savefig( path + '/' + str(plot_gene[0][0][:7]) + '.png')
                    else:
                        plt.savefig( path + '/' + str(plot_gene[0][0][:8]) + '.png')
                else:
                    plt.savefig( path + '/' + str(plot_gene[0][0][:9]) + '.png')
            except:
                plt.savefig( path + '/' + 'rank_' + str(count) + '.png')
        count=count+1
        plt.close(fig)

def make_combined_plots(tsne, discrete_exp, marker_exp, plot_genes, path):
    """Plots discrete alongside continuous expression to PDF.

    For each gene/gene pair listed in plot_genes, make two scatterplots: a plot
    showing discrete expression, and a plot showing expression on a color
    spectrum.  For gene pairs, make these two plots separately for each gene.
    Save each gene/gene pair as a PDF page, then save to path.

    :param tsne: A DataFrame with 'cell', 'tSNE_1', and 'tSNE_2' columns.
    :param discrete_exp: A DataFrame whose rows are cell identifiers, columns
        are gene identifiers, and values are boolean values representing gene
        expression.
    :param marker_exp: A DataFrame whose rows are cell identifiers, columns are
        gene identifiers, and values are float values representing gene
        expression.
    :param plot_genes: A list of 3-tuples, where the first element of each
        tuple is another 2-tuple containing the two titles to be used..
        The other 2 elements are the gene names to be plotted.
    :param path: The path to which the PDF will be saved.

    :returns: Nothing.
    """
    def make_pair_combined_page(fig, axes, titles, gene_1, gene_2):
        """Make page with two pairs of plots for given genes."""
        disc_coords = tsne.merge(discrete_exp[[gene_1, gene_2]], on='cell')
        cont_coords = tsne.merge(marker_exp[[gene_1, gene_2]], on='cell')
        for (graph_index, z_label) in ((0, gene_1), (1, gene_2)):
            make_plot(
                ax=axes[graph_index][0], title=titles[graph_index],
                coords=(
                    disc_coords['tSNE_1'].values,
                    disc_coords['tSNE_2'].values,
                    disc_coords[z_label].values
                ),
                cmap=CMAP_DISCRETE
            )
            make_plot(
                ax=axes[graph_index][1], title=str(z_label),
                coords=(
                    cont_coords['tSNE_1'].values,
                    cont_coords['tSNE_2'].values,
                    cont_coords[z_label].values
                ),
                cmap=CMAP_CONTINUOUS, draw_cbar=True
            )

    def make_single_combined_page(fig, title, axes, gene):
        """Make page with single pair of plot of given gene."""
        disc_coords = tsne.merge(discrete_exp[[gene]], on='cell')
        cont_coords = tsne.merge(marker_exp[[gene]], on='cell')
        make_plot(
            ax=axes[0], title=title,
            coords=(
                disc_coords['tSNE_1'].values,
                disc_coords['tSNE_2'].values,
                disc_coords[gene].values
            ),
            cmap=CMAP_DISCRETE
        )
        #if str(gene)[-8:] == 'NEGATION':
        make_plot(
            ax=axes[1], title=str(gene),
            coords=(
                cont_coords['tSNE_1'].values,
                cont_coords['tSNE_2'].values,
                cont_coords[gene].values
            ),
            cmap=CMAP_CONTINUOUS, draw_cbar=True
        )

    #with PdfPages(path) as pdf:
    counter = 1
    try:
        os.makedirs(path)
    except:
        os.system('rm -r ' + path)
        os.makedirs(path)
    for plot_gene in plot_genes:
        if pd.isnull(plot_gene[2]):
            fig, axes = plt.subplots(ncols=2, figsize=(10, 5))
            make_single_combined_page(
                fig, plot_gene[0][0], axes, plot_gene[1]
            )
        else:
            fig, axes = plt.subplots(
                nrows=2, ncols=2, figsize=(10, 10)
            )
            make_pair_combined_page(
                fig, axes, plot_gene[0], plot_gene[1], plot_gene[2]
            )
        #pdf.savefig(fig)
        plt.savefig( path + '/' + 'rank_' + str(counter) + '.png')
        counter=counter+1
        plt.close(fig)


def make_tp_tn_plot(plot_genes, sing_tp_tn, pair_tp_tn, path,sing_pair):
    """Plots TP/TN rates of genes/pairs to PDF.

    For each gene/gene pair listed in plot_genes, plot their TP/TN rate on a
    scatterplot, labeling the point with the gene/gene pair name.  When done,
    output this scatterplot to PDF and save to path.

    :param plot_genes: An array whose elements are tuples representing gene
        pairs.  If the second element is empty, it represents a singleton.
    :param sing_tp_tn: A DataFrame with ''TP', and 'TN' columns, with gene
        indices
    :param pair_tp_tn: A DataFrame with 'gene_1', 'gene_2', 'TP', and 'TN'
        columns.
    :param path: The path to which the PDF will be saved.

    :returns: Nothing.
    """
    PADDING = 0.002

    fig = plt.figure(figsize=[15, 15])
    plt.xlabel("True positive")
    plt.ylabel("True negative")
    plt.title("True positive/negative")
    plt.axis([0.0, 1.1, 0.0, 1.1])

    def get_data(genes):
        if pd.isnull(genes[1]):
            title = genes[0]
            data_row = sing_tp_tn.loc[genes[0]]
        else:
            title = genes[0] + "+" + genes[1]
            data_row = pair_tp_tn[
                ((pair_tp_tn['gene_1'] == genes[0]) &
                 (pair_tp_tn['gene_2'] == genes[1]))
                | ((pair_tp_tn['gene_1'] == genes[1]) &
                   (pair_tp_tn['gene_2'] == genes[0]))
            ]
        return [title, data_row['TP'], data_row['TN']]

    coords_df = pd.DataFrame()
    data = list(map(get_data, list(plot_genes)))
    coords_df[['title', 'TP', 'TN']] = pd.DataFrame(
        data, columns=['title', 'TP', 'TN']
    )

    if sing_pair == 1:
        new_list = [['zero',0,0]]
        for index, row in coords_df.iterrows():
            new_list.append([row['title'],row['TP'],row['TN']])
        new_df = pd.DataFrame(new_list, columns = ['title','TP','TN'])    
        #print(new_df)
        concat_new_df = new_df.head(50)
        #concat_new_df.drop(0,inplace=True)
        plt.scatter(concat_new_df['TP'], concat_new_df['TN'], s=3,c='red')
        #print(coords_df)
        #time.sleep(10000)
        texts = []
        for x,y,s in zip(concat_new_df['TP'],concat_new_df['TN'],concat_new_df['title']):
            texts.append(plt.text(x, y, s, size=7))
        adjust_text(texts,arrowprops=dict(arrowstyle="->", color='blue', lw=0.5))


        
    else:
        new_list = [['zero',0,0]]
        for index, row in coords_df.iterrows():
            TP = row['TP'].tolist()[-1]
            TN = row['TN'].tolist()[-1]
            new_list.append([row['title'],TP,TN])
        new_df = pd.DataFrame(new_list, columns = ['title','TP','TN'])    
        #print(new_df)
        concat_new_df = new_df.head(50)
        #concat_new_df.drop(0,inplace=True)
        plt.scatter(concat_new_df['TP'], concat_new_df['TN'], s=3,c='red')
        #print(coords_df)
        #time.sleep(10000)
        texts = []
        for x,y,s in zip(concat_new_df['TP'],concat_new_df['TN'],concat_new_df['title']):
            texts.append(plt.text(x, y, s, size=7))
        adjust_text(texts,arrowprops=dict(arrowstyle="->", color='blue', lw=0.5))
        
    #for index, row in coords_df.iterrows():
    #    plt.annotate(row['title'], (row['TP'] + PADDING, row['TN'] + PADDING))
    fig.savefig(path + '.png')
    plt.close(fig)
