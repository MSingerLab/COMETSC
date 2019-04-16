.. |cluster| image:: _static/cluster_format.png

.. |markers| image:: _static/markers_format.png

.. |tsne| image:: _static/tsne_format.png


Input formatting
-----------------

Before you use COMET with gene expression data, your data should be
formatted into 3 files as shown below. Once the CSV's look good,
simply rename the extension to '.txt'.

* ``markers.txt``: a table in CSV format in a TXT file. The first row
  of the table lists cell names, while the first column
  lists genes (all caps) , the cell at the very top left should be
  blank. Each element in the rest of the table should contain a
  numerical gene expression value, corresponding to the row/cell and
  column/gene of the element. Comma delimited only.
  
  |markers|
  
* ``tsne.txt``: a table in CSV format, of three columns without column
  labels. The first column is cell name (the same as those in
  ``markers.txt``), the second is the tSNE_1 value for the cell, and
  the third is the tSNE_2 value for the cell. Comma delimited only.
  
  |tsne|
  
* ``cluster.txt``: a table in CSV format, of two columns without
  column labels. The first column is cell name (consistent with
  ``markers.txt`` and ``tsne.txt``), and the second is the cluster of
  which the cell is a member. Cluster numbers should start from 1 and
  Cluster names should avoid using punctuation marks. Comma delimited only.
  
  |cluster|
  

* ``gene_list``: A list of genes to use for filtering in your
  data. Should be comma delimited and in all caps.



.. toctree::
                                                                  
