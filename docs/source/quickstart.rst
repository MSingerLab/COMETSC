.. _Github: https://github.com/MSingerLab/COMETSC

.. _Python: https://www.python.org/downloads/

.. _website: http://www.cometsc.com/index

.. _matplotlib:
https://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python

.. _tutorial: https://docs.python.org/3.6/tutorial/venv.html

Installation/Quickstart
========================

To skip installation on your own machine and run COMET on your data through COMET's interface, check out our website_!

To install COMET's Python implementation on your own machine, you should first have Python_ version 3.6 installed.

As with all python usage, it is recommended to use python virtual environments to avoid conflicting package versions. Find out how at python's documentation tutorial_.

The easiest way to install the package is to simply run in the terminal:

.. code-block:: console

   $ pip install COMETSC

This will auto-download the depencencies necessary for running the
tool. If you prefer to download the source code directly, you can clone the COMET source from our Github_, navigate to the cloned directory, then install it using Python.

Now, run COMET on your data. Give the files of your data as the first
three arguments and your desired output directory as your second argument.

In this Beta stage of COMET's development, your data must be formatted into three files in your input directory: 'markers.txt', 'cluster.txt', and 'tsne.txt'. See the :doc:`Manual<manual>` for more information.

In this example, we have our data located in the current directory. ``output/`` is the directory where COMET's output will be stored. 



.. code-block:: console
		
   $ hgmd marker_file tsne_file cluster_file output/

After this command is entered, COMET will run in the terminal, processing your data. See :doc:`Examples<examples>` for details on what this should look like.

The optional statements are described as follows:

-h ()

    Help.

-g (.txt)

    Optional gene list. Pointing -g to a text file with a list of
    genes will only run COMET with these genes. By default, it will
    use our own curated list of surface protein relevant genes. Gene
    lists should be either comma delimited or line delimited.

-C (int)

    Multi-process option. This allows the user to choose how many
    clusters to run in parallel. This option is entirely dependent on
    your hardware, a good benchmark is the number of cores on your
    computer is a good choice for this option if you want things to
    run faster. Defaults to 1.

-X (int)

    X-param of the XL-mHG test. Determines the percentage of ones
    above the threshold. See the statistical section for more
    info. Should only be values between 0 and 1.

-L (int)

    L-param of the XL-mHG test. Determines the lowest threshold to
    consider.
    See the statistical section for more info.

-Abbrev (list)

    Abbrev is a variable to turn on the heuristics for a given
    run. After typing '-Abbrev' , give a list of the number geene
    combinations you would like to be treated with the heuristic. For
    example, if you want to turn on the heuristic for 3- gene
    combinations only, input [3]. If you would like them for both 2-
    and 3- gene combinations, input [2,3].

-K (int)

    The number of gene combinations you would like to consider. This
    defaults to 2, can be at most 4 (for the time being). Note that
    the 4- gene combinations do not have a tested heuristic available.


-Down (int)

    For turning on downsampling. Will turn the data set into a random
    sampling of this many cells. Preserves cluster percentages by
    taking cells of each cluster at the same rate as the full set.
    
-Trim (int)

    Integer input for determining the size of the output
    files. Currently default is set to 2000, but an increase here will
    give more combinations / singletons in the final ranked output files.

-Count (Boolean)

    For turning on count data graphs. Will graph the log(count+1)
    instead of simply the count. Leads to nicer continuous expression
    plots.

-tenx (Boolean)

    For turning on tenx file integration. To submit a 10X format
    expression matrix, set this variable to True and in place of the
    expression matrix file in the command, use the folder containing
    the 10X data. Keep the default 10X names (e.g. barcodes.tsv,
    genes.tsv,matrix.mtx)

-online (Boolean)

    Online version of COMET, turning this on will produce a run as if
    it were submitted to the interface. This limits the run to looking
    at the 15 largest clusters.

    
Troubleshooting
========================    

There has been a known issue with the compatibility of matplotlib with
certain installations. A fix that has been successful in most test
cases can be found at matplotlib_ , consisting of changing the
matplotlibrc 'backend' variable to either 'Agg' or 'TkAgg'.

    
.. toctree::
