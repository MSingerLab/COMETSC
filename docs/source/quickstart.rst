.. _Github: https://github.com/MSingerLab/COMETSC

.. _Python: https://www.python.org/downloads/release/python-362/

.. _website: http://www.cometsc.com/index

.. _matplotlib: https://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python

.. _tutorial: https://virtualenv.pypa.io/en/latest/

Installation/Quickstart
========================

To skip installation on your own machine and run COMET on your data through COMET's interface, check out our website_!

To install COMET's Python implementation on your own machine, you
should have Python_ version 3.6 downloaded & installed. If you do not, you will not
be able to set up the proper environment. 

As with all python usage, it is recommended to use python virtual
environments to avoid conflicting package versions. The following
commands can be used to set up the virtual environment assuming there
are no active environments running on the computer.

The first line will install the virtual environment software to your
computer (assuming you already have pip installed).

.. code-block:: console

   $ pip install virtualenv

The second line then sets up a new virtual environment directory of your choice with
the proper version of python

.. code-block:: console

   $ virtualenv new_dir --python=python3.6

The following line will activate the virtual environment, just be sure that the
path used is the correct one with respect to the current working directory.
   
 .. code-block:: console

   $ source new_dir/bin/activate

Now, install the proper COMET version to the virtual environment.

.. code-block:: console
 
   $ pip install COMETSC

Now run COMET on your data!  For example:
   
.. code-block:: console

   $ Comet markers.txt tsne.txt cluster.txt output/

When you are finished using the virtual environment,
simply type 'deactivate' and you will return to your usual python.
   
.. code-block:: console

   $ deactivate

	
The easiest way to install the package without the use of virtual
environments is to simply run in the terminal:

.. code-block:: console

   $ pip install COMETSC

This will auto-download the depencencies necessary for running the
tool. If you prefer to download the source code directly, you can
clone the COMET source from our Github_, navigate to the cloned
directory, then install it using Python. If you are having issues
setting up the virtual environment, have a look at the documentation
at python's tutorial_ .

Now, run COMET on your data. Give the files of your data as the first
three arguments and your desired output directory as your second argument.

Your data must be formatted into three files in your input directory: 'markers.txt', 'cluster.txt', and 'vis.txt'. See the :doc:`Manual<usage>` for more information.

In this example, we have our data located in the current directory. ``output/`` is the directory where COMET's output will be stored. 



.. code-block:: console
		
   $ Comet marker_file vis_file cluster_file output/

After this command is entered, COMET will run in the terminal,
processing your data. See some :doc:`examples<Output>` for details on
what this should look like.

Download some example inputs :download:`here <_static/example_ins.zip>` !

The optional statements are described as follows:

-h ()

    Help.

-g (.txt)

    Optional gene list. Pointing -g to a text file with a list of
    genes will only run COMET with these genes. By default, it will
    run on all genes. Gene lists should be either comma delimited or line delimited.

-C (int)

    Multi-process option. This allows the user to choose how many
    clusters to run in parallel. This option is entirely dependent on
    your hardware, a good benchmark is the number of cores on your
    computer is a good choice for this option if you want things to
    run faster. Defaults to 1.

-X (int)

    X-param of the XL-mHG test. Determines the percentage of ones
    above the threshold. See the :doc:`statistical<details>` section for more
    info. Should only be values between 0 and 1. Defaults to .15 (15%)

-L (int)

    L-param of the XL-mHG test. Determines the lowest threshold to
    consider.
    See the :doc:`statistical<details>` section for more info. Defaults to 2 * size of cluster

-Abbrev (list)

    Abbrev is a variable to turn on the heuristics for a given
    run. After typing '-Abbrev' , give a list of the number gene
    combinations you would like to be treated with the heuristic. For
    example, if you want to turn on the heuristic for 3- gene
    combinations only, input [3]. If you would like them for both 2-
    and 3- gene combinations, input [2,3]. Defaults to an empty list.

-K (int)

    The number of gene combinations you would like to consider. This
    defaults to 2, can be at most 4 (for the time being). For example,
    if 4 is chosen, COMET will run 2, 3, and 4 gene combinations.


-Down (int)

    For turning on downsampling. Will turn the data set into a random
    sampling of this many cells. Preserves cluster percentages by
    taking cells of each cluster at the same rate as the full
    set. Defaults to no downsampling.
    
-Trim (int)

    Integer input for determining the size of the output
    files. Currently default is set to 2000 lines, but an increase here will
    give more combinations / singletons in the final ranked output files.

-Count (Boolean)

    Will graph the log(expression+1). Leads to
    nicer continuous expression plots when count data is used.
    Defaults to False, will plot the inputted expression as is.

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
