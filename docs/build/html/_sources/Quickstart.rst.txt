.. _Github: https://github.com/aaronyaosmith/HG_marker_detection

.. _Python: https://www.python.org/downloads/

.. _website: http://www.cometsc.com/login

Installation/Quickstart
========================

To skip installation on your own machine and run COMET on your data on COMET's servers, check out our website_!

To install COMET's Python implementation on your own machine, you should first have Python_ version 3.6 installed.

The easiest way to install the package is to simply run in the terminal:

.. code-block:: console

   $ pip install HGMarkerDetection

This will auto-download the depencencies necessary for running the
tool. If you prefer to download the source code directly, you can clone the COMET source from our Github_, navigate to the cloned directory, then install it using Python.

Now, run COMET on your data. Give the files of your data as the first
three arguments and your desired output directory as your second argument.

In this Beta stage of COMET's development, your data must be formatted into three files in your input directory: 'markers.txt', 'cluster.txt', and 'tsne.txt'. See the :doc:`Manual<manual>` for more information.

In this example, we have our data located in the current directory. ``output/`` is the directory where COMET's output will be stored. 



.. code-block:: console
		
   $ hgmd marker_file tsne_file cluster_file output/

   $ hgmd [-h] [-g [G]] [-C [C]] [-X [X]] [-L [L]] [-Abbrev [ABBREV]]
            [-K [K]]
After this command is entered, COMET will run in the terminal, processing your data. See :doc:`Output<Output>` for details on what this should look like.

The optional statements are described as follows:

-h
    Help.

-g (.txt)

    Optional gene list. Pointing -g to a text file with a list of
    genes will only run COMET with these genes. By default, it will
    use our own curated list of surface protein relevant genes.

-C (int)

    Multi-process option. This allows the user to choose how many
    clusters to run in parallel. This option is entirely dependent on
    your hardware, a good benchmark is the number of cores on your
    computer is a good choice for this option if you want things to
    run faster. Defaults to 1.

-X (int)

    X-param of the XL-mHG test. Determines the percentage of ones
    above the threshold. See the statistical section for more info.

-L (int)

    L-param of the XL-mHG test. Determines the lowest threshold to consider. See the statistical section for more info.

-Abbrev (Boolean)

    Abbrev is only used when dealing with 3-gene combinations. The
    full treatment of 3-gene combinations is computationally
    expensive, so if you are using anything other than a computing
    cluster this should be set to True (True by default). Allows
    3-gene combinations to be run on a laptop, but drastically reduces
    the space.

-K (int)

    The number of gene combinations you would like to consider. This
    defaults to 2, can be at most 3 (for the time being). Should only
    be either 2 or 3. Will be extended to 4 gene combinations in the
    future.



Common Installation Problems:
========================
Matplotlib failing -> Sometimes the bare installation of matplotlib
can cause problems with the 'backend' setting. Maneuvering to the file
'matplotlibrc' in the site packages and changing the backend variable
to 'Agg' will fix this problem. Unrelated to COMET but common issue we
ran into.

General python -> Most python problems are fixed by installing
virtualenv and setting up a virtual environment for your project with
python 3.6.X . Check out https://docs.python-guide.org/dev/virtualenvs/


.. toctree::
