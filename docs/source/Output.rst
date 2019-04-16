.. |pair_csv| image:: _static/pair_csv.png
.. |sing_csv| image:: _static/sing_csv.png
.. |combined| image:: _static/cluster_Spleen_1_combined-01.jpg
.. |discrete| image:: _static/cluster_Spleen_1_discrete-01.jpg
.. |sing_cmb| image:: _static/cluster_Spleen_1_singleton_combined-01.jpg
.. |sing_TPTN| image:: _static/cluster_Spleen_1_singleton_TP_TN-1.jpg
.. |TPTN| image:: _static/cluster_Spleen_1_TP_TN-1.jpg

Output
==========

Following are examples of COMET output. Descriptions of statistical values can be found in the :doc:`details<details>` section of the manual.

**Data:**

*CSV  pair output* : 
Gives the gene-pairs ordered based on our ranking system of
statistical relevance. Column 'rank' gives the final rank of the pair
compared to the others. See our discussion in the manual for more info
on ranking and the various statistics calculated for each pair.
|pair_csv|

*Singleton-only CSV output* : 
Gives the single-gene marker list, ordered based on statistical relevance.
|sing_csv|                                                                       


**Visualizations:**

Shown below is t-sne format, but any two-dimensional visualization method is fine to use as input to COMET.



*Combined continuous/discrete plots* : 
Gives the discrete and continuous plot for the single genes of the top
performing pairs.
|combined| 

*Discrete plots comparing combinations/singletons* : 
Discrete only plots that show a pair of genes each. Gives a
visualization of the discrete dual-expression alongside each single gene's
discrete expression.
|discrete|

*Singleton-only combined plots* : 
Plots that show the discrete and continuous versio  of the top
performing single genes.
|sing_cmb| 

*True positive/negative plot* : 
True positive/ True negative values for the top gene pairs.
|TPTN| 

*Singleton-only TP/TN plot* : 
true positive / True negatives values for the top single genes.
|sing_TPTN| 



    
.. toctree::
