.. _using:

Using RIA
=========

The program RIA takes a PLINK binary file as input (.bed/.bim/.fam) and produces a results file of the analysis. Basic usage of the program is given by typing:

.. code-block:: none

    ./ria -i mydata.bed


The most likely options that will need to be used are how to run the programs PLINK and KING:

.. code-block:: none

    ./ria -plink /path/to/plink/plink -king /path/to/king/king -i mydata.bed


Typing `ria` with no options will output usage details:


.. code-block:: none

    RIA: Regional IBD Analysis, v1.1
    ------------------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Research Software Engineering, Newcastle University

    Usage:
    ./ria [options] -i pedigree.bed
    or ./ria -pf parameterfile [pedigree.bed]

    Options:
    -window-size-cm n         -- set window size to n centimorgans
    -decreasing-cm            -- set decreasing centimorgan values as previous valid value
    -window-min-snps m        -- set required minimium number of SNPs in a window to m
    -step-size s              -- step size of windows, s
    -start-snp a              -- start analysis from SNP number a
    -end-snp b                -- end analysis at SNP number b
    -start-snp-name x         -- start analysis from SNP name x
    -end-snp-name y           -- end analysis at SNP name y
    -job d t                  -- job number d of t
    -i file.bed               -- input binary pedigree file, file.bed
    -o results.dat            -- output results file, results.dat
    -i-prior file             -- input prior IBDs, file
    -o-prior file             -- output prior IBDs, file
    -prior-only               -- calculate prior IBDs only
    -plink command            -- command used to run PLINK
    -i-posteriors-prefix fp   -- input posterior IBDs file prefix
    -o-posteriors-prefix fp   -- output posterior IBDs file prefix
    -posterior-start-window p -- start analysis from posterior window number p
    -posterior-end-window q   -- end analysis at posterior window number q
    -plink-options "ops"      -- PLINK pruning options used to calculate the prior
    -king command             -- command used to run KING
    -log results.log          -- log filename, results.log
    -ndv                      -- no dominance variance
    -so                       -- suppress output to screen

    Default Options in Effect:
    -window-size-cM 15
    -window-min-snps 100
    -step-size 50
    -plink plink
    -king king
    -o riaResults.dat



See :ref:`example` for an example of how to use RIA to analyse some data.

.. _parameterfile:

Parameter file
--------------

A parameter file, `.pf`, may be used with RIA instead of writing all of the options on the command line. To use a parameter file simply type:

.. code-block:: none

    ./ria -pf myparameters.pf


The parameter file should be a text file with one option written on each line. For example, to perform the analysis above the file `myparameters.pf` would be as follows:

.. code-block:: none

    -plink /path/to/plink/plink
    -king /path/to/king/king mydata.bed
    -window-size-cm 15
    -step-size 50
    -i mydata.bed
    -o myResults.dat


It is also possible to add comments to the file provided that the "-" character is not used, and to comment out any options by placing another character in front of any "-". For example, the above parameter file could be edited as follows: *


.. code-block:: none

    # Command used to run PLINK
    -plink /path/to/plink/plink

    # Command used to run KING
    -king /path/to/king/king

    # SNP window size
    -window-size-cm 15

    # Number of SNPs to move to the next SNP window
    -step-size 50

    # The all important data to analyse
    -i mydata.bed

    # My lovely analysis results
    -o myResults.dat

    # I might try this later
    #-i myOtherData.bed


.. _temporary:

Temporary Files
---------------

RIA requires the use of two other programs, PLINK and KING, and in order to use these programs several temporary files must be created. All of these files begin with "tempRIA" and, before the file extension, end with the process ID of the current RIA job. They will typically appear similar to the following:

.. code-block:: none

    tempRIA-priors1-57703.prune.in
    tempRIA-priors2-57703.bed
    tempRIA-priors2-57703.bim
    tempRIA-priors2-57703.fam
    tempRIA-priors2-57703.kin
    tempRIA-post-57703.fam
    tempRIA-posterior-57703.bed
    tempRIA-posterior-57703.bim
    tempRIA-posterior-57703.fam
    tempRIA-posterior-57703.kin



The use of the process ID ensures that several RIA jobs may be ran in the same location without any issues of interference. All of these temporary files should be deleted by RIA when they are no longer needed, however if RIA is forced to unexpectedly stop for some reason then these files may still be lying around. In which case, they should be carefully deleted if necessary.

Using Prior and Posterior Files
-------------------------------

As an alternative to calculating the priors and posteriors they may be supplied to RIA if they have been calculated and saved in advance. The following bash shell script has example commands on how to do this. 

.. code-block:: none

    #!/bin/bash

    # There are lots of posterior files, so save them in their own directory
    mkdir posteriors-data3

    # Run on example data and save prior
    ../ria -king /path/to/king/king -plink /path/to/plink/plink -i ../demo-data/exampleRIAData.bed -prior-only -o-prior prior3.dat -log results-prior3.log

    # Run on example data with calculated prior and save posteriors
    ../ria -king /path/to/king/king -plink /path/to/plink/plink -i ../demo-data/exampleRIAData.bed -i-prior prior3.dat -o-posteriors-prefix posteriors-data3/post -o resultsRIAExampleData3-1.dat -log results3-1.log

    # Run on example data using calculated prior and posteriors
    ../ria -king /path/to/king/king -plink /path/to/plink/plink -i ../demo-data/exampleRIAData.bed -i-prior prior3.dat -i-posteriors-prefix posteriors-data3/post -o resultsRIAExampleData3-2.dat -log results3-2.log

    # Run on example data using calculated prior and posteriors using a window
    ../ria -king /path/to/king/king -plink /path/to/plink/plink -i ../demo-data/exampleRIAData.bed -i-prior prior3.dat -i-posteriors-prefix posteriors-data3/post -o resultsRIAExampleData3-3.dat -log results3-3.log -posterior-start-window 10 -posterior-end-window 20

    # Run on example data using calculated prior and posteriors, no .bed file and with end posteriors window
    ../ria -king /path/to/king/king -plink /path/to/plink/plink -i-prior prior3.dat -i-posteriors-prefix posteriors-data3/post -o resultsRIAExampleData3-4.dat -log results3-4.log -posterior-end-window 557

    # Run on example data using calculated prior and posteriors, no .bed file and with start and end posteriors window
    ../ria -king /path/to/king/king -plink /path/to/plink/plink -i-prior prior3.dat -i-posteriors-prefix posteriors-data3/post -o resultsRIAExampleData3-5.dat -log results3-5.log -posterior-start-window 10 -posterior-end-window 20


**Note:** If the priors and posteriors are both given with no .bed file then the last posterior window number must be given using option "-posterior-end-window".


