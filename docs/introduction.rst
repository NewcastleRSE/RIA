.. _introduction:

Introduction
============

The program RIA is a C++ implementation of the method described in :cite:`nat:15` and uses calls to programs `PLINK <https://www.cog-genomics.org/plink/>`_ and `KING <https://www.kingrelatedness.com/history.shtml>`_ version 2.2.9.

.. _overview:


Overview
--------

The program RIA performs the following steps calling external PLINK and KING where required:


1. The first stage of RIA is the estimation of the IBD sharing probabilities between affected family members using a selection of SNPs across all available SNPs in the data, these probabilities form the *priors*. The program PLINK :cite:`purcell:etal:07` is used to prune the SNPs to give a representative selection of SNPs using only founders. Thus, it should be noted that founders are required in the SNP data. The cases are then used to estimate the IBD sharing probabilities using KING (version 2.2.9) :cite:`manichaikul:etal:10`. See :ref:`ibd-calc` on how the IBD sharing probabilities are estimated.

2. The next stage is to step across the genome, set to a default of 50 SNPs per step, and form a SNP window of around 500 to 2000 SNPs and use these SNPs to estimate the IBD sharing probabilities between affected family members using KING (version 2.2.9). These probabilities give the *posteriors*. See :ref:`ibd-calc` on how the IBD sharing probabilities are estimated.

3. The next step is to perform a non-parametric linkage analysis comparing the prior and posterior IBD sharing probabilities using the method as described in :cite:`cordell:etal:00` with minor modifications as described in :cite:`nat:15`. This produces a LOD score for each analysed SNP region together with parameter estimates for the (scaled) additive and dominance variances.


For full details of the methodology of RIA please see the accompanying paper :cite:`nat:15`.

.. _ibd_calc:

Estimation of IBD sharing probabilities
---------------------------------------


`KING <https://www.kingrelatedness.com/history.shtml>`_ version 2.2.9 is used with option `--homog` to estimate the IBD sharing probabilities. Let *IBD*2, *IBD*1 and *IBD*0 be the probabilities that two individuals share 2, 1 or 0 alleles IBD respectively and *K* the kinship coefficient, then *IBD*2, *IBD*1 and *IBD*0 are estimated for each SNP as follows:



1. KING version 2.2.9 is ran with the `--homog` option which returns estimates of *K* and *IBD*0.

2. *K* and *IBD*0 are truncated to plausible values between 0 to 0.5 and between 0 to 1 respectively if necessary.

3. *IBD*2 is estimated as 4*K* - (1 - *IBD*0).

4. *IBD*2 is truncated to a plausible value between 0 to 1 if necessary.

5. *IBD*1 is estimated as 1 - *IBD*0 - *IBD*2.

6. *IBD*1 is truncated to a plausible value between 0 to 1 if necessary.

7. For the posterior probabilities only: if any of the priors for *IBD*2, *IBD*1 or *IBD*0 are equal to zero then the corresponding posterior values are also set to zero. If this results in all IBD probabilities being equal to zero then an error is reported.

8. *IBD*2, *IBD*1 and *IBD*0 are scaled so that *IBD*2 + *IBD*1 + *IBD*0 = 1 if necessary.
