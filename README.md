## Welch-weighted Egger regression (WWER)
This github repo contains the code for the analyses described in
``Welch-weighted Egger regression reduces false positives in Mendelian 
randomization due to correlated pleiotropy''.

The primary method functions are `WWER::welch_test()` which is used to calculate
the weights, and `WWER::WWER()` which is used to estimate the effect.
For example, if I have a set of SNP effect esimates for LD-pruned SNPs
in two samples for two phenotypes:
```
# beta_exp_1, se_exp_1: exposure effects and SEs, sample 1.
# beta_exp_2, se_exp_2: exposure effects and SEs, sample 2.
# beta_out_1, se_out_1: outcome effects and SEs, sample 1.
# beta_out_2, se_out_2: outcome effects and SEs, sample 2.
weights = WWER::welch_test(beta_exp_1, beta_out_1, se_exp_1, se_out_1)$t
result = WWER::WWER(beta_exp_2, beta_out_2, se_exp_2, se_out_2, weights)
```

In practice, if you are applying this in the intended setting of exploratory
analysis of many phenotypes, the analysis code will be more complicated.

We will assume your data are in two objects `sumstats_1`
and `sumstats_2`. Each of these two objects should have two members
`$beta_hat` and `$se_hat` which are genotypes by phenotypes matrices
of effect estimates and standard errors for your dataset. In the file
`sumstats.R` we provide some helper functions for parsing your data into this
format from the Neale lab UKBB sumstats file format. See
`WWER::read_sumstats_neale()` as well as its use in
`WWER/Rmd/ukbb_analysis.Rmd`.

Once you have your data matrices, the first step is usually to choose 
instruments and calculate their weights using the first data matrix. This is done
using `selected_snps <- WWER::select_snps(sumstats_1)`.
The default settings correspond to the parameters used
in the analysis in the manuscript, however you can adjust the parameters to
your liking. For example, to use a more liberal SNP inclusion threshold you
can set `p_thresh=5e-06`, or to filter SNPs with similar scores more stringently
you can set e.g. `filter=3.2`. Please note this function does not take into
account LD. If your genotypes have already been pruned to be approximately
independent, you can use this directly. If they are not independent, this
function takes a `snps_to_use` argument which accepts a named list mapping
phenotype names to lists of SNPs that are independent and should be used for
that phenotype (for example coming from running PLINK clump on that sumstats
file.)

Once you have your instruments selected you can calculate a matrix of effects
using `WWER::fit_tce(sumstats_2, selected_snps)` which will return a matrix
of phenotype - phenotype effects as well as standard errors for all pairs
of phenotypes with valid instruments.

For a complete example of running this pipeline on the entire UKBB phenotype
set, see `WWER/Rmd/ukbb_preprocessing.Rmd` and `WWER/Rmd/ukbb_analysis.Rmd`.


### Frequenttly Asked Questions:
1. What is WWER?

    WWER is an approach to Mendelian randomization. It uses two datasets to
    attempt to downweight genetic instruments that may be acting through a
    third (unobserved) factor. This can reduce the false positives due to
    correlated pleiotropy.

2. What are the advantages of WWER over other methods for doing MR?

   There are two main advantages of WWER. The first is that it reduces false
   positives due to correlated pleiotropy. The second is that it is very fast.
   Compared to other methods that control false positives due to correlated
   pleiotropy, it is much faster. Compared to other fast methods, it is much
   better at controlling for correlated pleiotropy.

3. When should I use WWER instead of other methods?

   WWER is ideal for settings that require both speed and control of the FPR
   in the presence of correlated pleiotropy. It was created specifically for
   exploratory analysis in biobank-style data. However it may be helpful in
   any settings where you have more phenotype pairs than bandwidth to manually
   validate MR method assumptions for every pair, or if you don't have time
   to run slower mixture models on every pair.

4. What data do I need to use WWER?

    To estimate an effect for one pair of phenotypes, you need four sets of
    GWAS summary statistics: two sets with non-overlapping samples for each
    phenotype. In practice, this method is best suited for application to
    many pairs of phenotypes simultaneously. In that setting you need two
    non-overlapping sample sets for each phenotype.


For more information please see:
https://www.biorxiv.org/content/10.1101/2021.04.09.439229v1
