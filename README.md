# GWAS meta-analysis of anti-Müllerian hormone (AMH) levels in 7,049 premenopausal women

This repository includes the code used to perform analyses for the following publication:

Verdiesen, Renée MG, et al. "Genome-wide association study meta-analysis identifies three novel loci for circulating anti-Müllerian hormone levels in women." Human Reproduction 37.5 (2022): 1069-1082. doi: 10.1093/humrep/deac028. 


## Part 1: Code

#### 1. 1_EasyQC_FilelevelQC.ecf
Script for performing file-level QC of all study-specific files, using EasyQC. This script is based on the example script at https://homepages.uni-regensburg.de/~wit59712/easyqc/1000g/fileqc_1000G.ecf. Information about EasyQC can be found at: https://homepages.uniregensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf.

1_EasyQC_FilelevelQC.ecf includes the following steps:
1. Removal of rows with missing data on alleles, pvalue, beta, se, and/or effect allele frequency
2. Removal of rows with unrealistic values (e.g. PVAL < 0 | > 1)
3. Reduction of number of significant digits
4. Exclusion of monomorphic SNPs
5. Harmonize alleles
6. Remove INDELS (insertions and deletions)
7. Create CHRPOS column (CHR:POS format)
8. Create Markername (CHR:POS:EFFECT_ALLELE:OTHER_ALLELE) column
9. Filter duplicate SNPs
10. Write final input files for meta level QC and METAL analyses


#### 2. 2_EasyQC_MetalevelQC.ecf
Script for performing meta-level QC of all study-specific files, using EasyQC. This script is based on the example script at https://homepages.uni-regensburg.de/~wit59712/easyqc/1000g/fileqc_1000G.ecf.

2_EasyQC_MetalevelQC.ecf includes the following steps:
1. Create SE-N plot
2. Create P-Z scatter plot
3. Alelle frequency check
4. Create QQ-plot 
5. Identification of population stratification


#### 3. 3_Metal_metaanalysis_phet.txt
Script for performing the main AMH GWAS meta-analysis in METAL (https://genome.sph.umich.edu/wiki/METAL_Documentation).


#### 4. 4_GCTA_conditionalanalysis.sh
Script for performing approximate conditional and joint association analysis (cojo-slct) using GCTA software (https://cnsgenomics.com/software/gcta/ ).


#### 5. 5_DEPICT_pathwayanalysis.cfg
Scipt for performing pathway analysis using DEPICT (https://github.com/perslab/depict). For this analysis we included all suggestive significant SNPs (p < 5 x 10-6), which were clumped at LD r2 < 0.1 and a physical distance of 500kb using PLINK v.1.9 as part of the DEPICT pipeline.


####6. 6_R_twosampleMR.R
Script for performing two-sample Mendelian Randomization analyses for breast cancer and polycystic ovary syndrome (PCOS). For these analyses we used summary-level from the breast cancer GWAS by K Michailidou et al. (2017) (doi: 10.1038/nature24284) and the PCOS GWAS by F Day et al. (2018) (doi: 10.1371/journal.pgen.1007813). Mendelian Randomization analyses were performed the R package TwoSampleMR (version 0.5.1) (https://mrcieu.github.io/TwoSampleMR/).



## Part 2: Summary-level GWAS results

The summary-level data from the AMH GWAS meta-analysis will be posted to the GWAS catalog (https://www.ebi.ac.uk/gwas/). The corresponding link will be included here once our data are available online.
