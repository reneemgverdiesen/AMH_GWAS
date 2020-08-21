# GWAS meta-analysis of anti-Mullerian hormone (AMH) levels in 7049 premenopausal women

This repository includes code used to perform analyses for the manuscript 'Genome-wide association study meta-analysis identifies three novel loci for circulating  anti-Müllerian hormone levels in women' (available on bioRxiv: xxx). In addition, the summary-level GWAS results from this project can be found in this respository.



## Part 1: Code

#### 1. 1_EasyQC_FilelevelQC.ecf
Script for performing file-level QC of all study-specific files, using EasyQC (referentie). This script is based on the example script at https://homepages.uni-regensburg.de/~wit59712/easyqc/1000g/fileqc_1000G.ecf . Information about EasyQC can be found at: https://homepages.uniregensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf .

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
Script for performing meta-level QC of all study-specific files, using EasyQC (referentie). This script is based on ..

2_EasyQC_MetalevelQC.ecf includes the following steps:
1. Create SE-N plot
2. Create P-Z scatter plot
3. Alelle frequency check
4. Create QQ-plot 
5. Identification of population stratification


#### 3. 3_Metal_metaanalysis_phet.txt
Script for performing the main AMH GWAS meta-analysis in METAL (xx).


#### 4. 4_GCTA_conditionalanalysis.sh
Script for performing approximate conditional and joint association analysis using GCTA (referentie).


#### 5. 5_DEPICT_pathwayanalysis.xx
Scipt for performing pathway analysis using DEPICT (ref). For this analysis we included all suggestive significant SNPs (p < 5 x 10-6), which were clumped at LD r2 < 0.1 and a physical distance of 500kb using PLINK v.1.9 as part of the DEPICT pipeline (ref).


#### 6. 6_R_MRbase.R



## Part 2: Summary-level GWAS results

