#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -m beas
#$ -pe threaded 4
#$ -l h_rt=04:00:00


# This script performs an approximate conditional stepwise association analysis to select independently associated SNPs for AMH - we used GCTA for this step
# The input file gcta.input.AMHGWAS.ma was created using the following code:
# cat ../Final_METALresults/AMH_METAANALYSIS_main_final1.TBL | awk '{print $1, $2, $3, $4, $8, $9, $10, $12}' > gcta.input.AMHGWAS.ma
# sed -i -e '1s/MarkerName/SNP/' -e '1s/Allele1/A1/' -e '1s/Allele2/A2/' -e '1s/Freq1/freq/' -e '1s/Effect/b/' -e '1s/StdErr/se/' -e '1s/P-value/p/' -e '1s/TotalSampleSize/N/' gcta.input.AMHGWAS.ma
#
# We performed these steps only for chromosome 2, 19 and 20 since our 4 significant AMH loci reside at those chromosomes



gcta=/path/to/GCTA/gcta_1.93.1beta/gcta64

for chr in 2 19 20

do

$gcta --bfile epicnl.ped.norelatedness.for.gcta_chr$chr \
	    --chr $chr \
	    --cojo-file gcta.input.AMHGWAS.ma \
          --cojo-slct \
	    --cojo-p 5e-8 \
	    --thread-num 4 \
	    --out /path/to/AMH_GWAS/MetaAnalysis/Output/Meta_analysis/Paper/Final_GCTAresults/cojo.gcta.chr$chr
     
done
