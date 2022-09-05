#! /bin/bash 

# Go to folder where bash script is. Coherence with relative paths
cd $(dirname $0)

# Filter VCF for our accessions
vcftools \
--vcf ./source_SV/CNVnator_Q10_goodRD_noCN1-3.vcf \
--keep ./utility/accessions.txt \
--recode \
--out ./output_CNV/CNVnator_Q10_goodRD_noCN1-3_accs

# With filtered VCF, compute missing rate per marker
vcftools \
--vcf ./output_CNV/CNVnator_Q10_goodRD_noCN1-3_accs.recode.vcf \
--missing-site \
--out ./output_CNV/CNVnator_Q10_goodRD_noCN1-3_accs.recode

# With filtered VCF, compute missing rate per accession
vcftools \
--vcf ./output_CNV/CNVnator_Q10_goodRD_noCN1-3_accs.recode.vcf \
--missing-indv \
--out ./output_CNV/CNVnator_Q10_goodRD_noCN1-3_accs.recode
