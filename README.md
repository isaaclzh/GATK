###############################

This script is developed by the Lee Kong Chian School of Medicine (LKCMedicine) under Assoc. Prof Joanne Ngeow's lab, operating in the server in the Centre of Bioinformatics (CBI). It utilizes the BWA-GATK workflow which converts paired short reads fastq files to variant annotation used for germline variant studies.

This pipeline uses the following libraries
- bcftools 1.16
- bwa 0.7.17
- htslib 1.16
- python 3.10.7
- samtools 1.16.1
- snakemake 7.18.2
- vep 107

This pipeline was benchmarked on Ashkenazim Son HG002 (https://www.nist.gov/programs-projects/genome-bottle) as the truth set, yielding 99.5% within confidence regions versus 91.9% outside confidence regions for single nucleotide polymorphisms; and 98.4% versus 72.0% for insertion-deletion mutations.

Please email to zhenhanisaac.lin@ntu.edu.sg for any enquires.

###############################
