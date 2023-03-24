###############################

## BWA-GATK pipeline

This script is developed by the Lee Kong Chian School of Medicine (LKCMedicine) under Assoc. Prof Joanne Ngeow's lab, operating in the server in the Centre of Bioinformatics. It utilizes the **BWA-GATK** pipeline which converts paired short reads fastq to variant annotated files used for germline variant studies.

This pipeline uses the following libraries:
- bcftools 1.16
- bwa 0.7.17
- gatk 4.3.0
- htslib 1.16
- python 3.10.7
- samtools 1.16.1
- snakemake 7.18.2
- vep 107

This pipeline was benchmarked on Ashkenazim Son HG002 (https://www.nist.gov/programs-projects/genome-bottle) as the truth set, yielding 99.5% within confidence regions versus 91.9% outside confidence regions for single nucleotide polymorphisms; and 98.4% versus 72.0% for insertion-deletion mutations.

### Below are the instructions to using this script:

1. Before running, make sure that all reference libraries and additional files are already generated  
```source config/reference.sh```

2. Load the libraries  
```source config/lib.sh```

3. Upload your paired fastq inside the fastq folder  

* {file_name}_1{file_extension} for forward reads
* {file_name}_2{file_extension} for reverse reads

4. Change the number of chromosomes you wish to subsample inside the config.yml file (default: chr1-22, X, Y).   
```nano config/config.yml```

> Snakemake runs on 12 cores per sample (default)  
> Change the extension of fastq files inside the config.yml file (optional)  
> Change the path of the reference folder to the location of your reference folder (optional)  

5. Do a dry run  
```snakemake -c 24 --config bed=yes -np all```  

> Running snakemake on more than 4 samples at once would drastically reduce its performance  
> If running on 3 samples, change the -c 24 to -c 36  
> If you wish to run it without a bed file, set --config bed=no (applicable to steps 5 - 7)

6. Run snakemake by removing the dry run (-n) flag  
```snakemake -c 24 --config bed=yes -p all```

7. After running snakemake, you can generate an interactive HTML report for inspection of results inside the browser  
```snakemake --config bed=yes --report report.zip```
<br/><br/>
############################### 
<br/><br/>
Author: Lin Zhenhan Isaac  
Please email to zhenhanisaac.lin@ntu.edu.sg for any enquires
