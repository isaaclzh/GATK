import os
import re
import sys
import pandas as pd

files = [f for f in os.listdir("fastq") if not f.startswith('.')]
df = pd.DataFrame(files)

# Generating samples.tsv file
df['Sample'] = df.iloc[:,0].str.extract('(.+?)_\d')
df['fq'] = df.iloc[:,0].str.extract('.+_(\d)')
df = df.pivot(index = 'Sample', columns='fq', values = 0).reset_index().rename_axis(columns=None)
os.system('mkdir config/misc')
os.system('mkdir config/TMP')
df.to_csv('config/misc/samples.tsv', sep="\t", index=False)

# Generating cohort.txt file
sys.stdout = open('config/misc/cohort.txt','wt')
for i in range(len(df['Sample'])):
    print("sample " + str(i+1) + "\tresults/" + str(df['Sample'][i]) + "/gvcf/" + str(df['Sample'][i]) + ".g.vcf.gz")
sys.stdout.close()
