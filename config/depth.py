import os
import pandas as pd

sampl = []
depth = []
    
files = os.listdir('results')

for sample in files:
    seq_depth = pd.read_csv('results/' + sample + '/qc/coverage.txt', sep="\t", comment='#', engine='python', skipfooter=253)['MEAN_COVERAGE'].values[0]
    sampl.append(sample)
    depth.append(seq_depth)
    
df = pd.DataFrame({'Sample': sampl, 'Depth': depth})
df.to_csv('config/misc/depth.txt', sep="\t", index=False)
