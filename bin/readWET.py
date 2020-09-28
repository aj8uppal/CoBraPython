import pandas as pd
import sys, os, re

l0ed = []
l1ed = []
fld = []
l0wet = []
l1wet = []

for file in os.listdir('./output/'):

    match = re.search('l0ed_(\w+.\w+)_l1ed_(\w+.\w+)_fld_(\w+.\w+)_highpower_staveL0_A.csv', file)
    if match is not None:
        df0 = pd.read_csv('./output/{}'.format(file))
        df1 = pd.read_csv('./output/{}'.format(file).replace('staveL0', 'staveL1'))
        l0ed.append(float(match.group(1)))
        l1ed.append(float(match.group(2)))
        fld.append(float(match.group(3)))
        l0wet.append(df0.iloc[df0['heatFlux'].gt(10000).idxmax()]['Temperature'])
        l1wet.append(df1.iloc[df1['heatFlux'].gt(10000).idxmax()]['Temperature'])
        del df0
        del df1

df = pd.DataFrame(list(zip(l0ed, l1ed, fld, l0wet, l1wet)), columns=['l0ed', 'l1ed', 'fld', 'l0wet', 'l1wet'])
df.to_csv('./highpower_ed_fld.csv')
