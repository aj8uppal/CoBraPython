import pandas as pd
import sys
df = pd.read_csv(sys.argv[1])
print(df.iloc[df['heatFlux'].gt(10000).idxmax()]['Temperature'])
