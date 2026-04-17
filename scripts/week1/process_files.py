import sys #<1>
import pandas as pd
file = sys.argv[1] #<2> 
df = pd.read_csv(file)
summary = df.tumor_stage.value_counts()
summary.to_csv(file+".summary.csv")