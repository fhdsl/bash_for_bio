import sys #<1>
import pandas as pd
file = sys.argv[1] #<2> 
df = pd.read_csv(file)
print(df.tumor_stage.nunique())