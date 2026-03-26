import pandas as pd
from pathlib import Path

data_dir = Path(r'c:\Users\hamna pc\OneDrive\Desktop\GraphComm_Lite_Project\data\prostate\raw')
count_file = list(data_dir.glob('*.txt'))[0]

# Count lines
with open(count_file) as f:
    lines = sum(1 for _ in f)

# Read just header and first few rows to understand structure
df_small = pd.read_csv(count_file, sep='\t', nrows=10, index_col=0, low_memory=False)
print(f'Total lines in file: {lines}')
print(f'Sample shape: {df_small.shape}')
print(f'First 5 columns: {df_small.columns[:5].tolist()}...')
print(f'First 5 index values: {df_small.index.tolist()}')
print(f'\nFile structure: {df_small.shape[0]} genes x {df_small.shape[1]} cells (in sample)')
print(f'Estimated full file: ~{lines} genes')
