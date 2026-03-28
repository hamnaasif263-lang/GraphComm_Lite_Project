import gzip

path = 'data/GSE131907_series_matrix.txt.gz'
with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
    for i, line in enumerate(f):
        if i < 200:
            print(f"{i+1:04d}: {line.rstrip()}")
        else:
            break
