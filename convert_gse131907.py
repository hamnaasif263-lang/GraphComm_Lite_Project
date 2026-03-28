import gzip
import io
import os
import sys
from urllib.request import urlopen

import pandas as pd


URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz"
OUT_GZ = "data/GSE131907_series_matrix.txt.gz"
OUT_TXT = "data/GSE131907_series_matrix.txt"
OUT_CSV = "data/scRNA_data.csv"


def download_file(url, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    print(f"Downloading {url} -> {out_path}")
    try:
        with urlopen(url) as r, open(out_path, 'wb') as f:
            f.write(r.read())
        print("Download complete")
    except Exception as e:
        print("Download failed:", e)
        raise


def extract_matrix_block(gz_path):
    print(f"Reading gz file {gz_path}")
    with gzip.open(gz_path, 'rt', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            start_idx = i
        if line.startswith('!series_matrix_table_end'):
            end_idx = i
            break

    if start_idx is None or end_idx is None:
        raise RuntimeError('Could not find matrix begin/end markers in series matrix')

    # matrix lines are between start_idx+1 and end_idx
    matrix_lines = lines[start_idx + 1:end_idx]
    text = ''.join(matrix_lines)
    return text


def convert_and_save(matrix_text, out_csv):
    print("Parsing matrix into DataFrame")
    df = pd.read_csv(io.StringIO(matrix_text), sep='\t', header=0, index_col=0)
    print(f"Original matrix shape (genes x samples): {df.shape}")
    df_t = df.T
    print(f"Transposed matrix shape (samples x genes): {df_t.shape}")
    df_t.to_csv(out_csv)
    print(f"Wrote CSV to {out_csv}")


def main():
    try:
        download_file(URL, OUT_GZ)
    except Exception:
        print("Failed to download. If you already have the file, ensure it's at:", OUT_GZ)
        # continue if file exists
        if not os.path.exists(OUT_GZ):
            sys.exit(1)

    try:
        matrix_text = extract_matrix_block(OUT_GZ)
    except Exception as e:
        print("Failed to extract matrix block:", e)
        sys.exit(1)

    convert_and_save(matrix_text, OUT_CSV)


if __name__ == '__main__':
    main()
