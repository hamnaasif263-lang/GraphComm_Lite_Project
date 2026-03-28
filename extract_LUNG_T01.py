# extract_LUNG_T01.py
# Streams the GEO series_matrix.gz and extracts only columns (samples) belonging to LUNG_T01
# Produces data/scRNA_LUNG_T01_full.csv (cells x genes)

import os, gzip, csv, sys
import pandas as pd

ANNOT = "data/GSE131907_Lung_Cancer_cell_annotation.txt"
SERIES = "data/GSE131907_series_matrix.txt.gz"
OUT_TEMP = "data/_tmp_genes_selected.csv"   # genes as rows, selected sample columns
OUT_FINAL = "data/scRNA_LUNG_T01_full.csv"  # final: rows=cells, cols=genes

# ------------- helper: find barcodes for LUNG_T01 in annotation file -------------
def find_barcodes_for_patient(annotation_path, patient_tag="LUNG_T01"):
    df = pd.read_csv(annotation_path, sep="\t", comment='#', low_memory=False)
    # Try to find column that contains the patient id
    found_col = None
    for c in df.columns:
        try:
            if df[c].astype(str).str.contains(patient_tag).any():
                found_col = c
                break
        except Exception:
            continue
    if found_col is None:
        print("Could not find patient tag column automatically. Inspect annotation file and set correct column name.")
        print("Columns:", df.columns.tolist())
        raise SystemExit(1)
    # Now find rows matching the patient tag and get the cell IDs (barcodes)
    # Try common barcode-like columns: 'cell', 'barcode', 'cell_barcode', 'sample'
    barcode_col = None
    # heuristics to choose barcode column
    for c in df.columns:
        low = c.lower()
        if any(x in low for x in ["cell", "barcode", "cell_barcode", "cellid", "barcode_id", "barcodes"]):
            barcode_col = c
            break
    if barcode_col is None:
        # fallback to first column
        barcode_col = df.columns[0]
    barcodes = df.loc[df[found_col].astype(str).str.contains(patient_tag), barcode_col].astype(str).tolist()
    # Clean possible spaces
    barcodes = [b.strip() for b in barcodes if isinstance(b, str) and len(b.strip())>0]
    print(f"Found {len(barcodes)} barcodes for {patient_tag} (annotation col: '{found_col}', barcode col: '{barcode_col}')")
    return barcodes

# ------------- parse series_matrix header to map sample columns to barcodes -------------
def extract_selected_columns_from_series(series_gz, selected_barcodes, out_temp):
    # Open compressed series matrix
    with gzip.open(series_gz, "rt") as fh:
        reader = csv.reader(fh, delimiter='\t')
        # find header row (the one that starts with !Sample_title or sample IDs)
        header = None
        # We read lines until we encounter a header line (the first line that does not start with '!').
        # The standard GEO series_matrix has a header line with sample IDs as column names.
        for line in reader:
            if len(line)==0:
                continue
            if line[0].startswith("!"):
                continue
            # this is header (first non '!' line)
            header = line
            break
        if header is None:
            raise RuntimeError("Could not find header row in series matrix.")
        # header is something like: ID_REF    GSMxxxx   GSMyyyy ...
        # we assume sample names in header[1:]
        sample_names = header[1:]
        # Some GEO matrices have sample names like GSM3827114_LUNG_T01_... or only GSM ids
        # We will match selected_barcodes against sample_names
        # But often selected_barcodes are cell barcodes (not GSM ids) -> need different approach:
        # Many GEO series_matrix are gene x sample (sample = GSM), not single-cell barcodes.
        # If this series matrix is not per-cell but per-sample, we cannot extract single-cell columns here.
        # We'll attempt to match sample_names to barcodes; if no match, we will stop.
        # Check matches
        matches = []
        for i, s in enumerate(sample_names, start=1):
            for b in selected_barcodes:
                if b in s or s in b:
                    matches.append(i)
                    break
        if len(matches)==0:
            # No direct matches: warn and abort (this series file may be per-sample not per-cell)
            print("No sample header matched the barcodes. This series_matrix may be aggregated per sample (GSM) rather than per-cell.")
            print("In that case you need the raw filtered_feature_bc_matrix (MTX / h5) per-sample, or process the 10x h5.")
            raise RuntimeError("No header matches")
        # else, we will stream the file and write out columns: gene + selected sample columns
        print(f"Will extract {len(matches)} matching columns (first few indexes): {matches[:10]}")
        # Write temporary CSV: header gene + selected sample names
        with open(out_temp, "w", newline='') as outfh:
            writer = csv.writer(outfh)
            out_header = ["Gene"] + [sample_names[i-1] for i in matches]
            writer.writerow(out_header)
            # Process remaining lines (each line has gene + expression values)
            # Note: must continue reading from current position of 'reader' (we consumed header)
            for row in reader:
                if len(row)==0:
                    continue
                gene = row[0]
                # pick values at indexes in matches
                selected_vals = [row[i] if i < len(row) else "0" for i in matches]
                writer.writerow([gene] + selected_vals)
    print("Temporary gene x selected-samples CSV written:", out_temp)
    return out_temp

# ------------- transpose genes x samples -> samples x genes -------------
def transpose_to_cells_by_genes(tmp_csv, out_csv):
    print("Loading temporary CSV into pandas (this is much smaller).")
    df = pd.read_csv(tmp_csv, index_col=0)
    # df: genes x selected_samples -> transpose
    df_t = df.transpose()
    # set simple cell IDs (if header sample names are GSMs, we keep them)
    df_t.to_csv(out_csv, index=True)
    print("Final CSV saved:", out_csv)
    return out_csv

if __name__ == "__main__":
    # Accept patient tag as optional CLI argument, default kept for backwards compatibility
    if len(sys.argv) > 1:
        patient = sys.argv[1]
    else:
        patient = "LUNG_T01"   # default (may not exist in this annotation)
    if not os.path.exists(ANNOT):
        print("Annotation file not found at", ANNOT)
        sys.exit(1)
    if not os.path.exists(SERIES):
        print("Series matrix not found at", SERIES)
        print("Please download the series matrix from GEO and place it at that path.")
        sys.exit(1)

    # Step 1: find barcodes for this patient
    barcodes = find_barcodes_for_patient(ANNOT, patient_tag=patient)

    if len(barcodes)==0:
        print("No barcodes found for", patient)
        sys.exit(1)

    # Step 2: try extract selected columns from series matrix (stream)
    try:
        tmp = extract_selected_columns_from_series(SERIES, barcodes, OUT_TEMP)
    except Exception as e:
        print("Extraction failed:", e)
        print("This series_matrix format may not contain per-cell columns; you may need the per-sample 10x h5 files.")
        sys.exit(1)

    # Step 3: transpose and save final CSV
    transpose_to_cells_by_genes(tmp, OUT_FINAL)
    print("Done. Final file at:", OUT_FINAL)
