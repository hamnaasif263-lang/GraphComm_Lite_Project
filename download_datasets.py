import os
import GEOparse

# Root data directory
DATA_DIR = "data"

# GEO datasets for each cancer
datasets = {
    "breast_cancer": "GSE161529",
    "lung_cancer": "GSE131907",
    "prostate_cancer": "GSE141445"
}

for cancer_type, geo_id in datasets.items():

    save_path = os.path.join(DATA_DIR, cancer_type)

    # Create directory if not exists
    os.makedirs(save_path, exist_ok=True)

    print(f"\nDownloading {geo_id} into {save_path} ...")

    try:
        gse = GEOparse.get_GEO(
            geo=geo_id,
            destdir=save_path,
            annotate_gpl=True,
            silent=False
        )

        print(f"{geo_id} successfully downloaded into {save_path}")
    except Exception as e:
        print(f"Error downloading {geo_id}: {str(e)}")
        print(f"You may need to download {geo_id} manually from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}")

print("\nDataset download process completed.")
