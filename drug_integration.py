# drug_integration.py
# --------------------------------------------------------------
# Retrieve FDA-approved drugs for a given target gene
# --------------------------------------------------------------

import requests
import pandas as pd
import time
from typing import Optional


def _safe_get_json(resp: requests.Response) -> Optional[dict]:
    """Try to parse JSON; return None if payload isn't JSON.
    Provides a guard against HTML or empty responses from the API.
    """
    try:
        return resp.json()
    except ValueError:
        # not JSON
        return None


def query_fda_approved_drugs(target_gene: str, max_retries: int = 2) -> pd.DataFrame:
    """
    Query DGIdb API to get drugs interacting with a target gene.
    Returns a DataFrame of results. On API failure, attempts to load a local
    fallback CSV at `data/drug_db.csv` if present; otherwise returns empty DataFrame.
    """
    url = (
        f"https://dgidb.org/api/v2/interactions.json?genes={target_gene}"
        f"&interaction_sources=DrugBank"
    )

    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            data = _safe_get_json(resp)
            if data is None:
                raise RuntimeError("Non-JSON response from DGIdb")
            # success
            break
        except Exception as e:
            last_err = e
            wait = 1.0 * attempt
            print(f"Warning: attempt {attempt} failed for {target_gene}: {e}. Retrying in {wait}s...")
            time.sleep(wait)
    else:
        print(f"Error contacting DGIdb for {target_gene}: {last_err}")
        # fallback to local CSV if available
        try:
            local_path = "data/drug_db.csv"
            df_local = pd.read_csv(local_path)
            print(f"Loaded local drug DB from {local_path} as fallback.")
            return df_local[df_local['gene'].str.upper() == target_gene.upper()].reset_index(drop=True)
        except Exception:
            return pd.DataFrame()

    records = []
    for term in data.get('matchedTerms', []):
        for inter in term.get('interactions', []):
            records.append({
                "gene": target_gene,
                "drug_name": inter.get('drugName'),
                "interaction_source": inter.get('interactionClaimSource'),
                "is_approved": inter.get('isApproved', 'Unknown'),
                "interaction_types": inter.get('interactionTypes')
            })
    return pd.DataFrame(records)

def load_local_drug_list(csv_path: str) -> pd.DataFrame:
    """
    Load local CSV of FDA-approved drugs (for offline use).
    CSV columns: gene, drug_name, approval_status
    """
    return pd.read_csv(csv_path)

def map_cell_targets_to_drugs(target_genes: list, drug_df: pd.DataFrame) -> pd.DataFrame:
    """Filter local drug table by gene targets."""
    return drug_df[drug_df['gene'].isin(target_genes)].reset_index(drop=True)

# Demo test
if __name__ == "__main__":
    for g in ["IGF1R", "PIK3CA", "AKT1"]:
        df = query_fda_approved_drugs(g)
        print(f"\nDrugs for {g}:")
        print(df[["drug_name", "interaction_source"]].head())
