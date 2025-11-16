import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score


def get_relapse_labels(df: pd.DataFrame, signature: list = None, label_columns=None):
    """Return a pandas Series of binary relapse labels.
    Priority:
      - If one of `label_columns` exists in df, use it (cast to int)
      - Else if `signature` provided, compute mean signature expression per sample and median-split
      - Else raise ValueError
    """
    if label_columns is None:
        label_columns = ["relapse", "relapse_status", "label"]

    for col in label_columns:
        if col in df.columns:
            return df[col].astype(int)

    if signature is None or len(signature) == 0:
        raise ValueError("No label column found and no signature provided for relapse score.")

    present = [g for g in signature if g in df.columns]
    if not present:
        raise ValueError("None of the signature genes were found in the expression dataframe.")

    score = df[present].mean(axis=1)
    thresh = score.median()
    labels = (score > thresh).astype(int)
    return labels


def evaluate_classifier(embeddings: np.ndarray, labels: np.ndarray, n_splits: int = 5):
    """Train and evaluate a RandomForest classifier on embeddings.
    Returns a dict with accuracy and ROC-AUC (where computable).
    Uses stratified CV when possible; falls back to a single train/test split for very small n.
    """
    n = embeddings.shape[0]
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    results = {"accuracy_mean": None, "accuracy_std": None, "roc_auc_mean": None}

    # If too few samples for stratified KFold, do a single split
    if n < 2 or len(np.unique(labels)) < 2:
        print("Not enough samples or label variety to evaluate classifier properly.")
        clf.fit(embeddings, labels)
        preds = clf.predict(embeddings)
        acc = accuracy_score(labels, preds)
        results["accuracy_mean"] = acc
        return results

    splits = min(n, n_splits)
    if splits >= 2 and n >= splits:
        try:
            skf = StratifiedKFold(n_splits=splits, shuffle=True, random_state=42)
            acc_scores = cross_val_score(clf, embeddings, labels, cv=skf, scoring="accuracy")
            results["accuracy_mean"] = float(acc_scores.mean())
            results["accuracy_std"] = float(acc_scores.std())
            # compute ROC AUC with CV manually if possible
            aucs = []
            for train_idx, test_idx in skf.split(embeddings, labels):
                clf.fit(embeddings[train_idx], labels[train_idx])
                probs = clf.predict_proba(embeddings[test_idx])[:, 1]
                try:
                    aucs.append(roc_auc_score(labels[test_idx], probs))
                except Exception:
                    pass
            if aucs:
                results["roc_auc_mean"] = float(np.mean(aucs))
        except Exception as e:
            # fallback to single split
            Xtr, Xte, ytr, yte = train_test_split(embeddings, labels, test_size=0.33, stratify=labels, random_state=42)
            clf.fit(Xtr, ytr)
            preds = clf.predict(Xte)
            probs = clf.predict_proba(Xte)[:, 1] if hasattr(clf, "predict_proba") else None
            results["accuracy_mean"] = float(accuracy_score(yte, preds))
            if probs is not None:
                try:
                    results["roc_auc_mean"] = float(roc_auc_score(yte, probs))
                except Exception:
                    pass

    return results


def top_upregulated_genes(df: pd.DataFrame, labels: np.ndarray, top_k: int = 5, candidate_genes=None):
    """Return top_k genes most upregulated in label==1 vs label==0.
    If candidate_genes provided, restrict to that list.
    """
    if candidate_genes is not None:
        genes = [g for g in candidate_genes if g in df.columns]
    else:
        genes = list(df.columns)

    mean1 = df.loc[labels == 1, genes].mean(axis=0)
    mean0 = df.loc[labels == 0, genes].mean(axis=0)
    fc = (mean1 + 1e-9) / (mean0 + 1e-9)
    fc = fc.sort_values(ascending=False)
    return list(fc.head(top_k).index)
