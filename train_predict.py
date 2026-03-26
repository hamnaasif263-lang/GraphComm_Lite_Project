import torch
import torch.nn as nn
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from graph_model import GraphCommLite

def train_graph_model(features, adj, labels, epochs=50):
    device = torch.device("cpu")
    torch.set_num_threads(4)

    in_feats = features.shape[1]
    model = GraphCommLite(in_feats=in_feats, hidden=64, out_feats=2).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    loss_fn = nn.CrossEntropyLoss()

    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(features, adj)
        loss = loss_fn(outputs, labels)
        loss.backward()
        optimizer.step()
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{epochs} - Loss: {loss.item():.4f}")

    return model

def drug_response_prediction(cell_embeddings, drug_labels):
    """Train Random Forest on cell embeddings.
    
    Validates that labels have variance and that predictions are meaningful.
    """
    # Check label distribution
    unique_labels = np.unique(drug_labels)
    if len(unique_labels) < 2:
        print("⚠️ WARNING: Drug labels are not binary! All cells belong to single class.")
        print(f"   Unique labels: {unique_labels}")
        print("   This will result in trivial predictions (all 0 or all 1).")
    
    rf = RandomForestClassifier(n_estimators=100, random_state=42, min_samples_split=2, min_samples_leaf=1)
    rf.fit(cell_embeddings, drug_labels)
    
    # Check predictions for variance
    train_preds = rf.predict(cell_embeddings)
    unique_preds = np.unique(train_preds)
    if len(unique_preds) < 2:
        print("⚠️ WARNING: Random Forest is predicting only one class!")
        print(f"   All predictions are: {unique_preds[0]}")
        print("   This suggests label imbalance or insufficient signal in embeddings.")
    
    print("✅ Random Forest trained on cell embeddings for drug prediction.")
    return rf
