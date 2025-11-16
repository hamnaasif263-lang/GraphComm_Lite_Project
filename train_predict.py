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
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(cell_embeddings, drug_labels)
    print("Random Forest trained on cell embeddings for drug prediction.")
    return rf
