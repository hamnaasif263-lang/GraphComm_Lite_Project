import torch
import torch.nn as nn
import torch.nn.functional as F

class GraphCommLite(nn.Module):
    def __init__(self, in_feats, hidden, out_feats):
        super(GraphCommLite, self).__init__()
        self.fc1 = nn.Linear(in_feats, hidden)
        self.fc2 = nn.Linear(hidden, out_feats)

    def forward(self, x, adj):
        x = torch.matmul(adj, x)
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x
