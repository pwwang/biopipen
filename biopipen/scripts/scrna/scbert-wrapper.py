#!/usr/bin/env python3
"""scBERT wrapper for CellTypeAnnotation — runs inference and saves per-cell predictions."""
from argparse import ArgumentParser
import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn


class Identity(nn.Module):
    """Classification head matching scBERT's predict.py Identity module."""

    def __init__(self, dropout=0.0, h_dim=100, out_dim=10, seq_len=16907):
        super().__init__()
        self.conv1 = nn.Conv2d(1, 1, (1, 200))
        self.act1 = nn.ReLU()
        self.fc1 = nn.Linear(in_features=seq_len, out_features=512)
        self.act2 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(512, h_dim)
        self.act3 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        self.fc3 = nn.Linear(h_dim, out_dim)

    def forward(self, x):
        x = x[:, None, :, :]
        x = self.conv1(x)
        x = self.act1(x)
        x = x.view(x.shape[0], -1)
        x = self.fc1(x)
        x = self.act2(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = self.act3(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        return x


def main():
    parser = ArgumentParser(description="Run scBERT prediction")
    parser.add_argument(
        "-i", "--input", required=True, help="Input H5AD file (AnnData)"
    )
    parser.add_argument(
        "-m", "--model", required=True,
        help="Fine-tuned model checkpoint (.pth)"
    )
    parser.add_argument(
        "-l", "--label-dict", required=True,
        help="Label dictionary pickle file"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file for predictions"
    )
    parser.add_argument(
        "-r", "--ref", required=True,
        help="Path to scBERT repo (containing performer_pytorch/)"
    )
    parser.add_argument(
        "--bin-num", type=int, default=5,
        help="Number of bins (default: 5)"
    )
    parser.add_argument(
        "--gene-num", type=int, default=16906,
        help="Number of genes expected by model (default: 16906)"
    )
    parser.add_argument(
        "--seed", type=int, default=2021, help="Random seed (default: 2021)"
    )
    parser.add_argument(
        "--pos-embed", action="store_true", default=True,
        help="Use Gene2vec positional encoding (default: True)"
    )
    parser.add_argument(
        "--no-pos-embed", action="store_true",
        help="Disable Gene2vec positional encoding"
    )
    parser.add_argument(
        "--novel-type", action="store_true",
        help="Enable novel cell type detection"
    )
    parser.add_argument(
        "--unassign-thres", type=float, default=0.5,
        help="Confidence threshold for unassigned cells (default: 0.5)"
    )

    args = parser.parse_args()

    # Resolve pos_embed
    pos_embed = not args.no_pos_embed if args.no_pos_embed else args.pos_embed

    # Add scBERT repo to path
    sys.path.insert(0, args.ref)
    from performer_pytorch import PerformerLM

    CLASS = args.bin_num + 2
    SEQ_LEN = args.gene_num + 1
    UNASSIGN_THRES = args.unassign_thres if args.novel_type else 0.0
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Set seeds
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Load data
    adata = sc.read_h5ad(args.input)
    barcodes = adata.obs_names.tolist()
    data = adata.X

    # Load label dict
    with open(args.label_dict, "rb") as f:
        label_dict = pickle.load(f)
    # Ensure we have index -> name mapping
    # label_dict from training maps name -> index, reverse it
    if isinstance(next(iter(label_dict)), str):
        label_dict = {v: k for k, v in label_dict.items()}
    num_types = len(label_dict)

    # Build model
    model = PerformerLM(
        num_tokens=CLASS,
        dim=200,
        depth=6,
        max_seq_len=SEQ_LEN,
        heads=10,
        local_attn_heads=0,
        g2v_position_emb=pos_embed,
    )
    model.to_out = Identity(
        dropout=0.0, h_dim=128, out_dim=num_types, seq_len=SEQ_LEN
    )

    # Load checkpoint
    ckpt = torch.load(args.model, map_location=device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(device)
    model.eval()

    # Freeze parameters
    for param in model.parameters():
        param.requires_grad = False

    # Predict cell by cell (matching predict.py logic)
    num_cells = data.shape[0]
    pred_list = ["Unassigned"] * num_cells

    with torch.no_grad():
        for index in range(num_cells):
            full_seq = data[index].toarray()[0]
            full_seq[full_seq > (CLASS - 2)] = CLASS - 2
            full_seq = torch.from_numpy(full_seq).long()
            full_seq = torch.cat(
                (full_seq, torch.tensor([0]))
            ).to(device)
            full_seq = full_seq.unsqueeze(0)

            logits = model(full_seq)
            pred_prob = nn.Softmax(dim=-1)(logits)
            pred_final = pred_prob.argmax(dim=-1).item()

            if pred_prob.max().item() < UNASSIGN_THRES:
                pred_list[index] = "Unassigned"
            else:
                pred_list[index] = label_dict[pred_final]

    # Save results
    results = pd.DataFrame({
        "barcode": barcodes,
        "scbert_celltype": pred_list,
    })
    results.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
