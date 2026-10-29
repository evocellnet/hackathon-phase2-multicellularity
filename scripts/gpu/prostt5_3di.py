#!/usr/bin/env python
"""Predict Foldseek 3Di strings from amino-acid sequences with ProstT5 (encoder + CNN head) on GPU.

Re-implements the reference script `scripts/predict_3Di_encoderOnly.py` of the ProstT5 repository (Heinzinger et al.
2024, NAR Genom Bioinform 6:lqae150; github.com/mheinzinger/ProstT5; copy kept at data/raw/prostt5/
predict_3Di_encoderOnly.reference.py): ProstT5 T5 encoder in half precision, prefix "<AA2fold>", residues
space-separated, non-standard residues -> X, two-layer CNN (cnn_chkpnt/model.pt) predicting one of 20 3Di states per
residue, argmax. This is the same model family Foldseek uses for `createdb --prostt5-model` (there run on CPU here, because
the Foldseek GPU build fell back to CPU on this node).

Batches are formed from length-sorted sequences under a residue budget. Sequences longer than --window are split into
overlapping windows (--window residues, step --window minus --overlap); each residue takes its state from the window in which it is
farthest from a window edge, then the protein is written once complete. Sequences longer than --max-len are skipped (listed in
<out>.skipped). Resumable: IDs already present in <out> are not recomputed. Several input FASTAs may be given.
"""
import argparse
import os
import sys
import time

import torch
from torch import nn
from transformers import T5EncoderModel, T5Tokenizer

STATES = "ACDEFGHIKLMNPQRSTVWY"
STD = set(STATES)


class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.classifier = nn.Sequential(
            nn.Conv2d(1024, 32, kernel_size=(7, 1), padding=(3, 0)), nn.ReLU(), nn.Dropout(0.0),
            nn.Conv2d(32, 20, kernel_size=(7, 1), padding=(3, 0)))

    def forward(self, x):  # x: B x L x 1024 -> B x 20 x L
        return self.classifier(x.permute(0, 2, 1).unsqueeze(-1)).squeeze(-1)


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", help="output 3Di FASTA (appended; resumable)")
    ap.add_argument("fasta", nargs="+")
    ap.add_argument("--cnn", default="data/raw/prostt5/cnn_chkpnt_model.pt")
    ap.add_argument("--max-len", type=int, default=100000)
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--overlap", type=int, default=250)
    ap.add_argument("--residues-per-batch", type=int, default=40000)
    ap.add_argument("--max-batch", type=int, default=512)
    ap.add_argument("--model", default="Rostlab/ProstT5")
    args = ap.parse_args()

    dev = torch.device("cuda:0")
    tok = T5Tokenizer.from_pretrained(args.model, do_lower_case=False)
    enc_model = T5EncoderModel.from_pretrained(args.model).to(dev).half().eval()
    cnn = CNN()
    cnn.load_state_dict(torch.load(args.cnn, map_location="cpu")["state_dict"])
    cnn = cnn.to(dev).half().eval()

    seen = set()
    if os.path.exists(args.out):
        seen = {l[1:].strip() for l in open(args.out) if l.startswith(">")}
    seqs, skipped = [], []
    for fa in args.fasta:
        for name, s in read_fasta(fa):
            if name in seen:
                continue
            s = "".join(c if c in STD else "X" for c in s.upper().rstrip("*"))
            (skipped if len(s) > args.max_len else seqs).append((name, s))
    full_len, pending = {}, {}
    units = []  # (protein, start, subsequence)
    step = args.window - args.overlap
    for name, s in seqs:
        full_len[name] = len(s)
        if len(s) <= args.window:
            units.append((name, 0, s))
        else:
            starts = list(range(0, len(s) - args.window, step)) + [len(s) - args.window]
            for st in sorted(set(starts)):
                units.append((name, st, s[st:st + args.window]))
    n_prot = len(seqs)
    seqs = sorted(units, key=lambda x: -len(x[2]))
    with open(args.out + ".skipped", "a") as fh:
        for name, s in skipped:
            fh.write(f"{name}\t{len(s)}\n")

    t0, n, nres, i = time.time(), 0, 0, 0
    last = t0
    with open(args.out, "a") as out:
        while i < len(seqs):
            L = len(seqs[i][2])
            bs = max(1, min(args.max_batch, args.residues_per_batch // max(L, 1)))
            batch = seqs[i:i + bs]
            i += bs
            texts = ["<AA2fold> " + " ".join(s) for _, _, s in batch]
            te = tok(texts, add_special_tokens=True, padding="longest", return_tensors="pt").to(dev)
            with torch.no_grad():
                emb = enc_model(te.input_ids, attention_mask=te.attention_mask).last_hidden_state
                mask = te.attention_mask.clone()
                for k, (_, _, s) in enumerate(batch):
                    mask[k, len(s) + 1] = 0  # drop the end-of-sequence token, as in the reference
                emb = (emb * mask.unsqueeze(-1))[:, 1:]  # drop the prefix token
                pred = cnn(emb).argmax(dim=1).cpu().numpy()
            for k, (name, st, s) in enumerate(batch):
                states = pred[k, :len(s)]
                if st == 0 and len(s) == full_len[name]:
                    out.write(f">{name}\n{''.join(STATES[j] for j in states)}\n")
                    n += 1
                    continue
                buf = pending.setdefault(name, [[None] * full_len[name], [-1] * full_len[name], 0])
                for off, j in enumerate(states):
                    pos = st + off
                    margin = min(off, len(s) - 1 - off)
                    if margin > buf[1][pos]:
                        buf[0][pos], buf[1][pos] = STATES[j], margin
                buf[2] += 1
                if all(c is not None for c in buf[0]):
                    out.write(f">{name}\n{''.join(buf[0])}\n")
                    del pending[name]
                    n += 1
            nres += sum(len(s) for _, _, s in batch)
            if time.time() - last > 60:
                last = time.time()
                el = last - t0
                print(f"{n}/{n_prot} proteins {nres} res {el:.0f}s {nres/el:.0f} res/s", file=sys.stderr, flush=True)
    el = time.time() - t0
    if pending:
        print(f"ERROR: {len(pending)} windowed proteins incomplete", file=sys.stderr)
        sys.exit(1)
    print(f"DONE {n} proteins {nres} residues in {el:.0f}s ({nres/max(el,1e-9):.0f} res/s); skipped {len(skipped)} > {args.max_len} aa",
          file=sys.stderr)


if __name__ == "__main__":
    main()
