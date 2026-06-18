"""
G4 prediction on FASTA genome — run from terminal:

    python predict_g4.py \
        --fasta   /Users/dassagaripova/Downloads/Squid_export/GCF_001194135.2_ASM119413v2_genomic.fna \
        --tokenizer  /Users/dassagaripova/Downloads/Squid_export/6-new-12w-0 \
        --model      /Users/dassagaripova/Downloads/Squid_export/dnabert_mm_fold_0_kouzine_g4 \
        --out_dir    predictions_g4 \
        --device     cpu
"""

import argparse
import os
import gc
import pickle
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import torch
from torch.utils import data
from torch.utils.data import DataLoader
from tqdm import tqdm
from Bio import SeqIO
from transformers import BertConfig, BertForTokenClassification

# dna_tokenizer.py должен лежать рядом с этим скриптом
from dna_tokenizer import DNATokenizer, seq2kmer


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def split_seq(seq, length=512, pad=16):
    res = []
    n = len(seq)
    for st in range(0, n, length - pad):
        if st > 0 and st + pad >= n:
            break
        end = min(st + length, n)
        res.append(seq[st:end])
    return res


def stitch_np_seq(np_seqs, pad=16):
    total_length = sum(s.shape[-1] for s in np_seqs) - pad * (len(np_seqs) - 1)
    res = np.empty(total_length, dtype=np_seqs[0].dtype)
    pos = 0
    for i, seq in enumerate(np_seqs):
        seq_len = seq.shape[-1]
        if i > 0:
            pos -= pad
        res[pos:pos + seq_len] = seq[0, :]
        pos += seq_len
    return res


class PredDataset(data.Dataset):
    def __init__(self, sequence, tokenizer):
        self.pieces = split_seq(seq2kmer(sequence.upper(), 6).split(' '), length=512, pad=16)
        self.tokenizer = tokenizer

    def __len__(self):
        return len(self.pieces)

    def __getitem__(self, index):
        sequence = self.pieces[index]
        encoded = self.tokenizer.encode_plus(
            sequence, add_special_tokens=False, max_length=512
        )["input_ids"]
        return torch.LongTensor(encoded)


# -----------------------------------------------------------------------
# Checkpoint helpers
# -----------------------------------------------------------------------

NW_PREFIX = 'NW_'


def get_pickle_filename(out_dir, name):
    if name.startswith(NW_PREFIX):
        return os.path.join(out_dir, 'NW_all.pickle')
    return os.path.join(out_dir, f'{name}.pickle')


def load_existing_nw(out_dir):
    path = os.path.join(out_dir, 'NW__all.pickle')
    try:
        with open(path, 'rb') as f:
            d = pickle.load(f)
        return d if isinstance(d, dict) else {}
    except FileNotFoundError:
        return {}


def name_done(out_dir, name, nw_done, nw_buffer):
    if name.startswith(NW_PREFIX):
        return (name in nw_done) or (name in nw_buffer)
    return os.path.exists(get_pickle_filename(out_dir, name))


def save_pred(out_dir, name, prediction, nw_buffer):
    if name.startswith(NW_PREFIX):
        nw_buffer[name] = prediction
    else:
        with open(get_pickle_filename(out_dir, name), 'wb') as f:
            pickle.dump(prediction, f)


def flush_nw(out_dir, existing_nw, nw_buffer):
    if not nw_buffer:
        return
    merged = dict(existing_nw)
    merged.update(nw_buffer)
    with open(os.path.join(out_dir, 'NW__all.pickle'), 'wb') as f:
        pickle.dump(merged, f)
    print(f"Saved {len(nw_buffer)} NW sequence(s) to disk.")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta',      required=True,  help='Path to genome .fna / .fa file')
    parser.add_argument('--tokenizer',  required=True,  help='Path to tokenizer folder (6-new-12w-0)')
    parser.add_argument('--model',      required=True,  help='Path to model folder (dnabert_mm_fold_0_kouzine_g4)')
    parser.add_argument('--out_dir',    default='predictions_g4', help='Output directory for pickle files')
    parser.add_argument('--device',     default='cuda', help='"cuda", "cuda:1", "cpu", etc.')
    parser.add_argument('--batch_size', type=int, default=1)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Device
    if args.device.startswith('cuda') and not torch.cuda.is_available():
        print("CUDA not available, falling back to CPU.")
        args.device = 'cpu'
    device = torch.device(args.device)

    # Tokenizer
    print(f"Loading tokenizer from {args.tokenizer}")
    tokenizer = DNATokenizer.from_pretrained(args.tokenizer)

    # Model
    print(f"Loading model from {args.model}")
    config_path = os.path.join(args.tokenizer, 'config.json')
    config = BertConfig.from_pretrained(config_path)
    model = BertForTokenClassification.from_pretrained(args.model, config=config).to(device)
    model.eval()

    # Checkpointing state
    existing_nw = load_existing_nw(args.out_dir)
    nw_done = set(existing_nw.keys())
    nw_buffer = {}

    pred_ds = None
    pred_dataloader = None

    try:
        fasta_sequences = SeqIO.parse(open(args.fasta), 'fasta')
        for fasta in fasta_sequences:
            name, seq = fasta.id, str(fasta.seq)
            print(f"\nProcessing: {name}  (len={len(seq):,})")

            if name_done(args.out_dir, name, nw_done, nw_buffer):
                print(f"  Skipping — already processed.")
                continue

            del pred_ds, pred_dataloader
            gc.collect()

            pred_ds = PredDataset(seq, tokenizer)
            pred_dataloader = DataLoader(pred_ds, batch_size=args.batch_size, drop_last=False)

            cur_pred = []
            with torch.no_grad():
                for batch in tqdm(pred_dataloader, desc=f"  {name}"):
                    batch = batch.to(device)
                    logits = model(batch)['logits']
                    probs = torch.softmax(logits, dim=-1)[:, :, 1].cpu().numpy()
                    cur_pred.append(probs)

            final = stitch_np_seq(cur_pred)
            save_pred(args.out_dir, name, final, nw_buffer)
            print(f"  Saved prediction for {name}")

    finally:
        flush_nw(args.out_dir, existing_nw, nw_buffer)

    print("\nDone! All sequences processed.")


if __name__ == '__main__':
    main()
