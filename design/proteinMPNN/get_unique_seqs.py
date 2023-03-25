import argparse
from pathlib import Path
import pandas as pd
import torch

def get_seqs_score_from_fasta(fname): 
    seqs=[]
    sample=[]
    score=[]
    global_score=[]
    seq_recovery=[]
    with open(fname, 'r') as f: 
        for i,l in enumerate(f): 
            if i == 0 or i==1:
                continue
            if l.startswith('>'):
                line = l.split(',')
                sample.append(int(line[1].strip().lstrip('sample=')))
                score.append(float(line[2].strip().lstrip('score=')))
                global_score.append(float(line[3].strip().lstrip('global_score=')))
                seq_recovery.append(float(line[4].strip().lstrip('seq_recovery=')))
            else: 
                seqs.append(l.strip('\n'))

    df = pd.DataFrame(list(zip(sample, seqs, score, global_score, seq_recovery)), columns=['sample', 'seq','score','global_score','seq_recovery'])
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_fasta",
        help="Path to input folder file",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-o", "--output_dir", help="Path to output directory", type=Path, required=True
    )
    args = parser.parse_args()

    # get unique sequences here from the pMPNN output
    if not args.input_fasta.exists():
        raise FileNotFoundError(args.input_fasta)
    args.output_dir.mkdir(exist_ok=True)

    df2 = get_seqs_score_from_fasta(args.input_fasta)

    name = str(args.input_fasta).split('/')[-1].rstrip('.fa')
    df_ls = [] 
    gr_idx = []
    with open(args.output_dir / f'{name}_unique_seqs.fa', 'w') as f:
        for i, (s, gr) in enumerate(df2.groupby('seq')):
            gr_idx.extend([i for _ in range(len(gr))])
            df_ls.append(gr)
            f.write(f'> {i}\n{s}\n')
            #print(s)

    df_grouped = pd.concat(df_ls).reset_index().drop(columns=['index'])
    df_grouped['seq_idx']=gr_idx
    df_grouped.to_csv(args.output_dir / f'{name}_unique_seqs.csv', index=False)
