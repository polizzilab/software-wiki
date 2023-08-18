"""
simple script to load all pdb files in a folder

TODO add some more documentation on how to use
TODO read the hostname from hostname(1)
TODO randomly sample ports until you find one that's unused

Written by jmou
Modified by jchang 2023-08-18
"""

import os
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--folder", type=str, default='./')
parser.add_argument("--port", type=int, default=9999)
parser.add_argument("--host", type=str, default='npl1.in.hwlab')
args = parser.parse_args()

pdbs = sorted(os.listdir(args.folder))
for pdb in pdbs:
    if pdb.endswith('.pdb'):
        print(f'load http://{args.host}:{args.port}/{pdb}')

os.system(f'python -m http.server {args.port}')
