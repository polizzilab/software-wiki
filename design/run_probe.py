'''
Example script for running probe to determine atomic contacts

Input: .pdb structure
Output: pandas DataFrame of contacts, with columns:
            chain1 resnum1 resname1 name1 chain2 resnum2 resname2 name2 contact_type
        where contact_type is one of:
            - wc    wide contact
            - cc    close contact
            - so    small overlap
            - bo    bad overlap (aka clashing)
            - wh    weak hydrogen bond
            - hb    hydrogen bond

Optional output:
    - .tsv          contacts table saved to a file
    - .probe        raw probe output saved to a file (you never really need this)

Notes:
🎵 probe does not know about segment ID's, so make sure your pdb file has unique chains

- jchang 2023-12-05
'''
import subprocess
from pathlib import Path
import pandas as pd

PROBE_EXEC = '/programs/x86_64-linux/probe/2.16.130520/probe'
PROBE_FLAGS = '-U -CON -Explicit -NOFACE -WEAKH -DE32 -WAT2wat -4 -ON -MC ALL ALL'

def run_probe(input_pdb : str,
              flags : str = PROBE_FLAGS,
              probe_exec : str = PROBE_EXEC,
              output_tsv : str = None,
              output_raw : str = None,
              verbose : bool = False):
    '''
    Parameters:
        input_pdb           file to analyze with probe

    Optional parameters:
        flags               flags to pass into probe (see probe -h), see PROBE_FLAGS for default
        probe_exec          path to probe executable, see PROBE_EXEC for default
        output_tsv          if provided, the parsed probe output is written here as a .tsv
        output_raw          if provided, the raw probe output is written to this file

    Returns: 
        pd.DataFrame        contacts among atoms in input_pdb
    '''
    if not Path(input_pdb).exists():
        raise FileNotFoundError(input_pdb)
    command = [probe_exec] + flags.split() + [str(input_pdb)]
    if verbose:
        print(' '.join(command))
    probe_out = subprocess.run(command, text=True, capture_output=True, check=True).stdout
    if output_raw is not None:
        Path(output_raw).parent.mkdir(exist_ok=True, parents=True)
        with open(output_raw, 'w') as f:
            f.write(probe_out)

    contacts_df = parse_probe_output(probe_out)
    if output_tsv is not None:
        Path(output_tsv).parent.mkdir(exist_ok=True, parents=True)
        with open(output_tsv, 'w') as f:
            contacts_df.to_csv(output_tsv, index=False, sep='\t')

    return contacts_df


def parse_probe_output(probe_output):
    '''
    Parse the probe output into a pd.DataFrame where each row corresponds to a pair of contacting atoms
    '''
    rows = []
    for line in probe_output.split('\n'):
        try:
            _, _, contact_type, atom1, atom2, *_ = line.split(':')
            chain1   = atom1[:2].strip()
            resnum1  = int(atom1[2:6])
            resname1 = atom1[6:10].strip()
            name1    = atom1[10:15].strip()
            chain2   = atom2[:2].strip()
            resnum2  = int(atom2[2:6])
            resname2 = atom2[6:10].strip()
            name2    = atom2[10:15].strip()
        except ValueError as e:
            # sys.stderr.write(f"failed to parse line: {line}")
            continue

        rows.append(dict(
            chain1=chain1, resnum1=resnum1, resname1=resname1, name1=name1,
            chain2=chain2, resnum2=resnum2, resname2=resname2, name2=name2,
            contact_type=contact_type,
        ))
    return pd.DataFrame(rows)


if __name__=='__main__':
    run_probe('probe_input/3ry2_A_BTN.pdb', output_tsv='probe_output/3ry2_A_BTN.tsv', verbose=True)
