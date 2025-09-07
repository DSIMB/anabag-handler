"""
PDB Organizer Script

Usage:
    python organize_pdbs.py <path_to_mydataset_files> <pdb_raw_pdbs_dir> <path_data>

Arguments:
    path_to_mydataset_files   Directory containing dataset files 
                              (fetcher_pdb_guide.csv, per_residue_information_AG.tsv, etc.)
    pdb_raw_pdbs_dir          Directory with raw PDB files (downloaded previously)
    path_data                 Output directory for processed structures and per-residue files ()

Example:
    python organize_pdbs.py ./dataset_info ./pdb_downloads ./my_dataset
"""
import sys, os 
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore') # ignore mixed types for pandas

def get_per_res_info(path_per_res_AG,path_per_res_AB,selected_ids):
    print('Reading per residue files (will take ~2 minutes)')
    df = pd.read_csv(path_per_res_AG,sep='\t',index_col=0,engine='c')
    df = pd.concat([df,pd.read_csv(path_per_res_AB,sep='\t',index_col=0,engine='c')])
    df = df[df['One_digit_id'].isin(selected_ids)]
    print('Done reading')
    return df

def parse_pdb(pdbfile):
    """Return header and ATOM content lines."""
    with open(pdbfile,'r') as fin:
        fullcontent = fin.readlines()
    header = [l for l in fullcontent if l.startswith('REMARK')]
    content = [l for l in fullcontent if l.startswith('ATOM')]
    return header, content


def shift_resindex(resindex_chain, start=1):
    """Renumber residues starting from `start`."""
    new_currindex = start
    newindex = []
    for i, resindex in enumerate(resindex_chain):
        if i > 0 and resindex != resindex_chain[i-1]:
            new_currindex += 1
        newindex.append('{:>4}'.format(new_currindex))
    return newindex, new_currindex


def renumber_rechain(atomschain, chain=None, rechain=False, start=1, initial=False):
    """
    Renumber atoms and residues.
    - If rechain=True and chain is given, force all residues to that chain (e.g. 'A', 'B').
    - If continuous=True, numbering does not reset at chain boundaries (used for full complexes).
    - start controls the first residue number.
    ATOM      1  N   ALA G   1
    """
    resindex = [l[22:27] for l in atomschain]
    resname = [l[17:20] for l in atomschain]
    chaini = [l[21:22] for l in atomschain]

    # Continuous vs chain-reset numbering
    newindex, final_resid = shift_resindex(resindex, start)

    # Atom indices
    newAtmindex = ['{:>7}'.format(i) for i in np.arange(1, len(newindex)+1)]

    # Chain assignment

    if rechain and chain is not None:
        c = np.repeat(chain, len(newAtmindex))
    else:
        c = chaini

    # Rebuild ATOM lines
    newatoms = [
        f"ATOM{newAtmindex[i]}{l[11:21]}{c[i]}{newindex[i]}{l[26:]}"
        for i, l in enumerate(atomschain)
    ]

    if initial:
        # Dataframe mapping
        start_ri = int(resindex[0])
        df_index_chain = pd.DataFrame([
            {
                'Chain_pdbfi': chaini[i],
                'AA':resname[i],
                'Ix_pdbi': int(resindex[i]) - (start_ri-1),
                'Ix_pdbff_full': int(newindex[i]),
            }
            for i in range(len(resindex))
        ])
    else:
        # Dataframe mapping
        df_index_chain = pd.DataFrame([
            {
                'Chain_pdbff': c[i],
                'AA':resname[i],
                'Ix_pdbff': int(newindex[i])-start+1,
                'Ix_pdbff_full': int(newindex[i]),
            }
            for i in range(len(resindex))
        ])

    return newatoms, df_index_chain.drop_duplicates(), final_resid



def chain_splitter(content, chains_molecule):
    """Return only the ATOM lines for the given chains."""
    return [l for l in content if l[21:22] in chains_molecule]



def organizer(path_raw_pdb_file, digit_id, path_output_directory, chains_AG, chains_AB, chain_options="both"):
    """

    """
    _, content = parse_pdb(path_raw_pdb_file)

    dfs = []  # collect dataframes

    ag_atoms = chain_splitter(content, chains_AG)
    ab_atoms = chain_splitter(content, chains_AB)
    # --------------------------------------------------
    # 1. Original chains (no rechaining)
    # --------------------------------------------------
    if chain_options in ["original", "both"]:

        # Renumber separately (no rechaining)
        ag_new, df_ag, final_rix_ag = renumber_rechain(ag_atoms, rechain=False,  start=1, initial=True)
        ab_new, df_ab, _ = renumber_rechain(ab_atoms, rechain=False,  start=final_rix_ag+1, initial=True)

        combined_atoms = ag_new + ab_new
        dfs.append(pd.concat([df_ag, df_ab],ignore_index=True))

        outpath = os.path.join(path_output_directory, f"{digit_id}-initial_chains.pdb")
        with open(outpath, "w") as out:
            out.writelines(combined_atoms)

    # --------------------------------------------------
    # 2. Rechained version (antigen=A, antibody=B)
    # --------------------------------------------------
    if chain_options in ["rechain", "both"]:

        # Antigen first, then antibody
        ag_new, df_ag, final_rix_ag = renumber_rechain(
            ag_atoms, chain="A", rechain=True, start=1
        )
        ab_new, df_ab, _ = renumber_rechain(
            ab_atoms, chain="B", rechain=True, start=final_rix_ag+1
        )

        combined_atoms = ag_new + ab_new
        dfs.append(pd.concat([df_ag, df_ab],ignore_index=True))

        outpath = os.path.join(path_output_directory, f"{digit_id}-formated_chains.pdb")
        with open(outpath, "w") as out:
            out.writelines(combined_atoms)

    if chain_options != "both":
        return dfs[0]
    
    df_processed = dfs[0].merge(dfs[1], on=['Ix_pdbff_full','AA'])
    df_processed['One_digit_id'] = digit_id
    return df_processed

def process_bu(row,pdb_outuput_pdbi_dir,path_data):
    chains_AG = row['Chains_AG'].split()
    chains_AB = row['Chains_AB'].split()
    pdbi = row['PDB']
    digit_id = row['One_digit_id']

    path_raw_pdb_file = os.path.join(pdb_outuput_pdbi_dir,pdbi + '.pdb')
    if not os.path.exists(path_raw_pdb_file):
        print(f'Did not found pdb file at {path_raw_pdb_file}')

    path_output_directory = os.path.join(path_data,digit_id)
    path_output_structure_dir = os.path.join(path_output_directory,'structure')
    path_output_file_dir = os.path.join(path_output_directory,'files')
    os.makedirs(path_output_structure_dir,exist_ok=True)
    os.makedirs(path_output_file_dir,exist_ok=True)

    df_processed = organizer(path_raw_pdb_file, digit_id, path_output_structure_dir, chains_AG, chains_AB, chain_options="both")
    return df_processed,path_output_file_dir

def merge_with_per_residue_write_files(digit_id,per_res_df,df_processed,path_output_file_dir):
    sdf = per_res_df[(per_res_df['One_digit_id'] == digit_id) & (per_res_df['Stat_res_pdbm'] == 'Solved')].copy()
    sdf['Ix_pdbff_full'] = np.arange(1,sdf.shape[0] +1)
    cols_to_drop = [c for c in df_processed.columns if c not in ['One_digit_id','Chain_pdbfi','AA','Ix_pdbff_full']]
    sdf = sdf.drop(columns=cols_to_drop)
    df_features =  df_processed.merge(sdf,on=['One_digit_id','Chain_pdbfi','AA','Ix_pdbff_full'])
    # write to destination
    df_features[df_features['Chain_pdbff']=='A'].to_csv(os.path.join(path_output_file_dir,"per_residue_information_AG.tsv"),sep='\t')
    df_features[df_features['Chain_pdbff']=='B'].to_csv(os.path.join(path_output_file_dir,"per_residue_information_AB.tsv"),sep='\t')
    

def main():
    if len(sys.argv) != 4 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        sys.exit(1)

    path_to_mydataset_files = sys.argv[1]
    pdb_raw_pdbs_dir = sys.argv[2]
    path_data = sys.argv[3]

    # Validate inputs
    if not os.path.isdir(path_to_mydataset_files):
        print(f"Error: dataset directory not found: {path_to_mydataset_files}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(pdb_raw_pdbs_dir):
        print(f"Error: raw PDB directory not found: {pdb_raw_pdbs_dir}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(path_data):
        print(f"Warning: output directory {path_data} not found, creating it.")
        os.makedirs(path_data, exist_ok=True)

    # Load per-residue info
    path_per_res_AG = os.path.join(path_to_mydataset_files, 'per_residue_information_AG.tsv')
    path_per_res_AB = os.path.join(path_to_mydataset_files, 'per_residue_information_AB.tsv')

    fetcher_file = os.path.join(path_to_mydataset_files, 'fetcher_pdb_guide.csv')
    if not os.path.isfile(fetcher_file):
        print(f"Error: fetcher_pdb_guide.csv not found in {path_to_mydataset_files}", file=sys.stderr)
        sys.exit(1)

    biological_unit_df = pd.read_csv(fetcher_file)
    selected_ids = biological_unit_df['One_digit_id'].to_list()
    per_res_df = get_per_res_info(path_per_res_AG, path_per_res_AB, selected_ids)

    # Process each row
    for _, row in biological_unit_df.iterrows():
        try:
            df_processed, path_output_file_dir = process_bu(row, pdb_raw_pdbs_dir, path_data)
            merge_with_per_residue_write_files(row['One_digit_id'], per_res_df, df_processed, path_output_file_dir)
            print(f"Merged per-residue info for {row['One_digit_id']}")
        except Exception as e:
            print(f"Error formating BU/merging per-residue info for {row['One_digit_id']}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
