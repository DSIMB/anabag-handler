"""
PDB Fetcher Script

Usage:
    python fetch_pdbs.py <path_to_mydataset_files> <pdb_output_dir>

Arguments:
    path_to_mydataset_files   Directory containing dataset files (fetcher_pdb_guide.csv, per_residue_informations, etc.)
    pdb_output_dir            Directory where PDB files will be downloaded

Example:
    python fetch_pdbs.py ./dataset_info ./pdb_downloads
"""
from urllib.request import urlretrieve
import os,sys
import numpy as np
import pandas as pd


def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    source code: https://stackoverflow.com/questions/37335759/using-python-to-download-specific-pdb-files-from-protein-data-bank
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(pdbcode, str(err), file=sys.stderr)
        return None

def fetch_pdbs(biological_unit_df,pdb_otuput_dir):
    l_output_paths = []
    updbs = biological_unit_df["PDB"].unique()
    for pdbid in updbs:
        output_path = download_pdb(pdbid, pdb_otuput_dir)
        if output_path is not None:
            print(f'Downloaded {pdbid} at {output_path}')
            l_output_paths.append(output_path)
    return l_output_paths

def main():
    if len(sys.argv) != 3 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        sys.exit(1)

    path_to_mydataset_files = sys.argv[1]
    pdb_output_dir = sys.argv[2]

    # Validate paths
    if not os.path.isdir(path_to_mydataset_files):
        print(f"Error: dataset directory not found: {path_to_mydataset_files}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(pdb_output_dir):
        print(f"Error: output directory not found: {pdb_output_dir}", file=sys.stderr)
        sys.exit(1)

    # Load dataset info
    fetcher_file = os.path.join(path_to_mydataset_files, 'fetcher_pdb_guide.csv')
    if not os.path.isfile(fetcher_file):
        print(f"Error: fetcher_pdb_guide.csv not found in {path_to_mydataset_files}", file=sys.stderr)
        sys.exit(1)

    biological_unit_df = pd.read_csv(fetcher_file)

    # Download PDBs
    l_output_paths = fetch_pdbs(biological_unit_df, pdb_output_dir)

    print("\nFinished downloading.")
    print(f"Total PDBs downloaded: {len(l_output_paths)}")
    print('You can now use the script to extract biological units and format the structures:')
    print('python extract_bus.py')


if __name__ == "__main__":
    main()