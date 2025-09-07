import sys, os 
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore') # ignore mixed types for pandas

def get_per_res_info(path_per_res_AG,path_per_res_AB,selected_ids,usecols):
    print('Reading per residue files (will take ~2 minutes)')
    df = pd.read_csv(path_per_res_AG,sep='\t',usecols=usecols,engine='c')
    df = pd.concat([df,pd.read_csv(path_per_res_AB,sep='\t',usecols=usecols,engine='c')])
    df = df[df['One_digit_id'].isin(selected_ids)]
    return df

def get_guide_fetcher(path_per_res_AG,path_per_res_AB,selected_ids,output_fetcher=None):
    usecols = ['Chain_pdbff','Chain_pdbfi','One_digit_id']
    df = get_per_res_info(path_per_res_AG,path_per_res_AB,selected_ids,usecols)

    biological_unit_df = []
    for digit,ddf in df.groupby('One_digit_id'):
        cag = ddf[ddf['Chain_pdbff'] == 'A']['Chain_pdbfi'].unique()
        cab = ddf[ddf['Chain_pdbff'] == 'B']['Chain_pdbfi'].unique()
        biological_unit_df.append({
            "PDB":digit[:4],
            "One_digit_id":digit,
            "Chains_AG":' '.join(cag),
            "Chains_AB":' '.join(cab),
        })
    # generer le fichier guide pour ces identifiants
    biological_unit_df = pd.DataFrame(biological_unit_df)
    if output_fetcher is not None:
        print(f'Writting fetcher file in {output_fetcher} ...')
        biological_unit_df.to_csv(output_fetcher,index=False)
        print(f'Fetcher file written in {output_fetcher}')
    return biological_unit_df



# liste d'identifiant selected_ids

if __name__ == '__main__':
    selected_ids = ['8dis_0', '8kdr_0', '8sdg_0', '6vzi_0', '3ezj_1', '5lwf_1','8cxc_1', '8gou_0', '8da1_0', '6dbg_0']
    path_per_res_AG = '.../per_residue_information_AG.tsv'
    path_per_res_AB = '.../per_residue_information_AB.tsv'
    biological_unit_df = get_guide_fetcher(path_per_res_AG,path_per_res_AB,selected_ids)