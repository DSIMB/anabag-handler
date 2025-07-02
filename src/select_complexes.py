import pandas as pd
import sys, os
import numpy as np
import shutil


def dictionnary_of_selectors(df):
    for c in df.columns:
        if c in identifiers:
            continue
        elif c in ratios:
            print('{}\t=\t<{:.2f}# sign can be "<" or ">"'.format(c,df[c].mean()))
        elif c in counts:
            print('{}\t=\t{},{}# range: min, max'.format(c,df[c].min(),df[c].max()))
        elif c in lengths:
            print('{}\t=\t{},{}# range: min, max'.format(c,df[c].min(),df[c].max()))
        elif c in multiClasses:
            uclasses = list(pd.unique(df.loc[~pd.isnull(df[c]),c]))
            if len(uclasses) > 2:
                uclasses = uclasses[:2]
            print('{}\t=\t{}# classes separated by comas'.format(c,','.join(uclasses)))
        elif c in ranges:
            print('{}\t=\t{},{}# range: min, max'.format(c,df[c].min(),df[c].max()))
        elif c in thresholds:
            print('{}\t=\t<{:.2f}# sign can be "<" or ">"'.format(c,df[c].mean()))

# dictionnary_of_selectors(dfm)
# dictionnary_of_selectors(dfr)
# dictionnary_of_selectors(dfc)

# ======================== ======================== ======================== ======================== ======================== 
# ======================== ======================== ======================== ======================== ======================== 

internal = ['per_residue_info','rosetta_structures','formatted_structures','initial_structures','hetatm_structures']

def assign_values(c,param_field):
    if (c == 'SequenceIdentity') or (c =='InterfaceIdentity'):
        return ('redundancy',param_field.strip())
    
    elif (c in ratios) or (c in thresholds):
        signe = param_field.strip()[0]
        if '<>'.count(signe) == 0:
            print('Error in parameter field: did not find the sign for the threshold')
            print(c,param_field)
            return ('none',np.nan)
        value = float(param_field.strip()[1:])
        return ('threshold',signe,value)
    
    elif (c in counts) or (c in ranges) or (c in lengths): # range
        vals = param_field.strip().split(',')
        if len(vals) != 2:
            print('Error in parameter field: did not find two values speparated by a comma for the range')
            print(c,param_field)
            return  ('none',np.nan)
        valmin = float(vals[0])
        valmax = float(vals[1])
        return ('range',valmin,valmax)

    elif c in multiClasses:
        selclasses =  param_field.strip().split(',')
        if len(selclasses) == 0:
            print('Error in parameter field: class field empty (two or more classes must be separated by commas)')
            print(c,param_field)
            return  ('none',np.nan)
        return ('multiClass',selclasses)

    elif c in ['per_residue_info','rosetta_structures','formatted_structures','initial_structures','hetatm_structures']:
        if param_field.strip() == 'False':
            return ('internal',False)
        return ('internal',True)
    
    elif c == 'redundancy_mode':
        return ('internal',param_field.strip())


def read_configuration_complete(file_path):
    parameters_association = {'Antigen':'A','Antibody':'B','Complex':'C','Selection':'S'}
    parameters = {}
    with open(file_path, 'r') as fin:
        for l in fin:

            if l.startswith('\n') or l.startswith('#'):
                continue

            elif l.startswith('Parameters for:'):
                current_parameters_on = parameters_association[l.split(':')[1].strip().split()[0]]
                parameters[current_parameters_on] = {}
                continue
            else:
                lsplit = l.split('#')[0].split('=')
                parameter = lsplit[0].strip()
                param_field = lsplit[1].strip()

                parameters[current_parameters_on][parameter] = assign_values(parameter,param_field)

    if not ('S' in list(parameters.keys())):
        parameters['S'] = {}

    for c in internal + ['redundancy_mode']:
        if not (c in list(parameters['S'].keys())):
            if c == 'redundancy_mode':
                parameters['S'][c] = 'Intersection'
            else:
                parameters['S'][c] = False
    return parameters

def select_numerical_criterion_range(df,criterion,lowerLimit,higherLimit):
    """
    takes as input a numerical criterion with a lower and higher limit (like SASA_interface or UA_Intramembrane)
    select all complexes that satisfy the condition
    return a list of One_digit_id
    """
    lids = pd.unique(df.loc[(df[criterion] >=lowerLimit) & (df[criterion] <= higherLimit),'One_digit_id'])
    print(f'Found {len(lids)} complexes where {criterion} is between {lowerLimit} and {higherLimit}')
    return lids

def select_numerical_criterion_threshold(df,criterion,signe,threshold):
    """
    takes as input a numerical criterion with a lower and higher limit (like SASA_interface or UA_Intramembrane)
    select all complexes that satisfy the condition
    return a list of One_digit_id
    """
    if signe == '<':
        lids = pd.unique(df.loc[(df[criterion] <=threshold),'One_digit_id'])
        print(f'Found {len(lids)} complexes where {criterion} is <= {threshold} ')
    else:
        lids = pd.unique(df.loc[(df[criterion] > threshold),'One_digit_id'])
        print(f'Found {len(lids)} complexes where {criterion} is > {threshold} ')

    return lids

def select_multiClass_criterion(df,criterion,values):
    """
    takes as input a numerical criterion with a lower and higher limit (like SASA_interface or UA_Intramembrane)
    select all complexes that satisfy the condition
    return a list of One_digit_id
    """
    lids = pd.unique(df.loc[(df[criterion].isin(values)),'One_digit_id'])
    print(f'Found {len(lids)} complexes where {criterion} is in {" ".join(values)}')
    return lids

def select_group_criterion(df,group_criterion,pool_complexes):
    """
    takes as input a criterion that defines groups (like a sequence identity grouping)
    select One complex for each group
    return a list of One_digit_id
    """
    subdf = df.loc[df['One_digit_id'].isin(pool_complexes)]
    lpdbs = []
    for groupID, gdf in subdf.groupby(group_criterion):
        # random selection of one complex
        lpdbs.append(np.random.choice(gdf['One_digit_id'].to_list()))
    print(f'Selected {len(lpdbs)} complexes among the available groups')
    return lpdbs

def create_new_grouping(df,grouping_cols):
    grouping_col = '_'.join(grouping_cols)
    df[grouping_col] = df[grouping_cols[0]].astype(str)
    if len(grouping_cols) > 1:
        for gi in range(1,len(grouping_cols)):
            df[grouping_col] = df[grouping_col] + '_' + df[grouping_cols[gi]].astype(str)

    print(f'Found {pd.unique(df[grouping_col]).shape[0]} groups based on the following groupings: {", ".join(grouping_cols)}')    
    return df,grouping_col

def send_parameters(c,values,df_curr):
    mode_ = values[0]
    if mode_ == 'threshold':
        return select_numerical_criterion_threshold(df_curr,c,values[1],values[2])
    elif mode_ == 'range':
        return select_numerical_criterion_range(df_curr,c,values[1],values[2])
    elif mode_ == 'multiClass':
        return select_multiClass_criterion(df_curr,c,values[1])
    else:
        print('Error, found incorrect parameter')
        print(c,values)
        return []


def select_ids(parameters,dfm,dfr,dfc):
    grouping_cols = []
    selected_ids = []
    for chain,para in parameters.items():
        if chain == 'S':
            continue
        for column,values in para.items():
            if (column == 'SequenceIdentity') or (column == 'InterfaceIdentity'):
                grouping_cols.append(values[1])
                continue

            elif column in dfm.columns:
                df_curr = dfm.loc[dfm['Chain_pdbff'] == chain]
            elif column in dfr.columns:
                df_curr = dfr.copy(deep=True)
            elif column in dfc.columns:
                df_curr = dfc.copy(deep=True)
            else:
                print('Error, parameter used in config not found in dataframes')
                print(column,values)

            digit_ids = send_parameters(column,values,df_curr)

            if len(selected_ids) == 0:
                selected_ids.extend(digit_ids)
            else:
                selected_ids = list(np.intersect1d(selected_ids,digit_ids))
                print(f'Number of selected complexes that satisfies all criterions: {len(selected_ids)}',end='\n\n')


    print(grouping_cols)
    if (len(grouping_cols) > 0) :
        dfgrouping,grouping_col = create_new_grouping(dfc.copy(deep=True),grouping_cols)
        selected_ids = select_group_criterion(dfgrouping,grouping_col,selected_ids)
    
    return selected_ids




def copy_data(selected_ids,dfm,dfc,dfr,selectionParameters):
    path_to_mydataset_structures = os.path.join(PATHANABAG,'my_dataset','structures')
    path_to_mydataset_files = os.path.join(PATHANABAG,'my_dataset','files')
    labs = []
    lags = []
    print('Copying structure files')
    for digit_id in selected_ids:
        path_folder = os.path.join(PATHANABAG,'data',digit_id)
        if selectionParameters['per_residue_info'][1]:
            lags.append(os.path.join(path_folder,'files','per_residue_information_AG.tsv'))
            labs.append(os.path.join(path_folder,'files','per_residue_information_AB.tsv'))
        if selectionParameters['rosetta_structures'][1]:
            agr = os.path.join(path_folder,'structures',digit_id + '_AG-rosetta_relaxed.pdb')
            abr = os.path.join(path_folder,'structures',digit_id + '_AB-rosetta_relaxed.pdb')
            cplx = os.path.join(path_folder,'structures',digit_id + '-rosetta_relaxed.pdb')
            shutil.copy2(agr,path_to_mydataset_structures)
            shutil.copy2(abr,path_to_mydataset_structures)
            shutil.copy2(cplx,path_to_mydataset_structures)
        if selectionParameters['formatted_structures'][1]:
            cplx = os.path.join(path_folder,'structures',digit_id + '-formated_chains.pdb')
            shutil.copy2(cplx,path_to_mydataset_structures)
        if selectionParameters['initial_structures'][1]:
            cplx = os.path.join(path_folder,'structures',digit_id + '-initial_chains.pdb')
            shutil.copy2(cplx,path_to_mydataset_structures)
        if selectionParameters['hetatm_structures'][1]:
            cplx = os.path.join(path_folder,'structures',digit_id + '-initial_chains_hetatm.pdb')
            shutil.copy2(cplx,path_to_mydataset_structures)

    labdf,lagdf = [],[]
    if selectionParameters['per_residue_info'][1]:
        print('Copying per residue information')
        for i in range(len(labs)):
            labdf.append(pd.read_csv(labs[i],sep='\t',engine='c',index_col=0))
            lagdf.append(pd.read_csv(lags[i],sep='\t',engine='c',index_col=0))
    lagdf = pd.concat(lagdf,ignore_index=True)
    labdf = pd.concat(labdf,ignore_index=True)
    lagdf.to_csv(os.path.join(path_to_mydataset_files , 'per_residue_information_AG.tsv'),sep='\t')
    labdf.to_csv(os.path.join(path_to_mydataset_files , 'per_residue_information_AB.tsv'),sep='\t')

    dfm.loc[dfm['One_digit_id'].isin(selected_ids)].to_csv(os.path.join(path_to_mydataset_files , 'per_chain_pdbff_informations.tsv'),sep='\t')
    dfc.loc[dfc['One_digit_id'].isin(selected_ids)].to_csv(os.path.join(path_to_mydataset_files , 'cluster_informations.tsv'),sep='\t')
    dfr.loc[dfr['One_digit_id'].isin(selected_ids)].to_csv(os.path.join(path_to_mydataset_files , 'method_resolution.tsv'),sep='\t')

# ======================== ======================== ======================== ======================== ======================== 
# ======================== ======================== ======================== ======================== ======================== 


        
PATHANABAG = sys.argv[1] # '/dsimb/abbesses/grand/Documents/Articles/PPI_AGAB_V2/ANABAG'
PATHANABAG_dataset_info = os.path.join(PATHANABAG,'dataset_info')
dfm = pd.read_csv(os.path.join(PATHANABAG_dataset_info,'per_chain_pdbff_informations.tsv'),sep='\t',index_col=0,engine='c')
dfc = pd.read_csv(os.path.join(PATHANABAG_dataset_info,'cluster_informations.tsv'),sep='\t',index_col=0,engine='c')
dfr = pd.read_csv(os.path.join(PATHANABAG_dataset_info,'method_resolution.tsv'),sep='\t',index_col=0,engine='c')

identifiers = ['One_digit_id','Chain_pdbff',] + list(np.setdiff1d(dfc.columns,['AbClass','Sequence_length_Ag', 'Sequence_length_Ab']))
ratios = ['HelicesPercent[G,H,I]','StrandsPercent[E,B]','LoopsPercent[E,B]']
counts = ['Helices[G,H,I]','Strands[E,B]','Loops[E,B]',
        'UA_Active_site','UA_Binding_site','UA_Cross-link','UA_DNA_binding',
        'UA_Glycosylation','UA_Intramembrane','UA_Lipidation','UA_Natural_variant',
        'UA_Repeat','UA_Site','UA_Transmembrane','UA_Zinc_finger', 'Number_of_glycans','Number_of_modelled_residues']
lengths = ['Number_of_residues','Sequence_length_Ag', 'Sequence_length_Ab',]
multiClasses = [ 'AbClass', 'method', 'organism', 'affinity_method']
ranges = ['resolution','temperature', 'pH_pHcrystal',
        '⟨Z⟩_Complex_pH7','⟨Z⟩_Monomer_pH7','⟨Z⟩_Complex_pHcrystal','⟨Z⟩_Monomer_pHcrystal',
        '⟨Z²⟩-⟨Z⟩²_Complex_pH7','⟨Z²⟩-⟨Z⟩²_Monomer_pH7','⟨Z²⟩-⟨Z⟩²_Complex_pHcrystal','⟨Z²⟩-⟨Z⟩²_Monomer_pHcrystal',
        '⟨μ⟩_Complex_pH7','⟨μ⟩_Monomer_pH7','⟨μ⟩_Complex_pHcrystal','⟨μ⟩_Monomer_pHcrystal',
        '⟨μ²⟩-⟨μ⟩²_Complex_pH7','⟨μ²⟩-⟨μ⟩²_Monomer_pH7','⟨μ²⟩-⟨μ⟩²_Complex_pHcrystal','⟨μ²⟩-⟨μ⟩²_Monomer_pHcrystal', 
        'SASA_Complex','SASA_Monomer','SASA_interface',]
thresholds = ['affinity','delta_g']
rosetta_columns = ['yhh_planarity_Complex_crystal','yhh_planarity_Complex_relaxed','yhh_planarity_Monomer_crystal','yhh_planarity_Monomer_relaxed','rama_prepro_Complex_crystal','rama_prepro_Complex_relaxed','rama_prepro_Monomer_crystal','rama_prepro_Monomer_relaxed','ref_Complex_crystal','ref_Complex_relaxed','ref_Monomer_crystal','ref_Monomer_relaxed','p_aa_pp_Complex_crystal','p_aa_pp_Complex_relaxed','p_aa_pp_Monomer_crystal','p_aa_pp_Monomer_relaxed','dslf_fa13_Complex_crystal','dslf_fa13_Complex_relaxed','dslf_fa13_Monomer_crystal','dslf_fa13_Monomer_relaxed','fa_atr_Complex_crystal','fa_atr_Complex_relaxed','fa_atr_Monomer_crystal','fa_atr_Monomer_relaxed','fa_dun_Complex_crystal','fa_dun_Complex_relaxed','fa_dun_Monomer_crystal','fa_dun_Monomer_relaxed','fa_elec_Complex_crystal','fa_elec_Complex_relaxed','fa_elec_Monomer_crystal','fa_elec_Monomer_relaxed','fa_intra_rep_Complex_crystal','fa_intra_rep_Complex_relaxed','fa_intra_rep_Monomer_crystal','fa_intra_rep_Monomer_relaxed','fa_intra_sol_xover4_Complex_crystal','fa_intra_sol_xover4_Complex_relaxed','fa_intra_sol_xover4_Monomer_crystal','fa_intra_sol_xover4_Monomer_relaxed','fa_rep_Complex_crystal','fa_rep_Complex_relaxed','fa_rep_Monomer_crystal','fa_rep_Monomer_relaxed','fa_sol_Complex_crystal','fa_sol_Complex_relaxed','fa_sol_Monomer_crystal','fa_sol_Monomer_relaxed','hbond_bb_sc_Complex_crystal','hbond_bb_sc_Complex_relaxed','hbond_bb_sc_Monomer_crystal','hbond_bb_sc_Monomer_relaxed','hbond_lr_bb_Complex_crystal','hbond_lr_bb_Complex_relaxed','hbond_lr_bb_Monomer_crystal','hbond_lr_bb_Monomer_relaxed','hbond_sc_Complex_crystal','hbond_sc_Complex_relaxed','hbond_sc_Monomer_crystal','hbond_sc_Monomer_relaxed','hbond_sr_bb_Complex_crystal','hbond_sr_bb_Complex_relaxed','hbond_sr_bb_Monomer_crystal','lk_ball_wtd_Complex_crystal','lk_ball_wtd_Complex_relaxed','lk_ball_wtd_Monomer_crystal','lk_ball_wtd_Monomer_relaxed','omega_Complex_crystal','omega_Complex_relaxed','omega_Monomer_crystal','omega_Monomer_relaxed','pro_close_Complex_crystal','pro_close_Complex_relaxed','pro_close_Monomer_crystal','pro_close_Monomer_relaxed','score_Complex_crystal','score_Complex_relaxed','score_Monomer_crystal','score_Monomer_relaxed']


file_path = sys.argv[2]
parameters = read_configuration_complete(file_path)
selected_ids = select_ids(parameters,dfm,dfr,dfc)
copy_data(selected_ids,dfm,dfc,dfr,parameters['S'])


