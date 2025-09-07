# Antibody-Antigen Complex Database Features

This document describes all features available in the antibody-antigen complex database. Features are calculated per-residue and stored in TSV format.

## Data Overview

- **Total Features**: 197

## Table of Contents

- [Sequence Clustering](#sequence-clustering)
- [Identifiers and Indices](#identifiers-and-indices)
- [Structural Features](#structural-features)
- [Solvent Accessible Surface Area (SASA)](#solvent-accessible-surface-area-sasa)
- [UniProt Annotations](#uniprot-annotations)
- [Charge and pH](#charge-and-ph)
- [Network Centrality](#network-centrality)
- [Flexibility](#flexibility)
- [Detected Glycosylations Features](#detected-glycosylations-features)
- [Rosetta Energy Terms](#rosetta-energy-terms)

## Sequence Clustering

### `SG100AG`

Cluster ID for sequences sharing 100% identity among antigens

### `SG100AGc8`

Cluster ID for sequences sharing 100% identity among antigens and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG20AG`

Cluster ID for sequences sharing 20% identity among antigens

### `SG20AGc8`

Cluster ID for sequences sharing 20% identity among antigens and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG40AG`

Cluster ID for sequences sharing 40% identity among antigens

### `SG40AGc8`

Cluster ID for sequences sharing 40% identity among antigens and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG60AG`

Cluster ID for sequences sharing 60% identity among antigens

### `SG60AGc8`

Cluster ID for sequences sharing 60% identity among antigens and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG80AG`

Cluster ID for sequences sharing 80% identity among antigens

### `SG80AGc8`

Cluster ID for sequences sharing 80% identity among antigens and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG95AG`

Cluster ID for sequences sharing 95% identity among antigens

### `SG95AGc8`

Cluster ID for sequences sharing 95% identity among antigens and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG100AB`

Cluster ID for sequences sharing 100% identity among antibodies

### `SG100ABc8`

Cluster ID for sequences sharing 100% identity among antibodies and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG


### `SG60AB`

Cluster ID for sequences sharing 60% identity among antibodies

### `SG60ABc8`

Cluster ID for sequences sharing 60% identity among antibodies and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG80AB`

Cluster ID for sequences sharing 80% identity among antibodies

### `SG80ABc8`

Cluster ID for sequences sharing 80% identity among antibodies and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

### `SG95AB`

Cluster ID for sequences sharing 95% identity among antibodies

### `SG95ABc8`

Cluster ID for sequences sharing 95% identity among antibodies and >= 80% coverage
c8 denotes the default parameters for MMseq2. Briefly, two sequence sharing >X% identity and less than 80% coverage will be clustered appart. If you are after a size independant clustering, refer to the default clustering used in ANABAG SGxAB or SGxAG

## Identifiers and Indices

### `AA`

Amino acid type (three-letter code)

**Data Profile:**
- **Data Type**: `object`
- **Unique Values**: 20 (0.00% of non-null)

### `AlignIx_SG100AB`

Alignment position index within the SG100AB sequence cluster

### `AlignIx_SG100ABc8`

Alignment position index within the SG100ABc8 sequence cluster

### `AlignIx_SG60AB`

Alignment position index within the SG60AB sequence cluster

### `AlignIx_SG60ABc8`

Alignment position index within the SG60ABc8 sequence cluster

### `AlignIx_SG80AB`

Alignment position index within the SG80AB sequence cluster

### `AlignIx_SG80ABc8`

Alignment position index within the SG80ABc8 sequence cluster

### `AlignIx_SG95AB`

Alignment position index within the SG95AB sequence cluster

### `AlignIx_SG95ABc8`

Alignment position index within the SG95ABc8 sequence cluster

### `Align_index`

Alignment index with UniProt sequence

### `Annotated`

Flag indicating if residue has UniProt annotations

**Data Profile:**
- **Data Type**: `bool`
- **Binary Values**: False, True

### `Chain`

Chain identifier (initial chain)

### `Chain_pdbff`

Chain identifier in formatted PDB (B:Antibody, A:Antigen format)

**Data Profile:**
- **Data Type**: `object`
- **Values**: A,B

### `Chain_pdbfi`

Chain identifier in initial PDB file

**Data Profile:**
- **Data Type**: `object`
- **Unique Values**: 48 (0.00% of non-null)
- **Most Common**: H , L , B , E 
- **Other Values**: 44 additional unique values

### `Chain_pdbm`

Chain identifier in modeled PDB (same as initial chain)

### `Full_id`

Full unique identifier of the biological unit (PDB + chains, e.g., 9bu6_DCAB)

### `Identifier`

The user can skip
Internal: General identifier (may contain NaN)

### `Ix_fasta`

Residue index in FASTA sequence

### `Ix_pdbff`

Residue index in formated PDB file. Format pdbff is chains A,B for antigen, antibodies 
(formatted numbering) [1->N per chain]

### `Ix_pdbff_full`

Residue index in PDB file with all chains 
(formatted numbering)[1->N , N= #residues in complex]


### `Ix_pdbfi`

Residue index in formated PDB file with initial chains
(formatted numbering) [1->N per chain]

### `Ix_pdbfi_full`

Residue index in PDB file with all chains Same as Ix_pdbff_full

### `Ix_pdbi`

Residue index in initial PDB structure

### `Ix_pdbm`

Residue index in modeled PDB structure

### `Ix_unip`

Residue index in UniProt sequence

### `One_digit_id`

Single character identifier for the complex (PDBcode_digit). e.g 1a14_1

### `PDB`

PDB identifier of the structure

### `Res_fasta`

Residue at same index in FASTA sequence of PDB

### `Res_pdbi`

Residue name in initial PDB

### `Res_pdbm`

Residue name in modeled PDB

### `Res_unip`

Residue at same index in UniProt sequence

### `SeqUnip_SeqPdbm`

Sequence comparison between UniProt and PDB structure

**Data Profile:**
- **Data Type**: `object`
- **Unique Values**: 2 
- **Values**: Divergent, Identical

### `Stat_res_pdbi`

Status of residue in initial PDB 

**Data Profile:**
- **Data Type**: `object`
- **Values**: Missing, Solved

### `Stat_res_pdbm`

Status of residue in modeled PDB

**Data Profile:**
- **Data Type**: `object`
- **Values**: Modelled, Solved

### `UA_UAlign_index`

UniProt alignment index position

## Structural Features

### `phi`

Phi backbone dihedral angle (degrees)

### `psi`

Psi backbone dihedral angle (degrees)

### `structureS`

Secondary structure assignment by DSSP

**Data Profile:**
- **Data Type**: `object`
- **Values**: B, E, G, H, I ... (7 total)

### `Is_disulfide`

Binary flag indicating if residue forms disulfide bond

**Data Profile:**
- **Binary Values**: True,False

## Solvent Accessible Surface Area (SASA) and related features

### `ClassInterface_Model`

Interface residue classification (Core, Rim, Support, or Outside) based on SASA
Core: rSASA_Monomer > 25 , rSASA_Complex < 25
Support: rSASA_Monomer < 25 , rSASA_Complex < rSASA_Monomer
Rim: rSASA_Monomer > 25 , 25 < rSASA_Complex < rSASA_Monomer

**Data Profile:**
- **Values**: Core, Outside, Rim, Support

### `caDepth`

C-alpha atom depth from surface

### `rDepth`

Residue depth from surface

### `SASA_Complex`

Per-residue solvent accessible surface area in the bound complex

### `SASA_Difference`

Change in solvent accessible surface area upon complex formation (Complex - Monomer)

### `SASA_Monomer`

Per-residue solvent accessible surface area in the unbound monomer

### `SASA_apolar_Complex`

Per-residue apolar atom solvent accessible surface area in the bound complex

### `SASA_apolar_Monomer`

Per-residue apolar atom solvent accessible surface area in the unbound monomer

### `SASA_chain_Complex`

Per-residue backbone solvent accessible surface area in the bound complex

### `SASA_chain_Monomer`

Per-residue backbone solvent accessible surface area in the unbound monomer

### `SASA_polar_Complex`

Per-residue polar atom solvent accessible surface area in the bound complex

### `SASA_polar_Monomer`

Per-residue polar atom solvent accessible surface area in the unbound monomer

### `SASA_side_Complex`

Per-residue sidechain solvent accessible surface area in the bound complex

### `SASA_side_Monomer`

Per-residue sidechain solvent accessible surface area in the unbound monomer

## UniProt Annotations

### `UA_Active_site`

UniProt annotation marking catalytic/active site residues

**Data Profile:**
- **Values**: False, True

### `UA_Beta_strand`

UniProt secondary structure annotation for beta strands

**Data Profile:**
- **Values**: False, True

### `UA_Binding_site`

UniProt annotation marking ligand/substrate binding sites

**Data Profile:**
- **Values**: False, True

### `UA_Coiled_coil`

UniProt annotation for coiled-coil regions

**Data Profile:**
- **Values**: False, True

### `UA_Cross-link`

UniProt annotation for cross-linked residues

**Data Profile:**
- **Values**: False, True

### `UA_Cross-link_with`

UniProt annotation of partner residue in cross-link


### `UA_DNA_binding`

UniProt annotation for DNA-binding regions

**Data Profile:**
- **Values**: False, True

### `UA_Disulfide_bond`

UniProt annotation indicating if residue forms a disulfide bond

**Data Profile:**
- **Values**: False, True

### `UA_Disulfide_bond_with`

UniProt annotation of partner residue in disulfide bond

### `UA_Glycosylation`

UniProt annotation marking glycosylation sites

**Data Profile:**
- **Values**: False, True

### `UA_Helix`

UniProt secondary structure annotation for helices

**Data Profile:**
- **Values**: False, True

### `UA_Intramembrane`

UniProt annotation for intramembrane regions

**Data Profile:**
- **Values**: False, True

### `UA_Lipidation`

UniProt annotation for lipidation sites

**Data Profile:**
- **Values**: False, True

### `UA_Modified_residue`

UniProt annotation of post-translational modifications

**Data Profile:**
- **Values**: False, True

### `UA_Mutagenesis`

UniProt annotation of experimentally studied mutations

**Data Profile:**
- **Values**: False, True

### `UA_Natural_variant`

UniProt annotation of naturally occurring sequence variants

**Data Profile:**
- **Values**: False, True

### `UA_Repeat`

UniProt annotation for repeat regions

**Data Profile:**
- **Values**: False, True

### `UA_Site`

UniProt annotation for functionally important sites

**Data Profile:**
- **Values**: False, True

### `UA_Transmembrane`

UniProt annotation for transmembrane regions

**Data Profile:**
- **Values**: False, True

### `UA_Turn`

UniProt secondary structure annotation for turns

**Data Profile:**
- **Values**: False, True

### `UA_Zinc_finger`

UniProt annotation for zinc finger regions

**Data Profile:**
- **Values**: False, True

## Charge and pH

### `Charge_Complex_pH7`

Residue charge in complex at pH 7

### `Charge_Complex_pHcrystal`

Residue charge in complex at crystallization pH

### `Charge_Monomer_pH7`

Residue charge in monomer at pH 7

### `Charge_Monomer_pHcrystal`

Residue charge in monomer at crystallization pH

### `Charge_ideal_pH7`

Ideal behaviour residue charge at pH 7

### `Charge_ideal_pHcrystal`

Ideal behaviour residue charge at crystallization pH

### `pH_pHcrystal`

Experimental pH value (for solving the structure).

## Network Centrality

### `Betweenness_centrality_Complex`

Betweenness centrality in residue contact network of complex

### `Betweenness_centrality_Monomer`

Betweenness centrality in residue contact network of monomer

### `Betweenness_centrality_titratable_Complex`

Betweenness centrality in titratable residue network of complex

### `Betweenness_centrality_titratable_Monomer`

Betweenness centrality in titratable residue network of monomer

### `Closeness_centrality_Complex`

Closeness centrality in residue contact network of complex

### `Closeness_centrality_Monomer`

Closeness centrality in residue contact network of monomer

### `Closeness_centrality_titratable_Complex`

Closeness centrality in titratable residue network of complex

### `Closeness_centrality_titratable_Monomer`

Closeness centrality in titratable residue network of monomer

### `Degree_centrality_Complex`

Degree centrality in residue contact network of complex

### `Degree_centrality_Monomer`

Degree centrality in residue contact network of monomer

### `Degree_centrality_titratable_Complex`

Degree centrality in titratable residue network of complex

### `Degree_centrality_titratable_Monomer`

Degree centrality in titratable residue network of monomer

### `EigenCentrality_centrality_Complex`

Eigenvector centrality in residue contact network of complex

### `EigenCentrality_centrality_Monomer`

Eigenvector centrality in residue contact network of monomer

### `EigenCentrality_centrality_titratable_Complex`

Eigenvector centrality in titratable residue network of complex

### `EigenCentrality_centrality_titratable_Monomer`

Eigenvector centrality in titratable residue network of monomer

## Flexibility

### `mean_MEAN_LDDT`

mean of local distance difference test (lDDT) score
PREDICTED (by PEGASUS)

### `mean_RMSF`

Root mean square fluctuation mean
PREDICTED (by PEGASUS)

### `mean_STD_PHI`

mean of phi angle standard deviation
PREDICTED (by PEGASUS)

### `mean_STD_PSI`

mean of psi angle standard deviation
PREDICTED (by PEGASUS)

### `std_MEAN_LDDT`

standard deviation of local distance difference test (lDDT) score
PREDICTED (by PEGASUS)

### `std_RMSF`

Root mean square fluctuation standard deviation
PREDICTED (by PEGASUS)

### `std_STD_PHI`

standard deviation of phi angle standard deviation
PREDICTED (by PEGASUS)

### `std_STD_PSI`

standard deviation of psi angle standard deviation

## Detected Glycosylations Features

### `Sugar`

Flag indicating the sugar-residue type which is linked to the structure 
DETECTED (by Privateer)

**Data Profile:**
- **Values**: NAG, GLC, ..., NaN

### `SugarChain`

Chain identifier of attached sugar/glycan

**Data Profile:**
- **Values**: A, B, C, F, G ... (8 total)

### `SugarLinkedR`

Residue linked to sugar/glycan

**Data Profile:**
- **Values**: ASN, ..

### `Sugar_Ix`

Index of sugar/glycan modification


## Rosetta Energy Terms

### `dslf_fa13_Complex_crystal`

Disulfide geometry potential calculated for the bound complex using the crystal structure

### `dslf_fa13_Complex_relaxed`

Disulfide geometry potential calculated for the bound complex using the energy-minimized structure

### `dslf_fa13_Monomer_crystal`

Disulfide geometry potential calculated for the unbound monomer using the crystal structure

### `dslf_fa13_Monomer_relaxed`

Disulfide geometry potential calculated for the unbound monomer using the energy-minimized structure

### `fa_atr_Complex_crystal`

Lennard-Jones attractive energy between atoms in different residues calculated for the bound complex using the crystal structure

### `fa_atr_Complex_relaxed`

Lennard-Jones attractive energy between atoms in different residues calculated for the bound complex using the energy-minimized structure

### `fa_atr_Monomer_crystal`

Lennard-Jones attractive energy between atoms in different residues calculated for the unbound monomer using the crystal structure

### `fa_atr_Monomer_relaxed`

Lennard-Jones attractive energy between atoms in different residues calculated for the unbound monomer using the energy-minimized structure

### `fa_dun_Complex_crystal`

Internal energy of sidechain rotamers (Dunbrack 2010) calculated for the bound complex using the crystal structure

### `fa_dun_Complex_relaxed`

Internal energy of sidechain rotamers (Dunbrack 2010) calculated for the bound complex using the energy-minimized structure

### `fa_dun_Monomer_crystal`

Internal energy of sidechain rotamers (Dunbrack 2010) calculated for the unbound monomer using the crystal structure

### `fa_dun_Monomer_relaxed`

Internal energy of sidechain rotamers (Dunbrack 2010) calculated for the unbound monomer using the energy-minimized structure

### `fa_elec_Complex_crystal`

Coulombic electrostatic potential with distance-dependent dielectric calculated for the bound complex using the crystal structure

### `fa_elec_Complex_relaxed`

Coulombic electrostatic potential with distance-dependent dielectric calculated for the bound complex using the energy-minimized structure

### `fa_elec_Monomer_crystal`

Coulombic electrostatic potential with distance-dependent dielectric calculated for the unbound monomer using the crystal structure

### `fa_elec_Monomer_relaxed`

Coulombic electrostatic potential with distance-dependent dielectric calculated for the unbound monomer using the energy-minimized structure

### `fa_intra_rep_Complex_crystal`

Lennard-Jones repulsive energy between atoms in the same residue calculated for the bound complex using the crystal structure

### `fa_intra_rep_Complex_relaxed`

Lennard-Jones repulsive energy between atoms in the same residue calculated for the bound complex using the energy-minimized structure

### `fa_intra_rep_Monomer_crystal`

Lennard-Jones repulsive energy between atoms in the same residue calculated for the unbound monomer using the crystal structure

### `fa_intra_rep_Monomer_relaxed`

Lennard-Jones repulsive energy between atoms in the same residue calculated for the unbound monomer using the energy-minimized structure

### `fa_intra_sol_xover4_Complex_crystal`

Intra-residue LK solvation for atom-pairs beyond torsion-relationship calculated for the bound complex using the crystal structure

### `fa_intra_sol_xover4_Complex_relaxed`

Intra-residue LK solvation for atom-pairs beyond torsion-relationship calculated for the bound complex using the energy-minimized structure

### `fa_intra_sol_xover4_Monomer_crystal`

Intra-residue LK solvation for atom-pairs beyond torsion-relationship calculated for the unbound monomer using the crystal structure

### `fa_intra_sol_xover4_Monomer_relaxed`

Intra-residue LK solvation for atom-pairs beyond torsion-relationship calculated for the unbound monomer using the energy-minimized structure

### `fa_rep_Complex_crystal`

Lennard-Jones repulsive energy between atoms in different residues calculated for the bound complex using the crystal structure

### `fa_rep_Complex_relaxed`

Lennard-Jones repulsive energy between atoms in different residues calculated for the bound complex using the energy-minimized structure

### `fa_rep_Monomer_crystal`

Lennard-Jones repulsive energy between atoms in different residues calculated for the unbound monomer using the crystal structure

### `fa_rep_Monomer_relaxed`

Lennard-Jones repulsive energy between atoms in different residues calculated for the unbound monomer using the energy-minimized structure

### `fa_sol_Complex_crystal`

Lazaridis-Karplus solvation energy calculated for the bound complex using the crystal structure

### `fa_sol_Complex_relaxed`

Lazaridis-Karplus solvation energy calculated for the bound complex using the energy-minimized structure

### `fa_sol_Monomer_crystal`

Lazaridis-Karplus solvation energy calculated for the unbound monomer using the crystal structure

### `fa_sol_Monomer_relaxed`

Lazaridis-Karplus solvation energy calculated for the unbound monomer using the energy-minimized structure

### `hbond_bb_sc_Complex_crystal`

Sidechain-backbone hydrogen bond energy calculated for the bound complex using the crystal structure

### `hbond_bb_sc_Complex_relaxed`

Sidechain-backbone hydrogen bond energy calculated for the bound complex using the energy-minimized structure

### `hbond_bb_sc_Monomer_crystal`

Sidechain-backbone hydrogen bond energy calculated for the unbound monomer using the crystal structure

### `hbond_bb_sc_Monomer_relaxed`

Sidechain-backbone hydrogen bond energy calculated for the unbound monomer using the energy-minimized structure

### `hbond_lr_bb_Complex_crystal`

Backbone-backbone hydrogen bonds distant in primary sequence calculated for the bound complex using the crystal structure

### `hbond_lr_bb_Complex_relaxed`

Backbone-backbone hydrogen bonds distant in primary sequence calculated for the bound complex using the energy-minimized structure

### `hbond_lr_bb_Monomer_crystal`

Backbone-backbone hydrogen bonds distant in primary sequence calculated for the unbound monomer using the crystal structure

### `hbond_lr_bb_Monomer_relaxed`

Backbone-backbone hydrogen bonds distant in primary sequence calculated for the unbound monomer using the energy-minimized structure

### `hbond_sc_Complex_crystal`

Sidechain-sidechain hydrogen bond energy calculated for the bound complex using the crystal structure

### `hbond_sc_Complex_relaxed`

Sidechain-sidechain hydrogen bond energy calculated for the bound complex using the energy-minimized structure

### `hbond_sc_Monomer_crystal`

Sidechain-sidechain hydrogen bond energy calculated for the unbound monomer using the crystal structure

### `hbond_sc_Monomer_relaxed`

Sidechain-sidechain hydrogen bond energy calculated for the unbound monomer using the energy-minimized structure

### `hbond_sr_bb_Complex_crystal`

Backbone-backbone hydrogen bonds close in primary sequence calculated for the bound complex using the crystal structure

### `hbond_sr_bb_Complex_relaxed`

Backbone-backbone hydrogen bonds close in primary sequence calculated for the bound complex using the energy-minimized structure

### `hbond_sr_bb_Monomer_crystal`

Backbone-backbone hydrogen bonds close in primary sequence calculated for the unbound monomer using the crystal structure

### `hbond_sr_bb_Monomer_relaxed`

Backbone-backbone hydrogen bonds close in primary sequence calculated for the unbound monomer using the energy-minimized structure

### `lk_ball_wtd_Complex_crystal`

Weighted sum of anisotropic and isotropic solvation contributions calculated for the bound complex using the crystal structure

### `lk_ball_wtd_Complex_relaxed`

Weighted sum of anisotropic and isotropic solvation contributions calculated for the bound complex using the energy-minimized structure

### `lk_ball_wtd_Monomer_crystal`

Weighted sum of anisotropic and isotropic solvation contributions calculated for the unbound monomer using the crystal structure

### `lk_ball_wtd_Monomer_relaxed`

Weighted sum of anisotropic and isotropic solvation contributions calculated for the unbound monomer using the energy-minimized structure

### `omega_Complex_crystal`

Omega dihedral backbone constraint (planarity ~6° SD) calculated for the bound complex using the crystal structure

### `omega_Complex_relaxed`

Omega dihedral backbone constraint (planarity ~6° SD) calculated for the bound complex using the energy-minimized structure

### `omega_Monomer_crystal`

Omega dihedral backbone constraint (planarity ~6° SD) calculated for the unbound monomer using the crystal structure

### `omega_Monomer_relaxed`

Omega dihedral backbone constraint (planarity ~6° SD) calculated for the unbound monomer using the energy-minimized structure

### `p_aa_pp_Complex_crystal`

Probability of amino acid at Φ/Ψ angles calculated for the bound complex using the crystal structure

### `p_aa_pp_Complex_relaxed`

Probability of amino acid at Φ/Ψ angles calculated for the bound complex using the energy-minimized structure

### `p_aa_pp_Monomer_crystal`

Probability of amino acid at Φ/Ψ angles calculated for the unbound monomer using the crystal structure

### `p_aa_pp_Monomer_relaxed`

Probability of amino acid at Φ/Ψ angles calculated for the unbound monomer using the energy-minimized structure

### `pro_close_Complex_crystal`

Proline ring closure energy and psi angle energy of preceding residue calculated for the bound complex using the crystal structure

### `pro_close_Complex_relaxed`

Proline ring closure energy and psi angle energy of preceding residue calculated for the bound complex using the energy-minimized structure

### `pro_close_Monomer_crystal`

Proline ring closure energy and psi angle energy of preceding residue calculated for the unbound monomer using the crystal structure

### `pro_close_Monomer_relaxed`

Proline ring closure energy and psi angle energy of preceding residue calculated for the unbound monomer using the energy-minimized structure

### `rama_prepro_Complex_crystal`

Backbone torsion preference considering if preceding amino acid is proline calculated for the bound complex using the crystal structure

### `rama_prepro_Complex_relaxed`

Backbone torsion preference considering if preceding amino acid is proline calculated for the bound complex using the energy-minimized structure

### `rama_prepro_Monomer_crystal`

Backbone torsion preference considering if preceding amino acid is proline calculated for the unbound monomer using the crystal structure

### `rama_prepro_Monomer_relaxed`

Backbone torsion preference considering if preceding amino acid is proline calculated for the unbound monomer using the energy-minimized structure

### `ref_Complex_crystal`

Reference energy for amino acid (balances internal energy, important for design) calculated for the bound complex using the crystal structure

### `ref_Complex_relaxed`

Reference energy for amino acid (balances internal energy, important for design) calculated for the bound complex using the energy-minimized structure

### `ref_Monomer_crystal`

Reference energy for amino acid (balances internal energy, important for design) calculated for the unbound monomer using the crystal structure

### `ref_Monomer_relaxed`

Reference energy for amino acid (balances internal energy, important for design) calculated for the unbound monomer using the energy-minimized structure

### `score_Complex_crystal`

Total Rosetta energy score (sum of all energy terms) calculated for the bound complex using the crystal structure

### `score_Complex_relaxed`

Total Rosetta energy score (sum of all energy terms) calculated for the bound complex using the energy-minimized structure

### `score_Monomer_crystal`

Total Rosetta energy score (sum of all energy terms) calculated for the unbound monomer using the crystal structure

### `score_Monomer_relaxed`

Total Rosetta energy score (sum of all energy terms) calculated for the unbound monomer using the energy-minimized structure

### `yhh_planarity_Complex_crystal`

Sidechain hydroxyl group torsion preference (deprecated, see hxl_tors) calculated for the bound complex using the crystal structure

### `yhh_planarity_Complex_relaxed`

Sidechain hydroxyl group torsion preference (deprecated, see hxl_tors) calculated for the bound complex using the energy-minimized structure

### `yhh_planarity_Monomer_crystal`

Sidechain hydroxyl group torsion preference (deprecated, see hxl_tors) calculated for the unbound monomer using the crystal structure

### `yhh_planarity_Monomer_relaxed`

Sidechain hydroxyl group torsion preference (deprecated, see hxl_tors) calculated for the unbound monomer using the energy-minimized structure