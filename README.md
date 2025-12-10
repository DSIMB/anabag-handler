---
---
# ANABAG Handler

![Graphical presentation of ANABAG](images/banner_01.png)

This repository provides Python scripts to **filter and extract specific antibody-antigen complexes** and associated features from the **ANABAG dataset**.

---

## What Is ANABAG?

**ANABAG** (ANnotated AntiBody AntiGen) is a curated dataset of antibody–antigen complexes. It includes:

- 3D structural data (with various formats)
- Per-sequence and per-residue features
- Frequent updates (the zenodo record will be updated when a substantial number of new cases are added)

> Before using this repository, you must **manually download the ANABAG dataset** (see below). You can also build a subset of biological units without downloading the Zenodo dataset. If so, please see section Selecting and Extracting Biological Units Without a Pre-Downloaded Dataset

---

## Step 1: Download the ANABAG Dataset

1. Download from the following link:  
   **[https://zenodo.org/records/17065788](<!-- TODO: insert link -->)**

Last update: 29/08/2025 (Note 10/12/25: I had to defend my PhD, so I was not updating ANABAG. We are working now on an automatic solution for updates. Coming soon) 
You can download the data.tar.gz, or the light_version.tar.gz, both are compatible with the python scripts.  
You also need the per_residue_files if you wish to get the per residue features. 
In the ./dataset_info you can find per-BU and per-chain informations.

2. Extract the `.tar` archive:
   ```bash
   tar -xvf data.tar
   tar -xvf per_residue_files.tar
    ````

3. Move the extracted data or light version folder into this project directory (i.e., where `README.md` is located). The directory must be named 'data'.

4. Move the per_residue_information .tsv files (both) to the dataset_info/ directory.

---

## Step 2: Set Up the Environment

You need **Python 3.10+** and a few Python packages.

If using **Conda**:

```bash
conda create -n anabag_env python=3.10 pandas jupyter numpy matplotlib seaborn
conda activate anabag_env
```

---

## Step 3: Select Complexes or Features

Use the main script `select_complexes.py` to select complexes or features based on your criteria.

### Example usage:

```bash
python src/select_complexes.py path/to/ANABAG path/to/your_config.config
```

Example with provided template:

```bash
python src/select_complexes.py ./ dataset_info/selection_file_complete.config
```

📁 Output will be saved in the `/my_dataset/` directory:

* `/my_dataset/structures`: contains selected structures
* `/my_dataset/files`: contains selected feature files

---

## Create a Configuration File

The configuration file defines **how to filter complexes**. It is structured in four sections:

### Sections:

* `Parameters for: Antigen`
* `Parameters for: Antibody`
* `Parameters for: Complex`
* `Parameters for: Selection` (controls what is extracted)

### Syntax:

```ini
Parameters for: Antigen
SequenceIdentity = SG95AG          # SG20AG, SG40AG, SG60AG, SG80AG, SG95AG, SG100AG
UA_Active_site = 0,8               # Range (min, max)

Parameters for: Selection
per_residue_info = True            # Extract per-residue feature files
formatted_structures = True       # Extract formatted structures 
initial_structures = False        # Extract original chain label structures
rosetta_structures = False        # Extract Rosetta-relaxed structures
hetatm_structures = False         # Include hetero atoms
```

Note that if you use the light version - containing only the "formatted_structures" you should set initial_structures, rosetta_structures and hetatm_structures flags to False.

### Reference files:

* All possible parameters: `dataset_info/selection_file_complete.config`
* Example configuration: `dataset_info/selection_example.config`
* Explanation of features: `dataset_info/features.md`

---

## Selecting and Extracting Biological Units Without a Pre-Downloaded Dataset

If you don’t already have the dataset downloaded, you can still build your own subset of biological units using the scripts provided and the files in dataset_info. This workflow lets you select complexes according to your criteria, fetch them directly from the Protein Data Bank (PDB), and extract + format the relevant biological units. You won't have access to the modelled structures, the rosetta structures and the features computed on modelled missing regions.

The process involves three scripts:

### 1. Select complexes

Use select_complexes.py with a configuration file that specifies your selection parameters.
Make sure to set (in the .config file):

```ini
Parameters for: Selection
build_from_pdb = True
```

This option tells the script to generate a fetcher file containing:

PDB IDs, Chains to extract, Unique biological unit names (One_digit_id)

An example configuration can be found at:

dataset_info/selection_fetcher.config

### 2. Download PDB files

Once the fetcher file is generated, you can either:

Download the PDBs manually, or

Use the helper script fetch_pdbs.py:

```bash
python src/fetch_pdbs.py ./dataset_info ./pdb_downloads
```

This will download all required structures into the pdb_downloads/ directory.

### 3. Extract and format biological units (BU)

After downloading the raw PDB files, run extract_bus.py to parse and format them.
This script will:

Format the BU with its initial chains (renumber residues from 1->N, reorganize the order to: antigen, antibody) [ “initial_chains” version]
Format the BU in the format chain A (antigen) and B (antibody), [ “formatted_chains” version]

Generate per-residue files for each selected biological unit.

Organize outputs in the same structure as data/:

```ini
digit_id/structure/ # the .pdbs formated corresponding to the BU
digit_id/files/ # the per residue feature files corresponding to the residues in the extracted BUs 
```

Example:

```bash
python src/extract_bus.py ./dataset_info ./pdb_downloads ./subset_data
```

### ⚠️ Important Notes

If the extracted BUs contain modelled residues in the Zenodo version, these residues will be ignored in this workflow.

You can control the inclusion/exclusion of these structures with the Number_of_modelled_residues parameter in your config file.

If Number_of_modelled_residues is not set to "0,0" then you are selecting BUs that have modelled residues in the Zenodo version.

These residues are absent from the structure you just extracted, and will be ignored in the per residue files.

However, since the features have been calculated on the modelled structure, you may observe differences in features computed on the overall structure (e.g., net charge, percentage of secondary structures, etc...).

### Polyspecificity and Epitope Convergence Data

Direct access to the polyspecific antibodies and epitope convergence cases identified in the ANABAG paper is provided through two dedicated files:
Polyspecific Antibodies
Access polyspecific antibody cases via ./dataset_info/polyspecificity_labels.csv. This file contains antibodies that bind multiple distinct antigens, organized by different sequence identity thresholds.
Epitope Convergence
Access epitope convergence cases via ./dataset_info/epitope_convergence.csv. This file identifies different antibodies that target similar epitope regions.
File Structure

Each column represents a specific condition or grouping criterion
Values in cells represent group labels assigned to biological units (BUs)
Label -1 indicates no group assignment
BUs sharing the same label (other than -1) belong to the same polyspecific or convergent group

Usage Example
To identify polyspecific antibodies:

Select the appropriate column based on your sequence identity threshold of interest
Find all BUs with the same label (excluding -1)
These BUs represent the same antibody binding to different antigens

For epitope convergence, the same principle applies: BUs with matching labels represent different antibodies targeting similar epitopes.

## 📊 Visualize the Data (Optional)

You can preview and analyze selected data using the provided Jupyter notebook.

1. Start Jupyter:

   ```bash
   jupyter notebook
   ```

2. Open: `src/quick_analysis_example.ipynb`

3. Set your dataset path inside the notebook:

   ```python
   path_to_mydataset = 'path/to/anabag-handler/my_dataset/files'
   ```

### Reference & citation
ANABAG: Annotated Antibody–Antigen Data Set with Unique Features for Antibody Engineering Applications
Grandguillaume, Ilyas
Barroso da Silva, Fernando Luís
Etchebest, Catherine
Journal of Chemical Information and Modeling
doi: 10.1021/acs.jcim.5c01599
