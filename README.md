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
- Monthly updates

> Before using this repository, you must **manually download the ANABAG dataset** (see below).

---

## Step 1: Download the ANABAG Dataset

1. Download from the following link:  
   **[https://zenodo.org/records/15794632](<!-- TODO: insert link -->)**

Last update: 24/06/2025 

2. Extract the `.tar` archive:
   ```bash
   tar -xvf data.tar
    ````

3. Move the extracted folder into this project directory (i.e., where `README.md` is located). The directory must be named 'data'.

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

### References:

* Example config file: `dataset_info/example_configuration`
* All possible parameters: `dataset_info/complete_dictionnary_of_features.txt`
* Explanation of parameters: `dataset_info/parameters_dictionnary.md`

---

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

---

## 🗂️ Project Structure

```
ANABAG-handler/
├── src/
│   ├── select_complexes.py              # Main selection script
│   └── quick_analysis_example.ipynb     # Optional notebook for visualization
├── dataset_info/
│   ├── selection_file_complete.tsv
│   ├── cluster_informations.tsv
│   ├── per_chain_pdbff_informations.tsv
│   ├── method_resolution.tsv
│   ├── sequences_initial_chain.tsv
│   └── sequences_formated_chain.txt
├── images/
│   └── 3ulu_publi.png                   # Example visual / schema
├── README.md
└── (Place extracted ANABAG dataset here)
```

---
