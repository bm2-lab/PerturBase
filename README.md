## PerturBase

---

### Description

We present PerturBase—the most comprehensive database designed for the analysis and visualization of single-cell perturbation datasets (http://www.perturbase.cn/). PerturBase consolidates 122 datasets from 46 publicly accessible research studies, covering 115 single-modal and 7 multi-modal datasets that include 24254 genetic and 230 chemical perturbations from about 5 million cells.

We present a streamlined analysis workflow and provide demo data to help users better understand our analysis.

**Please ensure this GitHub repository is fully cloned to your local server**

### Demo Data

Due to the size exceeding GitHub's limitations, demo data is hosted at [Figshare](https://figshare.com/s/dddc4ddf91d0b100fd6c). If you wish to replicate our results, please retain the raw data files. Specifically,

For RNA, retain:

* raw.h5ad

For ATAC, retain:

* fragment.tsv.gz
* fragment.tsv.gz.tbi
* cbc_perturbation_dict.tsv

Following the instructions, the "Demo_data" folder would look like:

```bash
Demo_data
├── RNA
│   ├── raw.h5ad
│   
│   
├── ATAC
    ├── fragment.tsv.gz
    ├── fragment.tsv.gz.tbi
    ├── cbc_perturbation_dict.tsv

```

### Environment

R_env.yaml, Python_env.yaml are the `R_env` and `Python_env` packages used in the `pipeline_bash_files`. Please install them in your environment.

### MsigDB

Molecular Signatures Database used in PerturBase

### pipeline_bash_files

This is whole pipeline script. For each data set, select the corresponding script, input the data path and species, and complete the entire analysis.

For our demo ATAC fragment data set, run:

```bash
bash ATAC_fragment_all.sh /path_to_your_file/Demo_data/ATAC Hs
```

For our demo RNA data set, run:

```bash
bash RNA_all.sh  /path_to_your_file/Demo_data/RNA  Hs
```


To enhance user experience, we reorganized the code by replacing the absolute paths from the development machine with user-definable paths. Before executing the script, please replace `Python_env`, `R_env`, `Script_base_path`, `msigdb_Signature_path`, etc., with the paths on your local machine, and source the `conda.sh` under your Anaconda installation.

### scripts

Scripts that is invoked by the main program, i.e, pipeline_bash_files.

### Reproducibility

This manuscript [reproduce_the_result_of_mixscape.py](Reproducibility/reproduce_the_result_of_mixscape.py) is the code for reproducing the extended Data of Fig.9a of [Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens](https://www.nature.com/articles/s41588-021-00778-2). And We use the manuscript [evaluate_the_impact_of_hvg.py](Reproducibility/evaluate_the_impact_of_hvg.py) to evaluate the impact of the selection of hvg in data preprocessing on the perturbation-specific differentially expressed genes.

### Bug Report

Any question is welcomed in the [issue](https://github.com/bm2-lab/PerturBase/issues).

## Citation

Zhiting Wei, Duanmiao Si, Qi Liu et al. *PerturBase: a comprehensive database for single-cell perturbation data analysis and visualization*, 2024.

## Contacts

bm2-lab@tongji.edu.cn
1810546@tongji.edu.cn
2311470@tongji.edu.cn
