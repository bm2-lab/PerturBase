### PerturBase
---------------------------
### Description
We present PerturBase—the most comprehensive database designed for the analysis and visualization of scPerturbation data (http://www.perturbase.cn/). PerturBase consolidates 122 datasets from 46 publicly accessible research studies, covering 115 single-modal and 7 multi-modal datasets that include 24254 genetic and 230 chemical perturbations from about 5 million cells.




### DataProcessing_RNAseq
This folder [DataProcessing_RNAseq](DataProcessing_RNAseq) contains the code for preprocessing and denoising the raw RNAseq data.

### DataProcessing_ATACseq
This folder [DataProcessing_ATACseq](DataProcessing_ATACseq) contains the code for preprocessing and denoising the raw ATAC data including Fragment tsv file and Peak Count matrix h5ad file.

### E-distance
This folder [E-distance](E-distance) contains the code for calculating the distance between the perturb group and control group.

### Enrichment
This folder [Enrichment](Enrichment) contains the code for GO and KEGG enrichment and visualization.

### Reproduce
This manuscript [reProduce_the_result_of_Mixscape.py](Reproduce/reProduce_the_result_of_Mixscape.py) is the code for reproducing the extended Data of Fig.9a of [Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens](https://www.nature.com/articles/s41588-021-00778-2). And We use this manuscript [evaluate_the_impact_of_hvg.py](Reproduce/evaluate_the_impact_of_hvg.py) to evaluate the impact of the selection of hvg in data preprocessing on the perturbation-specific differentially expressed genes. 

### calDEG
This folder [calDEG](calDEG) contains the code for calculating the differentially expressed genes with the fives softwares。

## Citation
Zhiting Wei, Duanmiao Si, Qi Liu et al. *PerturBase: a comprehensive database for single-cell perturbation data analysis and visualization*, 2024.
## Contacts
bm2-lab@tongji.edu.cn
1810546@tongji.edu.cn
2311470@tongji.edu.cn
