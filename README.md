![image](https://github.com/bm2-lab/PerturBase/assets/24773527/ef604bb0-69e3-4eee-bcbe-59b384b8c8b8)# PerturBase
---------------------------
### Description
We present PerturBase—the first and most comprehensive database designed for the analysis and visualization of scPerturbation data (http://www.perturbase.cn/). PerturBase consolidates 122 datasets from 46 publicly accessible research studies, covering 115 single-modal and 7 multi-modal datasets that include 24254 genetic and 230 chemical perturbations from about 6 million cells.




### DataProcessing_RNAseq
This folder [DataProcessing_RNAseq](DataProcessing_RNAseq) contains the code for preprocessing and denoising the raw RNAseq data.

### DataProcessing_ATACseq
This folder [DataProcessing_ATACseq](DataProcessing_ATACseq) contains the code for preprocessing and denoising the raw ATAC data.

### E-distance
This folder [E-distance](E-distance) contains the code for calculating the distance between the perturb group and control (non-perturbed group)

### Plot
This folder [Plot](Plot) contains the basic plot function.

### Reproduce
This manuscript [reProduce_the_result_of_Mixscape.py](Reproduce/reProduce_the_result_of_Mixscape.py) is the code for reproducing the extended Data of Fig.9a of [Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens](https://www.nature.com/articles/s41588-021-00778-2). And We use this manuscript [evaluate_the_impact_of_hvg.py](Reproduce/evaluate_the_impact_of_hvg.py) to evaluate the impact of the selection of hvg in data preprocessing on the perturbation-specific differentially expressed genes. 

### calDEG
This folder [calDEG](calDEG) contains the code for calculating the differentially expressed genes with the fives softwares。

