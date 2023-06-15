# ComparativeStudyofMMOP
**Paper:** Multimodal Multi-objective Optimization: Comparative Study of the State-of-the-Art

This project collects the sources of the existing multimodal multi-objective evolutionary algorithms and the test suites.

## Usage of ExperimentResult.mat

The average results of IGD, IGDX and 1/PSP are provided in this data file.
For MATLAB, you can use the following command:
``` matlab
ExperimentResult.MMF{i}.igd(runid,numalg)
ExperimentResult.MMF{i}.igdx(runid,numalg)
ExperimentResult.MMF{i}.psp(runid,numalg)
ExperimentResult.IDMP{i}.igdx(runid,numalg)
ExperimentResult.IDMPe{i}.igdx(runid,numalg)
ExperimentResult.POLYGON{i}.igdx(runid,numalg)
```

to get the result of the i-th MMF problem for the runid-th run by numalg-th algorihtm. The algorithms rank is
``` matlab
alglabel={'Omni-optimizer','DN-NSGAII','MO\_Ring\_PSO\_SCD','MO\_PSO\_MM','DNEA','Tri-MOEA&TAR', 'DNEA-L','CPDEA','MP-MMEA','MMOEA/DC','MMEA-WI','HREA'};
```

## Citation information
@article{LI2023101253,
title = {Multimodal multi-objective optimization: Comparative study of the state-of-the-art},
journal = {Swarm and Evolutionary Computation},
volume = {77},
pages = {101253},
year = {2023},
issn = {2210-6502},
doi = {https://doi.org/10.1016/j.swevo.2023.101253},
url = {https://www.sciencedirect.com/science/article/pii/S2210650223000275},
author = {Wenhua Li and Tao Zhang and Rui Wang and Shengjun Huang and Jing Liang},
keywords = {Multimodal multi-objective optimization, Evolutionary computation, Comparative study, Review},
abstract = {Multimodal multi-objective problems (MMOPs) commonly arise in the real world where distant solutions in decision space correspond to very similar objective values. To obtain more Pareto optimal solutions for MMOPs, many multimodal multi-objective evolutionary algorithms (MMEAs) have been proposed. For now, few studies have encompassed most of the representative MMEAs and made a comparative comparison. In this study, we first review the related works during the last two decades. Then, we choose 15 state-of-the-art algorithms that utilize different diversity-maintaining techniques and compared their performance on different types of the existing test suites. Experimental results indicate the strengths and weaknesses of different techniques on different types of MMOPs, thus providing guidance on how to select/design MMEAs in specific scenarios.}
}

## Tips
If you have any problem, please do not hesitate to contact me. Email:liwenhua1030@aliyun.com