# Comparative Study of Multimodal Multi-objective Optimization

[![Paper](https://img.shields.io/badge/DOI-10.1016%2Fj.swevo.2023.101253-blue)](https://doi.org/10.1016/j.swevo.2023.101253)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![MATLAB](https://img.shields.io/badge/MATLAB-PlatEMO-orange)

A collection of **multimodal multi-objective evolutionary algorithms (MMEAs)**, test suites, and full experimental data for the comparative-study paper:

> W. Li, T. Zhang, R. Wang, S. Huang, and J. Liang, "Multimodal multi-objective optimization: Comparative study of the state-of-the-art," *Swarm and Evolutionary Computation*, vol. 77, Article 101253, 2023.
> DOI: [10.1016/j.swevo.2023.101253](https://doi.org/10.1016/j.swevo.2023.101253)

The paper PDFs are included in this repository.

## What's Inside

| Folder | Content |
|--------|---------|
| `PleaEMOCodes/` | Source code of the compared MMEAs (Omni-optimizer, DN-NSGA-II, MO-Ring-PSO-SCD, MO-PSO-MM, DNEA, Tri-MOEA&TAR, DNEA-L, CPDEA, MP-MMEA, MMOEA/DC, MMEA-WI, HREA, ...) |
| `OtherCodes/` | Additional algorithm implementations |
| `Benchmark Problems/` | Test suites: MMF, IDMP, IDMP_e, POLYGON |
| `RefPaper/` | Reference papers of the compared algorithms |
| `ExperimentResult.mat` | Complete raw results (IGD, IGDX, 1/PSP) of all runs |
| `RunResult.xlsx` | Aggregated result tables |

This repository can serve as a **one-stop resource** for running and comparing multimodal multi-objective optimizers on standard benchmarks.

## Usage of ExperimentResult.mat

The average results of IGD, IGDX and 1/PSP are provided in this data file.
For MATLAB, you can use the following command:
```matlab
ExperimentResult.MMF{i}.igd(runid,numalg)
ExperimentResult.MMF{i}.igdx(runid,numalg)
ExperimentResult.MMF{i}.psp(runid,numalg)
ExperimentResult.IDMP{i}.igdx(runid,numalg)
ExperimentResult.IDMPe{i}.igdx(runid,numalg)
ExperimentResult.POLYGON{i}.igdx(runid,numalg)
```

to get the result of the i-th MMF problem for the runid-th run by numalg-th algorithm. The algorithms rank is
```matlab
alglabel={'Omni-optimizer','DN-NSGAII','MO_Ring_PSO_SCD','MO_PSO_MM','DNEA','Tri-MOEA&TAR','DNEA-L','CPDEA','MP-MMEA','MMOEA/DC','MMEA-WI','HREA'};
```

## Citation

```bibtex
@article{Li2023Comparative,
  title   = {Multimodal multi-objective optimization: Comparative study of the state-of-the-art},
  author  = {Li, Wenhua and Zhang, Tao and Wang, Rui and Huang, Shengjun and Liang, Jing},
  journal = {Swarm and Evolutionary Computation},
  volume  = {77},
  pages   = {101253},
  year    = {2023},
  doi     = {10.1016/j.swevo.2023.101253}
}
```

## Tips

If you have any problem, please do not hesitate to contact me. Email: liwenhua1030@aliyun.com
