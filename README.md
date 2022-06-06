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
## Tips
The running result is provided.