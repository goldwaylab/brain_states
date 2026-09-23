# brain_states

Analysis code for **Brain State Dynamics in Ketamine-Induced Dissociation Resemble Those in Posttraumatic Stress Disorder**
(Goldway et al., *Biological Psychiatry Global Open Science*, 2026; [DOI](https://doi.org/10.1016/j.bpsgos.2025.100655), [OSF project](https://osf.io/46zfm/)).

## Overview

The study asks whether brain-state dynamics linked to dissociation under intravenous ketamine resemble those in PTSD.
The code clusters fMRI time series into recurring brain states, quantifies state dynamics
(fractional occupancy, recurrence rate, dwell time, transition probabilities, and transition energy from network control theory),
and relates them to dissociation (CADSS) in healthy volunteers given ketamine and in patients with PTSD.

## Pipeline

**MATLAB** (clustering and dynamics):

1. `repeatkmeans.m`: k-means clustering; creates a `repkmeans` directory with the k-means output folders, `kmeans_meta.mat` and `concTS.mat`.
2. `getAssignmentspy.m`: creates the `clusterAssignments` folder with the best partitioning for each number of clusters *k*.
3. `elbow.m`: plots to help choose the number of clusters (states).
4. `systems_plot.m`: radial plots of the state centroids, saved in `analyses\centroids`.
5. `getDynamics.m`: fractional occupancy, recurrence rate and dwell time, saved in `analyses\transitionprobabilities`.
6. `dynamics_tables.m`: tables of the above under `analyses`.

**R** (statistics and figures): `Brain_states_main.R`. Set `root_dir` at the top of the script to the folder that holds the data files
(the default, `"."`, assumes you run it from the repository root).

## Data included

| File | Contents |
|---|---|
| `CADSS.csv` | Clinician-Administered Dissociative States Scale scores (ketamine study) |
| `demog_k.csv`, `demog_ptsd.csv` | Demographics for the ketamine and PTSD samples |
| `energy.csv` | Transition energies between brain states |
| `fo_results.xlsx` | Fractional-occupancy results |

## Requirements

- MATLAB for the clustering and dynamics scripts.
- R with: `readxl`, `tidyr`, `dplyr`, `afex`, `emmeans`, `Hmisc`, `ggplot2`, `pheatmap`, `corrplot`, `ggstatsplot`, `gridExtra`, `stringr`, `R.matlab`, `RColorBrewer`, `viridis`, `reshape2`, `lm.beta`, `plotrix`, `gsubfn`.

## Acknowledgements

Based on [ejcorn/brain_states](https://github.com/ejcorn/brain_states). Bundled third-party code (NIfTI tools, Colormaps, `calc_lz_complexity`, `violin`, energy-landscape code) keeps its own license, included in the respective folders.

## Citation

Goldway, N., Markovits, T., Fine, N., Fruchtman-Steinbok, T., Gurevitch, G., Deco, G., Sharon, H., & Hendler, T. (2026). Brain State Dynamics in Ketamine-Induced Dissociation Resemble Those in Posttraumatic Stress Disorder. *Biological Psychiatry Global Open Science, 6*(2), 100655. https://doi.org/10.1016/j.bpsgos.2025.100655
