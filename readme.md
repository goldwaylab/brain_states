#### Based on https://github.com/ejcorn/brain_states.git

#### 1. repeatkmeans.m - clustering, creates repkmeans directory with kmeans output folders and kmeans_meta.mat, concTS.mat
#### 2. getAssignmentspy.m - gets kmeans and creates clusterAssignments folder with the best partitioning for each number of k
#### 3. elbow.m - create assisting graphs to decide on the number of clusters/states
#### 4. systems_plot.m - gets clusterAssignments and creates the radial plots in analyses\centroids
#### 5. getDynamics.m - gets clusterAssignments and repkmeans and creates FO, RR and Dwell time results in analyses\transitionprobabilities
#### 6. dynamics_tables.m - creates FO, RR and Dwell time tables under analyses

## Citation

Goldway, N., Markovits, T., Fine, N., Fruchtman-Steinbok, T., Gurevitch, G., Deco, G., Sharon, H., & Hendler, T. (2026). Brain State Dynamics in Ketamine-Induced Dissociation Resemble Those in Posttraumatic Stress Disorder. *Biological Psychiatry Global Open Science, 6*(2), 100655. https://doi.org/10.1016/j.bpsgos.2025.100655

This repository is a snapshot of the code shared on the paper's OSF project (https://osf.io/46zfm/).
