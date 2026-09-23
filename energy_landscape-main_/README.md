#### 1. based https://github.com/singlesp/energy_landscape/blob/main/README.md
#### 2. repeatkmeans_sps.m - clustering
#### 3. repeatkmeans_sps.m -  create assisting graphs to decide on the number of clusters/states
#### 4. ami_calc.m - assess clustering stability to choose the partition with the highest amount of adjusted mutual information shared with all other partitions
#### 5. subcentroids.m - generate subject-specific centroids for energy calculations.
#### 6. T_sweep_sps.m - sweep T and calculate energies then correlate with the emprical tansistion probabilities. Choose T that maximizes the negative correlation.
#### 7. subj_energy.m - create energy graphs
