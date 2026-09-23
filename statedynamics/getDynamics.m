clear all
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\dynamicsfxns');
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code')

% TR_CONST - Same TR for all scans, rather than reading *.json
[test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path] = setGlobalVar();
numClusters = 6;

masterdir = [basedir, test];
load([masterdir, '\repkmeans\kmeans_meta.mat']) %'ts_file_list', 'SCAN_LEN', 'scan_name', 'nsplits', 'ses_1', 'ses_2'
load(fullfile(masterdir,['clusterAssignments/repkmeansPartitionscorrelation','.mat']));
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),'.mat']));
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;    
[list_of_files, all_TRs, numTRs_all, subjInd, scanInd, all_nobs, improved] = get_variables(test, scanlab, TR, masterdir, TR_path);
    
savedir = [masterdir,'/analyses/transitionprobabilities'];
if ~(exist(savedir, 'dir') == 7)
    mkdir(savedir);
end    
save(fullfile(savedir,'list_of_files.mat'),'list_of_files')    
save(fullfile(savedir,'scanlab.mat'),'scanlab')    

for i = 1:numel(scanlab)
    nobs = all_nobs(i);
    %% calculate transition prob

    [transitionProbability,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(kClusterAssignments(scanInd == (i-1)),subjInd(scanInd == (i-1)));    
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilityMatrices_k',num2str(numClusters),'.mat']),'transitionProbabilityMats')
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilities_k',num2str(numClusters),'.mat']),'transitionProbability')
    save(fullfile(savedir,[scanlab{i},'NumTransitions_k',num2str(numClusters),'.mat']),'numTransitions')

    [transitionProbability,transitionProbabilityMats] = GET_TRANS_PROBS_NO_PERSIST(kClusterAssignments(scanInd == (i-1)),subjInd(scanInd == (i-1)));    
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilitiesNoPersist_k',num2str(numClusters),'.mat']),'transitionProbability','transitionProbabilityMats');

    %% calculate fractional occupancy, i.e. percentage of time spent in each state, regardless of order
    
    FractionalOccupancy = zeros(nobs,numClusters);

    for N = 1:nobs
        for K = 1:numClusters    
            numTRs = numTRs_all{i, N};
            FractionalOccupancy(N,K) = sum(and(and(kClusterAssignments == K,subjInd' == N),scanInd' == (i-1))) / numTRs;%Taly: scanInd-->scanInd', numTRs(i)--> numTRs
        end
    end

    save(fullfile(savedir,[scanlab{i},'FractionalOccupancy_k',num2str(numClusters),'.mat']),'FractionalOccupancy')

    %% calculate dwell time, i.e. average number of subsequent TRs that each state lasts for
    % store both mean and median, because dwell time may not be normally distributed
    DwellTimeMean = zeros(nobs,numClusters);
    DwellTimeMedian = zeros(nobs,numClusters);
    RunRate = zeros(nobs,numClusters);
    for N = 1:nobs
        TR = all_TRs{i, N};
        [dt_mean,dt_median,~,~,n_runs] = CALC_DWELL_TIME(kClusterAssignments(subjInd' == N & scanInd' == (i-1)),numClusters); %Taly scanInd --> scanInd'
        DwellTimeMean(N,:) = dt_mean*TR;        % store dwell time in seconds
        DwellTimeMedian(N,:) = dt_median*TR;        % store dwell time in seconds
        % store rate of appearance of runs, i.e. how many DMN runs appear in 1 minute
        % first calculate runs/TR,then runs/sec, then runs/min
        numTRs = numTRs_all{i, N};
        RunRate(N,:) = 60*(1/TR)*(n_runs/numTRs);  
    end

    save(fullfile(savedir,[scanlab{i},'DwellTime_k',num2str(numClusters),'.mat']),'DwellTimeMean','DwellTimeMedian','RunRate')

end

