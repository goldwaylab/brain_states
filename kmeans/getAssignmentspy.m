%addpaths; %Taly
clear all;
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\')
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\miscfxns')


[test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path] = setGlobalVar();
maxNumClusters = 10; %REMOVE TEMP REMOVE TEMP REMOVE TEMP REMOVE TEMP REMOVE TEMP 
%Taly: bash file parameters 
opts.VariableNamingRule = 'preserve';
distanceMethod = 'correlation';
clusterRange = minNumClusters:maxNumClusters; numK = 1 + maxNumClusters-minNumClusters;
k_offset = minNumClusters - 1;

masterdir = [basedir, test];

%Taly load data
repkmeans = '\repkmeans';
load([masterdir, repkmeans, '\kmeans_meta.mat']) %'ts_file_list', 'SCAN_LEN', 'scan_name', 'nsplits'

savedir = [masterdir,'\clusterAssignments'];
if ~(exist(savedir, 'dir') == 7)
    mkdir(savedir);
end

% for some reason, k-means jobs will sometimes randomly die
% check which saved k-means solution files exist, store in logical matrix splits_by_K
% splits are just repetitions of clustering with full sample, not subsamples
d = [masterdir,repkmeans];
splits_by_K = zeros(nsplits,numK);
for S = 1:nsplits
  for K = 1:numK
      fname = [d,'\kmeans',num2str(S),distanceMethod,'k_',num2str(clusterRange(K)),'\kmeans',...
          num2str(S),'k_',num2str(clusterRange(K)),'rep1.mat'];
      splits_by_K(S,K) = logical(exist(fname,'file'));
  end
end

nsplits = min(sum(splits_by_K,1));

%Taly
load([masterdir, repkmeans, '\concTS.mat'])
totalNumTPs = length(concTS);    

combPartitions = cell(numK,1); combSumD = cell(numK,1);
disp('start loading k means partitions');
for K = minNumClusters:maxNumClusters
    existingSplits = find(splits_by_K(:,K-k_offset));
    combPartitions{K-k_offset} = int8(zeros(totalNumTPs,nsplits));
    combSumD{K-k_offset} = zeros(1,nsplits);
        for splitInd = 1:nsplits
            S = existingSplits(splitInd);
            startInd = sum(splits_by_K(1:S,K-k_offset)) - splits_by_K(S,K-k_offset);
            load([d,'/kmeans',num2str(S),distanceMethod,'k_',num2str(K),'/kmeans',num2str(S),'k_',num2str(K),'rep1.mat']);
                combPartitions{K-k_offset}(:,splitInd) = int8(partition);
            %Taly original code reverse eng.: sumd comes from repeatkmeans and it is the output of k-means. 
            %It is the sum of squared distances (SSD) between data points and cluster centroids, which is a common way to measure the quality of a clustering solution.
            combSumD{K-k_offset}(:,splitInd) = sum(sumd);
        end
    %Taly: cannnot find the definition of multislice_pair_labeling
    %combPartitions{K-k_offset} = int8(multislice_pair_labeling(combPartitions{K-k_offset}));
    disp(['load partitions: K=',num2str(K)]);
end

nreps = nsplits;
save(fullfile(savedir,['repkmeansPartitions',distanceMethod,'.mat']),'combPartitions','combSumD','nreps','splits_by_K','-v7.3');

%concTS =
%dlmread(fullfile(basedir,['data/ConcTSCSV_','.csv']),','); Taly 
disp('time series loaded');

disp('finding lowest MSE partition');
for numClusters = minNumClusters:maxNumClusters
    combBestClusterCentroids = zeros(nparc,numClusters);
    %Taly original code reverse eng.: choose the cluster with the smallest sumd, i.e the best clustering 
    combBestClusterInd = find(combSumD{numClusters-k_offset} == min(combSumD{numClusters-k_offset}));
    combBestClusterInd = combBestClusterInd(1); %if multiple partitions exactly the same, just use the first one
    % store partition with lowest error
    clusterAssignments.(['k',num2str(numClusters)]).partition = combPartitions{numClusters - k_offset}(:,combBestClusterInd);
    clusterAssignments.(['k',num2str(numClusters)]).scan_name = scan_name;
    clusterAssignments.(['k',num2str(numClusters)]).SCAN_LEN = SCAN_LEN;
    

    % compute centroids based on that partition, averaging across all data used for clustering
    for K = 1:numClusters
        %Taly original code reverse eng.: take previousely chosen best
        %result of kmeans and find centroid for each cluster by mean on all
        %volumes assigned to that cluster
        combBestClusterCentroids(:,K) = mean(concTS(clusterAssignments.(['k',num2str(numClusters)]).partition == K,:),1)';
    end
    
    clusterAssignments.(['k',num2str(numClusters)]).bestCentroid = combBestClusterCentroids;
    clusterAssignments.(['k',num2str(numClusters)]).bestClusterInd = combBestClusterInd;
        
    clusterAssignments.(['k',num2str(numClusters)]).sumD = combSumD{K-k_offset}(combBestClusterInd);
    clusterAssignments.(['k',num2str(numClusters)]).clusterNames = NAME_CLUSTERS_ANGLE(combBestClusterCentroids, atlas);
    [clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(combBestClusterCentroids, atlas);
    clusterAssignments.(['k',num2str(numClusters)]).clusterNamesUp = clusterNamesUp;
    clusterAssignments.(['k',num2str(numClusters)]).clusterNamesDown = clusterNamesDown;

    save(fullfile(savedir,['k',num2str(numClusters),'.mat']),'clusterAssignments');
    clear clusterAssignments
    disp(['K=',num2str(K)]);
end
