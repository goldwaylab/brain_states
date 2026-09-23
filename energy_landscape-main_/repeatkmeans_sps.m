%% this data goes into elbow_sps.m and is used to choose k
clear all; close all;clc

[test, basedir, tokens, atlas, nparc] = param();
cd(basedir);

%% load data
[ts_file_list, ordered_files_count] = get_files_list(basedir, tokens);
[concTS, scan_name, SCAN_LEN]  = get_data(ts_file_list);

savedir = fullfile(basedir, test, 'Energy\');
if ~exist(savedir, 'dir')
    mkdir(savedir); 
end
cd(savedir);
save([savedir, 'concTS.mat'], 'concTS');
save([savedir, 'concTS_meta.mat'], 'scan_name', 'SCAN_LEN');

%% set inputs

nreps = 50; % how many times to repeat clustering. will choose lowest error solution
distanceMethod = 'correlation';
maxI = 1000; % how many times you allow kmeans to try to converge
TR = min(SCAN_LEN);
maxk = round(sqrt(TR)) -1; %k^2 must be less than TR to capture all transitions

for numClusters = 2:maxk %because we have 220 frames, k^2 must be less than 220 to capture all transitions
    
    disp(['K = ',num2str(numClusters)]);

    disp('start k-means');
    
    [partition,~,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    save(fullfile(savedir,['kmeans','k_',num2str(numClusters),'.mat']),'partition','sumd')
   
    clear partition
    clear sumd
    
    
end

disp('complete');
