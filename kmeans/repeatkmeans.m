clear all; close all;clc
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\')

zdim = 0; %Taly: from bash
distanceMethod = 'correlation'; %Taly: from bash
nreps = 1; %Taly: from bash
repkmeans = 'repkmeans';

[test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path, PAIRED_RUNS] = setGlobalVar();

[ts_file_list, ordered_files_count] = get_files_list(basedir, tokens, PAIRED_RUNS);


masterdir = [basedir, test];
savedir = fullfile(masterdir, repkmeans);
if ~(exist(savedir, 'dir') == 7)
    mkdir(savedir);
end
cd(savedir);

[concTS, scan_name,SCAN_LEN]  = get_data(ts_file_list);        
save('kmeans_meta.mat', 'ts_file_list', 'SCAN_LEN', 'scan_name', 'nsplits', 'ordered_files_count');
save('concTS.mat', 'concTS')

disp('data loaded');
for numClusters = 1:maxNumClusters %Taly - instead of bash files %189 frames 13^2 = 169
    for split = 1:nsplits     
        disp(['K = ',num2str(numClusters),'Split = ',num2str(split)]);
        savedir = fullfile(masterdir,repkmeans);
        cd(savedir);

        savedir = ['kmeans',num2str(split),distanceMethod,'k_',num2str(numClusters)];
        if ~(exist(savedir, 'dir') == 7)
            mkdir(savedir);
        end
        %%
        cd(savedir);

        disp('start k-means');
        for R = 1:nreps
            [partition,~,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod, 'MaxIter', 1000);
            save(['kmeans',num2str(split),'k_',num2str(numClusters),'rep',num2str(R),'.mat'],'partition','sumd');
            clear partition
            clear sumd
            disp(['K-means ',num2str(R)]);
        end
        
        disp('complete');
    end
end

