%% step 1 after choosing k (or a range of k)

%%
clear all; close all;clc
addpath('C:\Users\marko\ketamine\energy_landscape-main\misc_code')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\assesscluster')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\plottingfxns')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\miscfxns')

[test, basedir, tokens, atlas, nparc] = param();
cd(basedir);

%% set inputs

savedir = fullfile(basedir, test, 'Energy');mkdir(savedir);		% set save directory
load([savedir, '\concTS.mat']);
load([savedir, '\concTS_meta.mat']);

distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution
maxI = 1000; % how many times to let kmeans try to converge before aborting rep


[T,nparc] = size(concTS);


%% generate N partitions that will be compared pair-wise for mutual information

for numClusters=[6] %if undecided, run loop. If confident about k, you can run over only 1 value for numClusters
    disp(['Starting clusters number: ',num2str(numClusters)]);
    N = 10; %number of partitions to create and compare
    parts = NaN(T,N); %partitions will be stored here
    D = NaN(N,T,numClusters); %distance matrices will be stored here
    
    for i=1:N
        disp(['Clusters: ',num2str(numClusters),'. Starting kmeans number: ',num2str(i)]);
        [parts(:,i),~,~,D(i,:,:)] = kmeans(concTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    end
    
    %% calculate adjusted mutual information for every pair of partitions
    
    ami_results = NaN(N,N);
    
    for i=1:N
        for j=1:N
            ami_results(i,j) = ami(parts(:,i),parts(:,j));
        end
    end
    
    % assess
    [m,ind] = max(sum(ami_results,1)); %ind corresponds to the partition which has the highest mutual information with all other partitions
    partition = parts(:,ind); % take partition that has most agreement with all other for further analysis
    
    % plot SI Figure 1
    f = figure;
    
    imagesc(ami_results); title(['Adjusted Mutal Information between Partitions k=',num2str(numClusters)]); colorbar;
    axis square; set(gca,'FontSize',8);
    f.PaperUnits = 'inches';
    f.PaperSize = [4 2];
    f.PaperPosition = [0 0 4 2];
    saveas(f,fullfile(savedir,['AMI_bp','_k',num2str(numClusters),'.pdf']));
    
    
    %% compute centroids and plot
    
    centroids = GET_CENTROIDS(concTS,partition,numClusters);
    % name clusters based on alignment with Yeo resting state networks
    % centroids = zscore(centroids);
    clusterNames = NAME_CLUSTERS_ANGLE(centroids, atlas);  % need to add a prior Yeo partition labels for your parcellation
    [clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids, atlas);  % need to add a prior Yeo partition labels for your parcellation
    
    %SI Figure 3
    f = figure;
    subplot(1,2,1); imagesc(centroids); title('Centroids'); xticks(1:numClusters); xticklabels(clusterNames);
    axis square; colorbar; set(gca,'FontSize',8); COLOR_TICK_LABELS(true,false,numClusters);
    subplot(1,2,2); imagesc(corr(centroids)); title('Centroid Similarity'); colorbar; caxis([-1 1]);
    axis square; set(gca,'FontSize',8); xticks(1:numClusters); yticks(1:numClusters);
    xticklabels(clusterNames); yticklabels(clusterNames); xtickangle(90);
    COLOR_TICK_LABELS(true,true,numClusters);
    f.PaperUnits = 'inches';
    f.PaperSize = [4 2];
    f.PaperPosition = [0 0 4 2];
    saveas(f,fullfile(savedir,['Centroids_bp','_k',num2str(numClusters),'.pdf']));
    
    %% save
    
    save(fullfile(savedir,['Partition_bp','_k',num2str(numClusters),'.mat']), 'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids');
end
