%% make radar plots
clear all; close all;clc
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master_1\code\')
addpath('C:\Users\marko\ketamine\energy_landscape-main\misc_code')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\assesscluster')
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master_1\code\miscfxns')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\plottingfxns')

[test, basedir, tokens, atlas, nparc] = param();
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = 6;

savedir = fullfile(basedir, test, 'Energy');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp','_k',num2str(numClusters),'.mat']))


%%

overallNames = clusterNames;
[nparc,numClusters] = size(centroids);
[~,~,~,net7angle] = NAME_CLUSTERS_ANGLE(centroids, atlas);

%% plot

YeoColors = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;];
YeoColors = [YeoColors;YeoColors];

[~,~,net7angle_Up,net7angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids, atlas);
YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
numNets = numel(YeoNetNames);


%% make radial plots

clusterColors = GET_CLUSTER_COLORS(numClusters);

clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; thetaNames{8} = '';

f=figure;
ordedred_k = [1, 6, 2, 5, 3, 4];
for K_ = 1:numClusters
    K = ordedred_k(K_);
    ax = subplot(1,numClusters,K_,polaraxes); hold on
    polarplot(netAngle,[net7angle_Up(K,:) net7angle_Up(K,1)],'k');
    polarplot(netAngle,[net7angle_Down(K,:) net7angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rlim([0.0 0.7]);
    rticks([0.2 0.4 0.8]); rticklabels({'','0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K_,:),'FontSize',12);
end
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];
