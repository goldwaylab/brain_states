%  make elbow plot (Fig S2a-b) of within cluster / total variance explained by clusters at different k values
% assumes correlation distance

%% init

%% SI Figure 13 (need to run repeatkmeans_sps.m first)

clear all; close all;clc
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master_1\code\')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\assesscluster')

[test, basedir, tokens, atlas, nparc] = param();
savedir = fullfile(basedir,test,'Energy\');mkdir(savedir);

load([savedir, 'concTS.mat']);
load([savedir, 'concTS_meta.mat']);
cd(basedir);


%% inputs

N = size(concTS,1); % number of observations
TR = min(SCAN_LEN);
maxk = round(sqrt(TR))-1; %<- check, may not need "-1" maxk^2 should be < total TR in scan

k_rng = 2:maxk;

VarianceExplained = zeros(length(k_rng),1);

for numClusters = k_rng
	disp(['K = ',num2str(numClusters)])
	load(fullfile(savedir,['kmeans','k_',num2str(numClusters),'.mat']));
	kClusterCentroids = GET_CENTROIDS(concTS,partition,numClusters);
	VarianceExplained(numClusters - 1) = VAREXPLAINED(concTS,partition,kClusterCentroids,numClusters);
end

save(fullfile(savedir,['VarianceExplained','.mat']),'VarianceExplained','k_rng');

% Fig S2a-b
f=figure;
plot(k_rng,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
%prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['ElbowPlot','.pdf']));

f = figure;
plot(k_rng(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
%prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['GainInVarianceExplained','.pdf']));
