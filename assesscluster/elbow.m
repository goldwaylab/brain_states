%  make elbow plot (Fig S2a-b) of within cluster / total variance explained by clusters at different k values
% assumes correlation distance
%{
Taly
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
%}
%Taly

clear all; close all;clc
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\plottingfxns');
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code')

[test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path] = setGlobalVar();
masterdir = [basedir, test];
distanceMethod = 'correlation';
nreps = 1;

%load(fullfile(basedir,['data/Demographics','.mat']));
repkmeans = '\repkmeans';
load([masterdir, repkmeans, '\kmeans_meta.mat']) %'ts_file_list', 'SCAN_LEN', 'scan_name', 'nsplits'
    
%addpaths;
%load(fullfile(basedir,['data/Demographics','.mat']));
%load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
%concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
%concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
first = true;
TS = 0;
for file_i = 1:length(ts_file_list)
    TS_i = readtable([ts_file_list(file_i).folder,'\', ts_file_list(file_i).name]); %  189 sample in AICHA parcellation
    TS_i = TS_i(2:end, :);
    disp(ts_file_list(file_i).name)


    if(first == true)
        first = false;
        TS = TS_i;
    else
        TS = vertcat(TS, TS_i);
    end
end
concTS = table2array(TS);
%Taly   

N = size(concTS,1); % number of observations
savedir = fullfile(masterdir,'analyses','choosing_k'); mkdir(savedir);

k_rng = minNumClusters:maxNumClusters;
VarianceExplained = zeros(length(k_rng),1);

for numClusters = k_rng
    disp(['K = ',num2str(numClusters)])
    load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),'.mat']));
    partition = clusterAssignments.(['k',num2str(numClusters)]).partition;
    kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
    VarianceExplained(numClusters - 1) = VAREXPLAINED(concTS,partition,kClusterCentroids,numClusters);
end

save(fullfile(savedir,['VarianceExplained.mat']),'VarianceExplained','k_rng');

% Fig S2a-b
f=figure;
plot(k_rng,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['ElbowPlot.pdf']));

f = figure;
plot(k_rng(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['GainInVarianceExplained.pdf']));
