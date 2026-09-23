%% compute subject-specific energy matrices (Figure 4a/b ii-iii) after choosing T with T_sweep_sps.m
clear all; close all;
addpath('C:\Users\marko\ketamine\energy_landscape-main\')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\control')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\plottingfxns')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\miscfxns')

sub_test = 'PTSD_ALL'; %'KET';'PTSD_CH_AC'; 'PTSD_COMP'; 'NIH' 'PTSD_ALL'; 'PTSD_ALL
title_1 = ' ses1';
title_2 = ' ses2';


[test, basedir, tokens, atlas, nparc] = param();

cd(basedir);
[ts_file_list, ordered_files_count] = get_files_list(basedir, tokens);
savedir = fullfile(basedir,test,'Energy');mkdir(savedir);		% set save directory

%% load BOLD data

c = 0;
T = 0.001; % set time scale parameters based on values from T_sweep_sps.m
load([savedir, '\concTS.mat']);
load([savedir, '\concTS_meta.mat']);

numClusters=6;

load(fullfile(savedir,['\Partition_bp','_k',num2str(numClusters),'.mat']),'partition','centroids','clusterNames')
[nsubjs, SES1_stop, SES2_start, SES2_stop,subjInd, partition, ses1_start_file, ses2_start_file]...   
    = get_ses1_ses_2(ts_file_list, partition, concTS, SCAN_LEN, scan_name, sub_test);

if contains(atlas,'SchafferTian-Yeo.xlsx') %Schaefer2018 don't have connectivity data for subcortical - drop it
    centroids = centroids(55:end, :);
    nparc = 400;
    sc = readmatrix('C:\Users\marko\ketamine\atlas\cs400_7NW.csv');
elseif contains(atlas,'C:\Users\marko\ketamine\energy_landscape-main\data\Lausanne_463_subnetworks.mat')
    indices_to_delete = [14, 463]; %Taly: from NAME_CLUSTERS_ANGLE.m file
    centroids(indices_to_delete, :) = [];
    nparc = 461;
else
    disp('Unsupported atlas')
end

overallCentroids=centroids;

load(fullfile(savedir,['\subjcentroids_split','_k',num2str(numClusters),'.mat']),'centroids')
if contains(atlas,'SchafferTian-Yeo.xlsx')
    centroids = centroids(:,55:end, :);
elseif contains(atlas,'C:\Users\marko\ketamine\energy_landscape-main\data\Lausanne_463_subnetworks.mat')
    centroids(:, indices_to_delete, :) = [];
end


%% energy matrices

Anorm = NORMALIZE(sc,c); 

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

E_full=NaN(nsubjs*2,numClusters^2);

e_i = 0;
for i=1:nsubjs*2 
    if i <= nsubjs
        shift = ses1_start_file-1;
    else
        shift = ses2_start_file-nsubjs-1;
    end
    ts_file_name = ts_file_list(i+shift).name;
    e_i = e_i+1;
    x0 = squeeze(centroids(i+shift,:,Xo_ind));
    xf = squeeze(centroids(i+shift,:,Xf_ind)); % now each column of x0 and xf represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full(e_i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition        
end


%% plot ses_1 vs ses_2 full 
Energy = E_full;
[~,pavg,~,t]=ttest(Energy(1:nsubjs,:),Energy(nsubjs+1:nsubjs*2,:));

fdravg = mafdr(pavg,'BHFDR',1);
fdravg = reshape(fdravg,[numClusters numClusters])';
pavg = reshape(pavg,[numClusters numClusters])';

grpAvgLSD = reshape(mean(Energy(1:nsubjs,:)),[numClusters numClusters])';
grpAvgPL = reshape(mean(Energy(nsubjs+1:nsubjs*2,:)),[numClusters numClusters])';


grpDiff = reshape(squeeze(t.tstat),[numClusters numClusters])';% .* -log(fdravg); (add sign() around tstat if want

maxVal = max(max([grpAvgLSD,grpAvgPL])); % sync color scales
minVal = min(min([grpAvgLSD,grpAvgPL]));

figure;
subplot(1,3,1);
imagesc(grpAvgLSD);
xticks(1:numClusters); yticks(1:numClusters); colormap
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title([strrep(sub_test, '_', ' '), title_1]);
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,2);
imagesc(grpAvgPL);
xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title([strrep(sub_test, '_', ' '), title_2]);
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,3);
LSDMinusPLTP = (grpDiff);%fdrpv1t; %((grpAvgPL - grpAvgLSD)); %switching order for manuscript figures
imagesc(LSDMinusPLTP); colormap('parula');%colormap('viridis');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Initial State'); xlabel('Final State');
sig_thresh = 0.05;
[y,x] = find(pavg < sig_thresh);
text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
[y,x] = find(fdravg < sig_thresh);
text(x-.15,y+.18,'**','Color','w','Fontsize', 36);
u_caxis_bound = max(max(LSDMinusPLTP));
l_caxis_bound = min(min(LSDMinusPLTP));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound]; 
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,2,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title([strrep(sub_test, '_', ' '), title_2, '<',  title_1]);
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');



energyTable = table('Size', [size(Energy, 1)*numel(clusterNames), 3], 'VariableTypes', {'string', 'double', 'string'}, 'VariableNames', {'sub', 'tp', 'state'});

% Add 36 empty columns with names [clusterNames{i}, clusterNames{j}]
for i = 1:numel(clusterNames)
    energyTable.(clusterNames{i}) = cell(height(energyTable), 1);
end

Energy_1 = Energy(1:nsubjs,:);
Energy_2 = Energy(nsubjs+1:nsubjs*2,:);
for sub = 1:nsubjs
    e_tp1 = reshape(Energy_1(sub,:),[numClusters numClusters])';
    e_tp2 = reshape(Energy_2(sub,:),[numClusters numClusters])';
    for i = 1:numel(clusterNames)
        row_tp1 = (sub-1)*numel(clusterNames)+i;
        row_tp2 = (nsubjs + sub - 1)*numel(clusterNames)+i;
        energyTable{row_tp1, 'sub'} = {ts_file_list(sub).name(5:9)};
        energyTable{row_tp2, 'sub'} = {ts_file_list(sub).name(5:9)};
        energyTable{row_tp1, 'tp'} = 1;
        energyTable{row_tp2, 'tp'} = 2;
        energyTable{row_tp1, 'state'} = {clusterNames(i)};
        energyTable{row_tp2, 'state'} = {clusterNames(i)};

        for j = 1:numel(clusterNames)
            energyTable{row_tp1, clusterNames{j}} = num2cell(e_tp1(i, j));
            energyTable{row_tp2, clusterNames{j}} = num2cell(e_tp2(i, j));
        end
    end
end

writetable(energyTable, ['C:\Users\marko\ketamine\document\','_', sub_test, '_full_energy.csv']);