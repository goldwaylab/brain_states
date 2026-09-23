%% compute subject-specific energy matrices (Figure 4a/b ii-iii) after choosing T with T_sweep_sps.m
clear all; close all;
addpath('C:\Users\marko\ketamine\energy_landscape-main')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\control')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\plottingfxns')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\miscfxns')

[test, basedir, tokens, atlas, nparc] = param();
sub_test = 'all';
cd(basedir);
savedir = fullfile(basedir,test,'Energy');mkdir(savedir);		% set save directory

%% load BOLD data

c = 0;
T = 0.001; % set time scale parameters based on values from T_sweep_sps.m

numClusters=6;
load([savedir, '\concTS.mat']);
load([savedir, '\concTS_meta.mat']);

load(fullfile(savedir,['\Partition_bp','_k',num2str(numClusters),'.mat']),'partition','centroids','clusterNames')
overallCentroids=centroids(55:end, :);
nparc = 400;
%sc = readmatrix('C:\Users\marko\ketamine\atlas\cs400_7NW.csv');
sc = readmatrix('C:\Users\marko\ketamine\atlas\Schaefer2018_400_cs.csv');

load(fullfile(savedir,['\subjcentroids_split','_k',num2str(numClusters),'.mat']),'centroids')
centroids = centroids(:,55:end, :);
[ts_file_list, ordered_files_count] = get_files_list(basedir, test, tokens);
[nsubjs, SES1_stop, SES2_start, SES2_stop,subjInd, partition, ses1_start_file, ses2_start_file]...   
    = get_ses1_ses_2(ts_file_list, partition, concTS, SCAN_LEN, scan_name, sub_test);


%% energy matrices

Anorm = NORMALIZE(sc,c); 

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

E_full=NaN(nsubjs*2,numClusters^2);
E_weighted=NaN(nsubjs*2,numClusters^2);

for i=1:nsubjs*2 
    x0 = squeeze(centroids(i,:,Xo_ind));
    xf = squeeze(centroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
    
    % compute weighted control energy:
    
end

%% plot PL weighted vs PL full



%% plot LSD vs PL 

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
xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,2);
imagesc(grpAvgPL)
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('PL');
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
title('PL > LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');



