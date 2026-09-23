%% generate subject specific centroids to make subj-specific comparisons and an TE comparison matrix

clear all; close all;
addpath('C:\Users\marko\ketamine\energy_landscape-main\misc_code')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\assesscluster')
addpath('C:\Users\marko\ketamine\energy_landscape-main\ejc_bs_code\plottingfxns')
addpath('C:\Users\marko\ketamine\energy_landscape-main')

sub_test = 'all';

[test, basedir, tokens, atlas, nparc] = param();
savedir = fullfile(basedir, test, 'Energy');
cd(basedir);

%% load BOLD data
[ts_file_list, ordered_files_count] = get_files_list(basedir, tokens);
load([savedir, '\concTS.mat']);
load([savedir, '\concTS_meta.mat']);

numClusters=6;
numNets=7;

clusterNames={};
clusterNamesUp={};
clusterNamesDown={};

%%

%load partition to calculate centroids on

load([savedir,'\Partition_bp','_k',num2str(numClusters),'.mat'],'partition','centroids', 'clusterNames'); %make sure file matches centroids you want to initialize from
sub_test = 'all';
[nsubjs, SES1_stop, SES2_start, SES2_stop,subjInd, partition, ses1_start_file, ses2_start_file]...   
    = get_ses1_ses_2(ts_file_list, partition, concTS, SCAN_LEN, scan_name, sub_test);
nscans = nsubjs*2;
SES1subjInd = subjInd;
SES1subjInd(:, SES1_stop+1:end) = 0;
SES2subjInd = subjInd;
SES2subjInd(:, 1:SES1_stop) = 0;
SES2subjInd(:, SES2_stop+1: end) = 0;

clusterNames_all = clusterNames;
clear centroids;
centroids=NaN(nscans,nparc,numClusters);

for i=1:nsubjs
    
    disp(['Generating SES1 centroids for subj ',num2str(i)]);
    centroids(i,:,:) = GET_CENTROIDS(concTS(SES1subjInd==i,:),partition(SES1subjInd==i),numClusters);
    [clusterNames{i},~,~,net7angle(i,:,:)] = NAME_CLUSTERS_ANGLE(squeeze(centroids(i,:,:)), atlas);
    [clusterNamesUp{i},clusterNamesDown{i},net7angle_Up(i,:,:),net7angle_Down(i,:,:)] = NAME_CLUSTERS_UP_DOWN(squeeze(centroids(i,:,:)), atlas);
   
    disp(['Generating SES2 centroids subj ',num2str(i)]);
    centroids(i+nsubjs,:,:) = GET_CENTROIDS(concTS(SES2subjInd==i,:),partition(SES2subjInd==i),numClusters);
    [clusterNames{i+nsubjs},~,~,net7angle(i+nsubjs,:,:)] = NAME_CLUSTERS_ANGLE(squeeze(centroids(i+nsubjs,:,:)), atlas);
    [clusterNamesUp{i+nsubjs},clusterNamesDown{i+nsubjs},net7angle_Up(i+nsubjs,:,:),net7angle_Down(i+nsubjs,:,:)] = NAME_CLUSTERS_UP_DOWN(squeeze(centroids(i+nsubjs,:,:)), atlas);
    
end


%% ttest on radial values


for i=1:numClusters
    [~,pUP(i,:),~,statUP{i}]=ttest(net7angle_Up(1:nsubjs,i,:),net7angle_Up(nsubjs+1:nsubjs*2,i,:));
    [~,pDOWN(i,:),~,statDOWN{i}]=ttest(net7angle_Down(1:nsubjs,i,:),net7angle_Down(nsubjs+1:nsubjs*2,i,:));
end

pfdrDN=reshape(mafdr(reshape(pDOWN,1,[]),'BHFDR',1),numClusters,numNets);
pfdrUP=reshape(mafdr(reshape(pUP,1,[]),'BHFDR',1),numClusters,numNets);


UPstat = NaN(numClusters,numNets);
DNstat = NaN(numClusters,numNets);
for i=1:numClusters
    UPstat(i,:) = squeeze(statUP{i}.tstat(:,:,1:7));
    DNstat(i,:) = squeeze(statDOWN{i}.tstat(:,:,1:7));
end

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
clusterNames = clusterNames_all;

% SI Figure 6b
figure;
subplot(1,2,1);
imagesc(UPstat); 
xticks(1:numNets); xticklabels(YeoNetNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Centroids'); xlabel('Networks');
sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
[y1,x1] = find(pUP < sig_thresh); 
text(x1-.12,y1+.12,'*','Color','w','Fontsize', 24);
[y,x] = find(pfdrUP < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 24);
u_caxis_bound = max(max(4));
l_caxis_bound = min(min(-4));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound u_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(u_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('UP Radial values from subj centroids ses1 > ses2');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');


subplot(1,2,2);
imagesc(DNstat); 
xticks(1:numNets); xticklabels(YeoNetNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Centroids'); xlabel('Networks');
sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
[y1,x1] = find(pDOWN < sig_thresh); 
text(x1-.12,y1+.12,'*','Color','w','Fontsize', 24);
[y,x] = find(pfdrDN < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
text(x-.12,y+.12,'**','Color','w','Fontsize', 24);
u_caxis_bound = max(max(4));
l_caxis_bound = min(min(-4));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound u_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(u_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('DOWN Radial values from subj centroids ses1 > ses2');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');



%% make a correlation plot comparing LSD average centroids with PL average centriods
load([savedir,'\Partition_bp','_k',num2str(numClusters),'.mat'],'clusterNames'); 
LSDcentroids = squeeze(mean(centroids(1:nsubjs,:,:))); %run these two pairs of centroids on systems_plot_sps.m to get SI Figure 6a/b
PLcentroids = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,:)));

err = NaN(numClusters,numClusters);
for i=1:numClusters
    for k=1:numClusters
        err(i,k)=immse(LSDcentroids(:,i),PLcentroids(:,k));
    end
end

[rcent,pcent]=corr(LSDcentroids,PLcentroids);
%SI Figure 6d
figure;
imagesc(rcent);
xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('SES1 Centroids'); xlabel('SES2 Centroids');
title('Correlation between Condition-Specific Centroids');
caxis_bound = 1;
h = colorbar; ylabel(h,'R'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

figure;
imagesc(err);
xticks(1:numClusters); yticks(1:numClusters); colormap('hot'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('SES1 Centroids'); xlabel('SES2 Centroids');
title('MSE between Condition-Specific Centroids');
u_caxis_bound = min(max(err));
l_caxis_bound = min(min(err));
h = colorbar; ylabel(h,'MSE'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound u_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(u_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');


%% plot radial plots for each subj/condition (can skip unless interested)

% YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};
% YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
% numNets = numel(YeoNetNames);
% clusterColors = GET_CLUSTER_COLORS(numClusters);
% clusterColors = hex2rgb(clusterColors);
% %YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;247 253 205; 218 166 86; 199 109 117] / 255;
% YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117] / 255;
% YeoColors = [YeoColors;YeoColors];
% netAngle = linspace(0,2*pi,numNets+1);
% thetaNames = YeoNetNames; thetaNames{8} = '';
% 
% 
% for i=1:nsubjs*2
% 
%     overallNames = clusterNames{i};
%     icentroids = squeeze(centroids(i,:,:));
%     inet7angle = squeeze(net7angle(i,:,:));
%     inet7angle_Up = squeeze(net7angle_Up(i,:,:));
%     inet7angle_Down = squeeze(net7angle_Down(i,:,:));
%     f(i)=figure;
%     for K = 1:numClusters
%         ax = subplot(1,numClusters,K,polaraxes); hold on
%         polarplot(netAngle,[inet7angle_Up(K,:) inet7angle_Up(K,1)],'k');
%         polarplot(netAngle,[inet7angle_Down(K,:) inet7angle_Down(K,1)],'r');
%         thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
%         rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
%         for L = 1:numNets
%             ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
%                 YeoColors(L,:), ax.ThetaTickLabel{L});
%         end
%         set(ax,'FontSize',10);
%         title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
%     end
%     f(i).PaperUnits = 'inches';
%     f(i).PaperSize = [8 1.5];
%     f(i).PaperPosition = [0 0 8 1.5];
%     
% end



%% save


save(fullfile(savedir,['subjcentroids_split','_k',num2str(numClusters),'.mat']),'centroids')
