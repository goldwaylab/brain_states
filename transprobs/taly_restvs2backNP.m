%{
Taly
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

scanlab = {'RestComb','TwoBackComb'}; 

savedir = [masterdir,'/analyses/transitionprobabilities'];

rng('shuffle');
%}
%Taly

clear all
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\nullfxns')
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\plottingfxns')
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\miscfxns')
basedir_ = 'C:\Users\marko\ketamine\Results_No_Physio\';
weighted_results = {'Raw'};  %; 'NMDA'; 'opioid_kappa_1'; 'opioid_mu_1'; 'dopamine_D2'};
scanlab = {'No_Ket', 'Ket'}; 
numClusters = 6;
name_root = '';
for i = 1:size(weighted_results, 1)
    basedir = [basedir_, weighted_results{i}];
    masterdir = [basedir, '\results'];
    savedir =  [masterdir, '\transitions'];

    %% load transition probabilities
    disp('loading transition probabilities')
    
    NoKet = load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilitiesNoPersist_k',num2str(numClusters),name_root,'.mat']);
    NoKetTransitionProbabilityMats = NoKet.transitionProbabilityMats;
    Ket = load([masterdir,'/analyses/transitionprobabilities/',scanlab{2},'TransitionProbabilitiesNoPersist_k',num2str(numClusters),name_root,'.mat']);
    KetTransitionProbabilityMats = Ket.transitionProbabilityMats;
    
    load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
    clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

    [meta_state_names, meta_state_pairs] = get_meta_states(clusterNames);
    numMetaStates = size(meta_state_names, 2);
    NoKetTransitionProbabilityMats_meta = zeros([size(NoKetTransitionProbabilityMats,1), numMetaStates, numMetaStates]);
    KetTransitionProbabilityMats_meta = zeros([size(KetTransitionProbabilityMats,1), numMetaStates, numMetaStates]);
    for p_i = 1:size(NoKetTransitionProbabilityMats,1)
        NoKetTransitionProbabilityMats_meta(p_i, :, :) = ...
            unite_states(squeeze(NoKetTransitionProbabilityMats(p_i,:,:)), meta_state_pairs, meta_state_names);
        KetTransitionProbabilityMats_meta(p_i, :, :) = ...
            unite_states(squeeze(KetTransitionProbabilityMats(p_i,:,:)), meta_state_pairs, meta_state_names);
    end    
    
    %% permute rest and 2-back within subjects
    
    nperms = 100000;
    pvals_twotail = PERM_TEST(KetTransitionProbabilityMats_meta,NoKetTransitionProbabilityMats_meta,nperms);
    p_values = sort(pvals_twotail(:));
    non_zero_index = find(p_values ~= 0, 1);
    p_values = p_values(non_zero_index:end);
    q_values = mafdr(p_values, 'BHFDR', true);
    pvals_twotail_fdr = zeros(size(pvals_twotail));
    count = 1;
    for i_meta = 1:size(pvals_twotail_fdr,1)
        for j_meta = 1:size(pvals_twotail_fdr,2)
            if i_meta ~= j_meta
                pvals_twotail_fdr(i_meta, j_meta) = q_values(count);
                count = count+1;
            end
        end
    end

    success = t_test_(KetTransitionProbabilityMats_meta,NoKetTransitionProbabilityMats_meta, meta_state_names, savedir);
    
    %% Plot 
    disp('plot')
    
    cd(savedir);
    
    grpAvgNoKet = squeeze(nanmean(NoKetTransitionProbabilityMats_meta,1));
    grpAvgNoKet = grpAvgNoKet .* ~eye(numMetaStates);
    grpAvgKet = squeeze(nanmean(KetTransitionProbabilityMats_meta,1)); %restGrpAvg(logical(eye(numMetaStates))) = 0;
    grpAvgKet = grpAvgKet .* ~eye(numMetaStates);
    maxval = max(max([grpAvgNoKet,grpAvgKet]))
    
    f=figure;
    subplot(1,3,1);
    imagesc(grpAvgNoKet); 
    caxis([0 maxval]); colorbar
    xticks(1:numMetaStates); yticks(1:numMetaStates); %colormap('plasma');
    xticklabels(meta_state_names); xtickangle(90); yticklabels(meta_state_names); axis square;
    COLOR_TICK_LABELS(true,true,numMetaStates);
    ylabel('Current State'); xlabel('Next State');
    title(scanlab{1});
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    
    subplot(1,3,2); 
    imagesc(grpAvgKet); 
    caxis([0 maxval]); colorbar
    xticks(1:numMetaStates); yticks(1:numMetaStates); %colormap('plasma');
    xticklabels(meta_state_names); xtickangle(90); yticklabels(meta_state_names); axis square;
    COLOR_TICK_LABELS(true,true,numMetaStates);
    ylabel('Current State'); xlabel('Next State');
    title(scanlab{2});
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    
    subplot(1,3,3);
    TwoBackMinusRestTPMat = (grpAvgKet-grpAvgNoKet);
    imagesc(TwoBackMinusRestTPMat.*~eye(numMetaStates)); %colormap('plasma');
    xticks(1:numMetaStates); yticks(1:numMetaStates)
    xticklabels(meta_state_names); xtickangle(90); yticklabels(meta_state_names); axis square;
    COLOR_TICK_LABELS(true,true,numMetaStates);
    ylabel('Current State'); xlabel('Next State');
    sig_thresh = 0.05 / (numMetaStates^2-numMetaStates);      % bonferroni correction, for two-tailed p-values so only
    [y,x] = find(pvals_twotail.*~eye(numMetaStates) < sig_thresh);
    text(x-.12,y+.12,'*','Color','w');
    %h=colorbar; ylabel(h,'-log_{10}(p)'); caxis([0 log10(nperms)]); h.Ticks = [0 ut lt log10(nperms)]; h.TickLabels = [0 round(ut,2,'significant') round(lt,2,'significant') log10(nperms)];
    caxis_bound = max(max(abs(TwoBackMinusRestTPMat.*~eye(numMetaStates))));
    h = colorbar; ylabel(h,[scanlab{2}, '-', scanlab{1}]); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
    COLOR_TICK_LABELS(true,true,numMetaStates);
    title([scanlab{2}, '>', scanlab{1}]);
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    
    f.PaperUnits = 'inches';
    f.PaperSize = [8 4];
    f.PaperPosition = [0 0 8 4];
    save('Fig4a-c,e-f__SourceData.mat','grpAvgNoKet','grpAvgKet','TwoBackMinusRestTPMat','pvals_twotail')
    saveas(f,['RestvsTwoBackNoPersist_nonpar_met',num2str(numMetaStates),'.pdf']);
    
    %% Persistence probabilities: rest vs n-back
    %{
    % retrieve diagonal again
    grpAvgRest = squeeze(nanmean(restTransitionProbabilityMats,1));    
    grpAvgTwoBack = squeeze(nanmean(TwoBackTransitionProbabilityMats,1));
    TwoBackMinusRestTPMat = (grpAvgTwoBack-grpAvgRest);
    
    [y,x] = find(diag(pvals_twotail)' < sig_thresh);
    f=figure;
    imagesc(diag(TwoBackMinusRestTPMat)');
    caxis_bound = max(max(abs(diag(TwoBackMinusRestTPMat))));
    xticks([]); yticks([]);
    text(x-.12,y+.12,'*','Color','w');
    caxis([-caxis_bound caxis_bound]); colormap('plasma'); %colorbar
    %h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
    disp(['Persistence probability 2-back - rest colorbar bounds: +/-',num2str(caxis_bound)])
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    f.PaperUnits = 'inches';
    f.PaperSize = [1 .2];
    f.PaperPosition = [0 0 1 .2];
    saveas(f,['PersistRestvsTwoBack_nonpar_k',num2str(numMetaStates),'.pdf']);
    
    %}
end

function meta_states = unite_states(states, meta_state_pairs, meta_state_names)
%unit states probabilities to move to another state to meta state
%probabilities
    meta_states = zeros(3, 3);
    for meta_1 = 1:3
        for meta_2 = 1:3
            meta_states(meta_1,meta_2) = ...
                (states(meta_state_pairs(meta_1, 1),meta_state_pairs(meta_2, 1)) + ...
                states(meta_state_pairs(meta_1, 1),meta_state_pairs(meta_2, 2)) + ...
                states(meta_state_pairs(meta_1, 2),meta_state_pairs(meta_2, 1)) + ...
                states(meta_state_pairs(meta_1, 2),meta_state_pairs(meta_2, 2)))/2;
        end
    end
    meta_states = meta_states.*~eye(size(meta_state_names,2));
    meta_states = meta_states./sum(meta_states,2);
end

function [meta_state_names, meta_state_pairs] = get_meta_states(clusterNames)
    clusterNames_ = clusterNames;
    meta_state_pairs = zeros(3, 2); 
    meta_state_names = cell(1, 3);
    pair = 0;
    for clusterName_i = 1:size(clusterNames_)
        clusterName = clusterNames_(clusterName_i);
        clusterName = clusterName{1};
        if clusterName == 'USED'
            continue
        end
        pair = pair + 1;
        meta_state_pairs(pair, 1) = clusterName_i;
        meta_state = clusterName(1:3);
        meta_state_names{pair} = meta_state;
        if meta_state == 'DMN'
            second_name = 'VAT';
        else
            second_name = 'xxx';
        end 
        for second_name_i = clusterName_i+1: size(clusterNames_)
            meta_2 = clusterNames(second_name_i);
            meta_2 = meta_2{1};
            meta_2 = meta_2(1:3);
            if or(meta_2 == meta_state,  meta_2 == second_name)
                meta_state_pairs(pair, 2) = second_name_i;
                clusterNames_{second_name_i} = 'USED';
                break
            end
        end
    end
end

function success = t_test_(KetTransitionProbabilityMats_meta,NoKetTransitionProbabilityMats_meta, meta_state_names, savedir)
    success = true;
    dataTable = table();
    for i = 1:size(KetTransitionProbabilityMats_meta, 2)
        for j = 1:size(KetTransitionProbabilityMats_meta, 2)
            if i == j
                continue
            end
            ket = KetTransitionProbabilityMats_meta(:, i, j);
            dataTable = [dataTable, array2table(ket, 'VariableNames', ...
                {['ket_', meta_state_names{i}, '_', meta_state_names{j}]})];
            not_ket = NoKetTransitionProbabilityMats_meta(:, i, j);
            dataTable = [dataTable, array2table(not_ket, 'VariableNames', ...
                 {['no_ket_', meta_state_names{i}, '_', meta_state_names{j}]})];
            [h, p, ci, stats] = ttest(ket, not_ket);
            text = ['from ', meta_state_names{i}, ' to ', meta_state_names{j},...
                ': no ket av=', num2str(mean(not_ket)), ', ket av=',num2str(mean(ket)), ', p:', num2str(p)];
            disp(text);
        end
    end
    file_name = [savedir, '\tr.xlsx'];
    writetable(dataTable, file_name, 'WriteVariableNames', true);
end

