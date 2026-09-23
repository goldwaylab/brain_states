clear all
addpath('C:\Users\marko\ketamine\brain_states\brain_states-master\code\miscfxns\');
numClusters = 6;
YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};

[test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path] = setGlobalVar();
load([basedir, test,'\repkmeans\kmeans_meta.mat'])
load([basedir, test,'\repkmeans\concTS.mat'])
load(fullfile(basedir, test,'clusterAssignments',['k',num2str(numClusters),'.mat']));
savedir = [basedir, test, '\analyses'];
RESULT = [savedir, '\all_FO_long_clinical_updated_nw.xlsx'];
summary_file = [savedir, '\nw_act.xlsx'];
summary = readtable(summary_file);
emptyCol = NaN(size(summary, 1), 1);
clusterNames = clusterAssignments.k6.clusterNames;
for j = 1:size(clusterNames,1)
    for i = 1:size(YeoNetNames,2)
        summary.([clusterNames{j}, '_', YeoNetNames{i}, '+']) = emptyCol;
        summary.([clusterNames{j}, '_', YeoNetNames{i}, '-']) = emptyCol;
    end
end

partition = clusterAssignments.k6.partition;
sub_test = 'PTSD_ALL';
[nsubjs, SES1_stop, SES2_start, SES2_stop,subjInd, partition, ses1_start_file, ses2_start_file, concTS, ts_file_list]...   
    = get_ses1_ses_2(ts_file_list, partition, concTS, SCAN_LEN, scan_name, sub_test);
start_subj = 1;
for scans = 1: nsubjs
    end_subj = start_subj + SCAN_LEN(scans)-1;
    sub_centroids_TP1 = zeros(nparc,numClusters);
    sub_centroids_TP2 = zeros(nparc,numClusters);
    concTS_sub_TP1 = concTS(start_subj:end_subj, :);
    concTS_sub_TP2 = concTS(SES1_stop+start_subj:SES1_stop+end_subj, :);
    partition_sub_TP1 = partition(start_subj:end_subj);
    partition_sub_TP2 = partition(SES1_stop+start_subj:SES1_stop+end_subj);
    for K = 1:numClusters
        sub_centroids_TP1(:,K) = mean(concTS_sub_TP1(partition_sub_TP1 == K,:),1)';
        sub_centroids_TP2(:,K) = mean(concTS_sub_TP2(partition_sub_TP2 == K,:),1)';
    end

    [~,~,net7angle_Up_TP1,net7angle_Down_TP1] = NAME_CLUSTERS_UP_DOWN(sub_centroids_TP1, atlas);
    [~,~,net7angle_Up_TP2,net7angle_Down_TP2] = NAME_CLUSTERS_UP_DOWN(sub_centroids_TP2, atlas);

    for state = 1:numClusters
        for nw_i = 1:size(YeoNetNames,2)
            summary.([clusterNames{state}, '_', YeoNetNames{nw_i}, '+'])(scans) = net7angle_Up_TP2(state, nw_i)-net7angle_Up_TP1(state, nw_i);
            summary.([clusterNames{state}, '_', YeoNetNames{nw_i}, '-'])(scans) = net7angle_Down_TP2(state, nw_i)-net7angle_Down_TP1(state, nw_i);
        end
    end
    start_subj = end_subj+1;
end
writetable(summary, RESULT);

