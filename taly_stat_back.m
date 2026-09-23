clear all

test = 'unified_28_1_23_Schaefer2018';
[weighted_results, maxNumClusters, minNumClusters, nsplits, nparc, basedir_, atlas] = setGlobalVar(test);
numClusters = 6;
for i = 1:size(weighted_results, 1)
    basedir = [basedir_, weighted_results{i}]; 

    variable_names = load([basedir, '\results\clusterAssignments\k',num2str(numClusters),'.mat']).clusterAssignments.(['k',num2str(numClusters)]).clusterNames  ;
    for ii = 1: size(variable_names,1)
        variable_names{ii} = [num2str(ii), '_', variable_names{ii} ];
    end

    subs_list = load([basedir, '\results\clusterAssignments\k',num2str(numClusters),'.mat']).clusterAssignments.(['k',num2str(numClusters)]).scan_name  ;
    half_sub = size(subs_list,1)/2;
    subs_list = subs_list(1:half_sub);
    subs = cell(1,size(subs_list,1));
    for s=1:size(subs_list,1)
        subs{1,s} = subs_list{s,1}(1:9);
    end
    subs = table(subs.', 'VariableNames', {'sub'});
    
        
    if contains(test, 'unified_28_1_23')
        list_of_files = load([basedir, '/results/analyses/transitionprobabilities/list_of_files.mat']).list_of_files;
        scanlab = {'Pl_ofir','Pl_Nitzan','Pl_Neomi', 'PTSD_AC_1', 'PTSD_AC_2', 'PTSD_Ch_1', 'PTSD_Ch2','Placebo_Noam', 'ket', 'PTSD_Naomi_1', 'PTSD_Naomi_Imp', 'PTSD_Naomi_NotImp'};
        all_tables_rr = [];
        all_tables_dwt = [];
        all_tables_fo = [];
        for group_i = 1:size(scanlab,2)
            sub_dataTable = cell2table(list_of_files(group_i,:)', 'VariableNames', {scanlab{group_i}});
            data_size = 80;

            data_rr = load([basedir, '\results\analyses\transitionprobabilities\', scanlab{group_i},'DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
            dataTable_rr = create_table(data_rr, data_size, scanlab{group_i}, variable_names);
            all_tables_rr = [all_tables_rr, dataTable_rr];

            data = load([basedir, '\results\analyses\transitionprobabilities\', scanlab{group_i},'DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
            dataTable = create_table(data, data_size, scanlab{group_i}, variable_names);
            all_tables_dwt = [all_tables_dwt, dataTable];

            data = load([basedir, '\results\analyses\transitionprobabilities\', scanlab{group_i},'FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
            dataTable = create_table(data, data_size, scanlab{group_i}, variable_names);
            all_tables_fo = [all_tables_fo, dataTable];

        end
        file_name = [basedir, '\results\analyses\rr_results.xlsx'];
        writetable(all_tables_rr, file_name, 'WriteVariableNames', true);
   
        file_name = [basedir, '\results\analyses\dwt_results.xlsx'];
        writetable(all_tables_dwt, file_name, 'WriteVariableNames', true);

        file_name = [basedir, '\results\analyses\fo_results.xlsx'];
        writetable(all_tables_fo, file_name, 'WriteVariableNames', true);

    elseif contains(test, 'NOAM_NAOMI_Dynamor_NIH')
        scannum = 0;
        load('C:\Users\marko\ketamine\NOAM_NAOMI_Dynamor_NIH\Raw\results\repkmeans\kmeans_meta.mat') %'ts_file_list', 'SCAN_LEN', 'scan_name', 'nsplits', 'ses_1', 'ses_2'
        improved = ["007", "008", "016", "018", "031", "034", "036", "047", "056",  "060"];
        all_subs = {};
        sub_i = 1;
        for file_i = 1:size(SCAN_LEN,1)
            if or(and(scannum == 0, contains(scan_name(file_i), 'rest1')),...
                  or(and(scannum == 1, contains(scan_name(file_i), 'ses-2')),...
                     or(and(scannum == 2, contains(scan_name(file_i), 'TP1')),...
                        or(and(scannum == 3, contains(scan_name(file_i), 'run-1')),...
                           and(scannum == 4, contains(scan_name(file_i), 'TP2'))))))
                scannum = scannum+1;
                sub_i = 1;                
            end
            if contains(scan_name(file_i), 'TP2')
                    if sub_i == 1
                        subjInd_i_ipmr = 1;
                        subjInd_i_not_impr = 1;
                        scannum_impr = scannum+1;
                        scannum_not_impr = scannum+2;
                    end
                sub_i_TP2 = scan_name(file_i);
                sub_i_TP2 = sub_i_TP2{1}(5:7);
                if sum(contains(improved, sub_i_TP2))
                    all_subs{scannum_impr, subjInd_i_ipmr} = scan_name{file_i,:};
                    subjInd_i_ipmr = subjInd_i_ipmr+1;
                else
                    all_subs{scannum_not_impr, subjInd_i_not_impr} = scan_name{file_i,:};
                    subjInd_i_not_impr = subjInd_i_not_impr+1;   
                end
            else
                all_subs{scannum+1, sub_i} = scan_name{file_i,:};
            end
            sub_i = sub_i + 1;
        end
        sub_dataTable_1 = cell2table(all_subs(1,:)', 'VariableNames', {'subs_rest_Noam'});
        sub_dataTable_2 = cell2table(all_subs(2,:)', 'VariableNames', {'subs_rest_Dynamore'});
        sub_dataTable_3 = cell2table(all_subs(3,:)', 'VariableNames', {'subs_ket'});
        sub_dataTable_4 = cell2table(all_subs(4,:)', 'VariableNames', {'subs_ptsd_naomi'});
        sub_dataTable_5 = cell2table(all_subs(5,:)', 'VariableNames', {'subs_ptsd_nih'});
        sub_dataTable_6 = cell2table(all_subs(6,:)', 'VariableNames', {'subs_ptsd_naomi_tp2_imroved'});
        sub_dataTable_7 = cell2table(all_subs(7,:)', 'VariableNames', {'subs_ptsd_naomi_tp2_not_imroved'});


        data_1 = load([basedir, '\results\analyses\transitionprobabilities\', 'rest_Noam','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_2 = load([basedir, '\results\analyses\transitionprobabilities\', 'rest_Dynamore','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_3 = load([basedir, '\results\analyses\transitionprobabilities\', 'ket','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_4 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_5 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_nih','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_6 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi_tp2_imroved','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_7 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi_tp2_not_imroved','DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        data_size = max(max(max(size(data_1,1), size(data_2,1)), max(size(data_3,1), size(data_4,1))), max(size(data_5,1), size(data_6,1)));
        dataTable_1 = create_table(data_1, data_size, 'rest_Noam', variable_names);
        dataTable_2 = create_table(data_2, data_size, 'rest_Dynamore', variable_names);
        dataTable_3 = create_table(data_3, data_size, 'ket', variable_names);
        dataTable_4 = create_table(data_4, data_size, 'ptsd_naomi', variable_names);
        dataTable_5 = create_table(data_5, data_size, 'ptsd_nih', variable_names);
        dataTable_6 = create_table(data_6, data_size, 'ptsd_naomi_tp2_imroved', variable_names);
        dataTable_7 = create_table(data_7, data_size, 'ptsd_naomi_tp2_not_imroved', variable_names);
        file_name = [basedir, '\results\analyses\rr_results.xlsx'];
        writetable([sub_dataTable_1, dataTable_1, sub_dataTable_2, dataTable_2, sub_dataTable_3, dataTable_3, ...
            sub_dataTable_4, dataTable_4, sub_dataTable_5, dataTable_5, sub_dataTable_6, dataTable_6,...
            sub_dataTable_7, dataTable_7], file_name, 'WriteVariableNames', true);

        data_1 = load([basedir, '\results\analyses\transitionprobabilities\', 'rest_Noam','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        data_2 = load([basedir, '\results\analyses\transitionprobabilities\', 'rest_Dynamore','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        data_3 = load([basedir, '\results\analyses\transitionprobabilities\', 'ket','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        data_4 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        data_5 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_nih','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        data_6 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi_tp2_imroved','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        data_7 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi_tp2_not_imroved','DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        dataTable_1 = create_table(data_1, data_size, 'rest_Noam', variable_names);
        dataTable_2 = create_table(data_2, data_size, 'rest_Dynamore', variable_names);
        dataTable_3 = create_table(data_3, data_size, 'ket', variable_names);
        dataTable_4 = create_table(data_4, data_size, 'ptsd_naomi', variable_names);
        dataTable_5 = create_table(data_5, data_size, 'ptsd_nih', variable_names);
        dataTable_6 = create_table(data_6, data_size, 'ptsd_naomi_tp2_imroved', variable_names);
        dataTable_7 = create_table(data_7, data_size, 'ptsd_naomi_tp2_not_imroved', variable_names);

        file_name = [basedir, '\results\analyses\dwt_results.xlsx'];
        writetable([sub_dataTable_1, dataTable_1, sub_dataTable_2, dataTable_2, sub_dataTable_3, dataTable_3, ...
            sub_dataTable_4, dataTable_4, sub_dataTable_5, dataTable_5, sub_dataTable_6, dataTable_6,...
            sub_dataTable_7, dataTable_7], file_name, 'WriteVariableNames', true);

        data_1 = load([basedir, '\results\analyses\transitionprobabilities\', 'rest_Noam','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        data_2 = load([basedir, '\results\analyses\transitionprobabilities\', 'rest_Dynamore','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        data_3 = load([basedir, '\results\analyses\transitionprobabilities\', 'ket','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        data_4 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        data_5 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_nih','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        data_6 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi_tp2_imroved','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        data_7 = load([basedir, '\results\analyses\transitionprobabilities\', 'ptsd_naomi_tp2_not_imroved','FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        dataTable_1 = create_table(data_1, data_size, 'rest_Noam', variable_names);
        dataTable_2 = create_table(data_2, data_size, 'rest_Dynamore', variable_names);
        dataTable_3 = create_table(data_3, data_size, 'ket', variable_names);
        dataTable_4 = create_table(data_4, data_size, 'ptsd_naomi', variable_names);
        dataTable_5 = create_table(data_5, data_size, 'ptsd_nih', variable_names);
        dataTable_6 = create_table(data_6, data_size, 'ptsd_naomi_tp2_imroved', variable_names);
        dataTable_7 = create_table(data_7, data_size, 'ptsd_naomi_tp2_not_imroved', variable_names);
        file_name = [basedir, '\results\analyses\fo_results.xlsx'];
        writetable([sub_dataTable_1, dataTable_1, sub_dataTable_2, dataTable_2, sub_dataTable_3, dataTable_3, ...
            sub_dataTable_4, dataTable_4, sub_dataTable_5, dataTable_5, sub_dataTable_6, dataTable_6,...
            sub_dataTable_7, dataTable_7], file_name, 'WriteVariableNames', true);    

    elseif dual
        data = load([basedir, '\results\analyses\transitionprobabilities\', str1,'DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        dataTable_ket = array2table(data, 'VariableNames', strcat(str1, variable_names));
        data = load([basedir, '\results\analyses\transitionprobabilities\', str2,'DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        dataTable_no_ket = array2table(data, 'VariableNames', strcat(str2, variable_names));
        file_name = [basedir, '\results\analyses\k',num2str(numClusters),'_rr_results.xlsx'];
        save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs, str1, str2, numClusters)

        data = load([basedir, '\results\analyses\transitionprobabilities\', str1,'DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        dataTable_ket = array2table(data, 'VariableNames', strcat(str1, variable_names));
        data = load([basedir, '\results\analyses\transitionprobabilities\', str2,'DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        dataTable_no_ket = array2table(data, 'VariableNames', strcat(str2, variable_names));
        file_name = [basedir, '\results\analyses\k',num2str(numClusters),'_dwt_results.xlsx'];
        save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs, str1, str2, numClusters)
        
        data = load([basedir, '\results\analyses\transitionprobabilities\', str1,'FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        dataTable_ket = array2table(data, 'VariableNames', strcat(str1, variable_names));
        data = load([basedir, '\results\analyses\transitionprobabilities\', str2,'FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        dataTable_no_ket = array2table(data, 'VariableNames', strcat(str2, variable_names));
        file_name = [basedir, '\results\analyses\k',num2str(numClusters),'_fo_results.xlsx'];
        save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs, str1, str2, numClusters)
    else
        data = load([basedir, '\results\analyses\transitionprobabilities\', str1,'DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        dataTable_ptsd = array2table(data, 'VariableNames', strcat(str1, variable_names));
        file_name = [basedir, '\results\analyses\rr_results.xlsx'];
        writetable([subs, dataTable_ptsd], file_name, 'WriteVariableNames', true);

        data = load([basedir, '\results\analyses\transitionprobabilities\', str1,'DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        dataTable_ptsd = array2table(data, 'VariableNames', strcat(str1, variable_names));
        file_name = [basedir, '\results\analyses\dwt_results.xlsx'];
        writetable([subs, dataTable_ptsd], file_name, 'WriteVariableNames', true);

        data = load([basedir, '\results\analyses\transitionprobabilities\', str1,'FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        dataTable_ptsd = array2table(data, 'VariableNames', strcat(str1, variable_names));
        file_name = [basedir, '\results\analyses\fo_results.xlsx'];
        writetable([subs, dataTable_ptsd], file_name, 'WriteVariableNames', true);    
    end
end    
    
function success = save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs, str1, str2, numClusters)
    emptyArray = cell(size(subs,1), numClusters+1);
    variable_names_ = ['Param';variable_names];
    stat = cell2table(emptyArray, 'VariableNames', variable_names_);
    stat{1,1}={'p-value '};
    stat{2,1}={[str1,'average']};
    stat{3,1}={[str2, 'average']};
    for i = 1:size(variable_names,1)
         p = [my_ttest(dataTable_ket{:, i}, dataTable_no_ket{:, i}, variable_names(i,1))];
         stat{1,i+1} = {p};
         stat{2,i+1} = {mean(dataTable_ket{:, i})};
         stat{3,i+1} = {mean(dataTable_no_ket{:, i})};
    end
    
    writetable([subs, dataTable_ket, dataTable_no_ket, stat], file_name, 'WriteVariableNames', true);
    success = true;
end
function [p] = my_ttest(group1, group2, title)

    [h, p, ci, stats] = ttest(group1, group2);
    % Display the results
    disp(title);
    fprintf('Paired t-test results:\n');
    fprintf('p-value: %f\n', p);
    fprintf('t-statistic: %f\n', stats.tstat);
    fprintf('Degrees of freedom: %d\n', stats.df);

end

function dataTable = create_table(data, targetRows, prefix, variable_names)
    numColumns = size(data, 2);
    expandedArray = nan(targetRows, numColumns);
    expandedArray(1:size(data, 1), :) = data;
    dataTable = array2table(expandedArray, 'VariableNames', strcat(prefix, variable_names));
end