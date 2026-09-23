clear all

[test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path] = setGlobalVar();

numClusters = 6;
masterdir = [basedir, test];
scanlab = load([masterdir, '/analyses/transitionprobabilities/scanlab.mat']).scanlab;
variable_names = load([masterdir, '/clusterAssignments/k',num2str(numClusters),'.mat']).clusterAssignments.(['k',num2str(numClusters)]).clusterNames  ;
for variable_name_i = 1:size(variable_names,1)
    variable_names{variable_name_i} = [variable_names{variable_name_i},'_',num2str(variable_name_i),'_'] ;
end
list_of_files = load([masterdir, '/analyses/transitionprobabilities/list_of_files.mat']).list_of_files;        

create_tables(masterdir, numClusters, scanlab, list_of_files, variable_names, true);



function success = create_tables(basedir, numClusters, scanlab, list_of_files, variable_names, enriched)
    for ii = 1: size(variable_names,1)
        variable_names{ii} = [num2str(ii), '_', variable_names{ii} ];
    end

    all_tables_rr = [];
    all_tables_dwt = [];
    all_tables_fo = [];

    for group_i = 1:size(scanlab,2)
        sub_dataTable = cell2table(list_of_files(group_i,:)', 'VariableNames', {scanlab{group_i}});
        data_size = size(sub_dataTable,1);

        data = load([basedir, '\analyses\transitionprobabilities\', scanlab{group_i},'DwellTime_k',num2str(numClusters),'.mat']).RunRate  ;
        dataTable_rr = create_table(data, data_size, scanlab{group_i}, variable_names);

        data = load([basedir, '\analyses\transitionprobabilities\', scanlab{group_i},'DwellTime_k',num2str(numClusters),'.mat']).DwellTimeMean  ;
        dataTable_dwt = create_table(data, data_size, scanlab{group_i}, variable_names);

        data = load([basedir, '\analyses\transitionprobabilities\', scanlab{group_i},'FractionalOccupancy_k',num2str(numClusters),'.mat']).FractionalOccupancy  ;
        dataTable_fo = create_table(data, data_size, scanlab{group_i}, variable_names);

        if enriched
            all_tables_rr = [all_tables_rr, sub_dataTable, dataTable_rr];
            all_tables_dwt = [all_tables_dwt, sub_dataTable, dataTable_dwt];
            all_tables_fo = [all_tables_fo, sub_dataTable, dataTable_fo];
        else
            all_tables_rr = [all_tables_rr, dataTable_rr];
            all_tables_dwt = [all_tables_dwt, dataTable_dwt];
            all_tables_fo = [all_tables_fo, dataTable_fo];
        end

    end
    file_name = [basedir, '\analyses\rr_results.xlsx'];
    writetable(all_tables_rr, file_name, 'WriteVariableNames', true);

    file_name = [basedir, '\analyses\dwt_results.xlsx'];
    writetable(all_tables_dwt, file_name, 'WriteVariableNames', true);

    file_name = [basedir, '\analyses\fo_results.xlsx'];
    writetable(all_tables_fo, file_name, 'WriteVariableNames', true);

    success = true;
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