clear all
basedir_ = 'C:\Users\marko\ketamine\Results_No_Physio\';
weighted_results = {'Raw'; 'NMDA'; 'opioid_kappa_1'; 'opioid_mu_1'; 'dopamine_D2'};
for i = 3:3%size(weighted_results, 1)
    basedir = [basedir_, weighted_results{i}];

    variable_names = load([basedir, '\results\clusterAssignments\k5.mat']).clusterAssignments.k5.clusterNames  ;
    for ii = 1: size(variable_names,1)
        variable_names{ii} = [num2str(ii), '_', variable_names{ii} ];
    end
    subs = {load([basedir, '\results\inpu_files_list_ses_1.mat']).ts_file_list_1.name};
    for s=1:size(subs,2)
        subs{1,s} = subs{1,s}(1:8);
    end
    subs = table(subs.', 'VariableNames', {'sub'});
    
    
    data = load([basedir, '\results\analyses\transitionprobabilities\\KetDwellTime_k5.mat']).RunRate  ;
    dataTable_ket = array2table(data, 'VariableNames', strcat('ket_', variable_names));
    data = load([basedir, '\results\analyses\transitionprobabilities\No_KetDwellTime_k5.mat']).RunRate  ;
    dataTable_no_ket = array2table(data, 'VariableNames', strcat('no_ket_', variable_names));
    file_name = [basedir, '\results\analyses\rr_results_k5.xlsx'];
    save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs)
    
    data = load([basedir, '\results\analyses\transitionprobabilities\KetDwellTime_k5.mat']).DwellTimeMean  ;
    dataTable_ket = array2table(data, 'VariableNames', strcat('ket_', variable_names));
    data = load([basedir, '\results\analyses\transitionprobabilities\No_KetDwellTime_k5.mat']).DwellTimeMean  ;
    dataTable_no_ket = array2table(data, 'VariableNames', strcat('no_ket_', variable_names));
    file_name = [basedir, '\results\analyses\dwt_results_k5.xlsx'];
    save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs)
    
    data = load([basedir, '\results\analyses\transitionprobabilities\KetFractionalOccupancy_k5.mat']).FractionalOccupancy  ;
    dataTable_ket = array2table(data, 'VariableNames', strcat('ket_', variable_names));
    data = load([basedir, '\results\analyses\transitionprobabilities\No_KetFractionalOccupancy_k5.mat']).FractionalOccupancy  ;
    dataTable_no_ket = array2table(data, 'VariableNames', strcat('no_ket_', variable_names));
    file_name = [basedir, '\results\analyses\fo_results_k5.xlsx'];
    save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs)
    
end    
    
function success = save_stat(dataTable_ket, dataTable_no_ket, file_name, variable_names, subs)
    emptyArray = cell(size(subs,1), 6);
    variable_names_ = ['Param';variable_names];
    stat = cell2table(emptyArray, 'VariableNames', variable_names_);
    stat{1,1}={'p-value '};
    stat{2,1}={'ket average'};
    stat{3,1}={'no ket average'};
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