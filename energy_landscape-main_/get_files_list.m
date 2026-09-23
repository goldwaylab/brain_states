function [ts_file_list, ordered_files_count] = get_files_list(basedir, tokens)
        ts_file_list_1_ = [];
        ts_file_list_2_ = [];
        ordered_files_count_1 = [];
        ordered_files_count_2 = [];
        for i = 1:2:size(tokens,1)
            if i+1 > size(tokens,1)
                break
            end
            ts_file_list_1 = dir(fullfile(basedir, [tokens{i}, '*.csv']));
            ts_file_list_2 = dir(fullfile(basedir, [tokens{i+1}, '*.csv']));
            [ts_file_list_comb, ordered_files_count_comb] = get2RunsData(ts_file_list_1, ts_file_list_2);
            ts_file_list_1_ = [ts_file_list_1_; ts_file_list_comb(1:ordered_files_count_comb(1))];
            ts_file_list_2_ = [ts_file_list_2_; ts_file_list_comb(ordered_files_count_comb(1)+1:end)];
            ordered_files_count_1 = [ordered_files_count_1, ordered_files_count_comb(1)];
            ordered_files_count_2 = [ordered_files_count_2, ordered_files_count_comb(1)];
            
        end
        ordered_files_count = [ordered_files_count_1, ordered_files_count_2];
        ts_file_list=[ts_file_list_1_; ts_file_list_2_];
end

function [ts_file_list, ordered_files_count] = get2RunsData(ts_file_list_1, ts_file_list_2)
    % Extract subjects from file structures
    subjects_1 = cell(size(ts_file_list_1));
    for subj_i = 1:numel(ts_file_list_1)
        index = strfind(ts_file_list_1(subj_i).name, '_ses-')-1;
        subjects_1{subj_i} = ts_file_list_1(subj_i).name(1:index);
    end
    
    subjects_2 = cell(size(ts_file_list_2));
    for subj_i = 1:numel(ts_file_list_2)
        index = strfind(ts_file_list_1(subj_i).name, '_ses-')-1;
        subjects_2{subj_i} = ts_file_list_2(subj_i).name(1:index);
    end
    
    % Find common subjects
    common_subjects = intersect(subjects_1, subjects_2);
    
    % Filter ts_file_list_1 and ts_file_list_2
    filtered_list_1 = ts_file_list_1(ismember(subjects_1, common_subjects));
    filtered_list_2 = ts_file_list_2(ismember(subjects_2, common_subjects));

    for subj_i=1:size(filtered_list_1,1)
        if ~strcmp(filtered_list_1(subj_i).name(1:index),filtered_list_2(subj_i).name(1:index))
            disp([filtered_list_1(subj_i).name(1:index)], filtered_list_2(subj_i).name(1:index))
        end
    end
    
    ts_file_list = [filtered_list_1; filtered_list_2];
    ses_1 = size(filtered_list_1,1);
    ses_2 = size(filtered_list_2,1);
    ordered_files_count = [ses_1, ses_2];
end

