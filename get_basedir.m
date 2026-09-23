function basedir_ = get_basedir(test)
    test_data_row = {
        'PTSD_KET_Schaefer2018_tables_subs', 'E:\ptsd_ketamine\';
        'PAIN_KET_Schaefer2018_tables_subs', 'E:\ketamine\pain\';
        'JOY_Schaefer2018_tables_subs', 'E:\joy\'
        'JOY_DOUBLE_Schaefer2018_tables_subs', 'E:\joy\'
        };
    count = 0;
    for test_i = 1:size(test_data_row, 1)
        if contains(test, test_data_row{test_i, 1})
            count = count + 1;
            basedir_ = test_data_row{test_i, 2};
            disp(basedir_)
        end
    end
    if count > 1
        disp('More than 1 entry')
    end
end