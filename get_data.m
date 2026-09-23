function [concTS, scan_name, SCAN_LEN]  = get_data(ts_file_list)
    first = true;
    TS = 0;
    SCAN_LEN = 0;
    scan_name = cell(length(ts_file_list), 1);
    for subj_i = 1:length(ts_file_list)
        TS_i = readtable([ts_file_list(subj_i).folder,'\', ts_file_list(subj_i).name]); %  189 sample in AICHA parcellation
        TS_i = TS_i(2:end, :);
        SCAN_LEN_i = size(TS_i,1);
        if SCAN_LEN_i~=192
            disp(ts_file_list(subj_i).name)
        end
        scan_name{subj_i} = ts_file_list(subj_i).name;
        if(first == true)
            first = false;
            TS = TS_i;
            SCAN_LEN = SCAN_LEN_i;
        else
            TS = vertcat(TS, TS_i);
            SCAN_LEN = vertcat(SCAN_LEN, SCAN_LEN_i);
        end
    end
    concTS = table2array(TS);
end
