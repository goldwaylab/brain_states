function [list_of_files, all_TRs, numTRs_all, subjInd, scanInd, all_nobs, improved] = get_variables(test, scanlab, TR, masterdir, TR_path)
    max_subj = 100;
    kmeans_meta = load([masterdir,'\repkmeans\kmeans_meta.mat']);
    SCAN_LEN = kmeans_meta.SCAN_LEN;
    scan_name = kmeans_meta.scan_name;
    ordered_files_count = kmeans_meta.ordered_files_count;

    improved=[];
    all_nobs = [];
    subjInd = [];
    scanInd = [];
    list_of_files = cell(size(scanlab,2), max_subj);
    all_TRs = cell(size(scanlab,2), max_subj);
    numTRs_all = cell(size(scanlab,2), max_subj);
    scannum = 0;
    subjInd_i = 1;
    if TR == -1
        read_TR = true;
    else
        read_TR = false;
    end
    for scan_name_i = 1:size(scan_name, 1)    
        name = scan_name{scan_name_i};
        scanInd = [scanInd, ones(1, SCAN_LEN(scan_name_i)) * scannum];
        subjInd = [subjInd, ones(1, SCAN_LEN(scan_name_i)) * subjInd_i];
        list_of_files(scannum+1,subjInd_i) = scan_name(scan_name_i);
        if and(read_TR, name(5)~='K')
            jsonString = fileread([TR_path, name(1:end-3), 'json']);
            jsonData = jsondecode(jsonString);
            TR = jsonData.RepetitionTime;
        elseif name(5) =='K'
            TR = 2.5;
        end
        all_TRs (scannum+1,subjInd_i)= {TR};
        numTRs_all (scannum+1,subjInd_i)= {SCAN_LEN(scan_name_i)};
        
        if (scannum+1) <= size(ordered_files_count, 2)
            if subjInd_i == ordered_files_count(scannum+1)
                all_nobs = [all_nobs, subjInd_i];
                subjInd_i = 1;
                scannum = scannum+1;
            else
                subjInd_i = subjInd_i+1;
            end
        else
            subjInd_i = subjInd_i+1;
        end
    end
    all_nobs = [all_nobs, subjInd_i-1];
    
end

