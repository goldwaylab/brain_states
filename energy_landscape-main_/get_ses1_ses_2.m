
function [nsubjs, SES1_stop, SES2_start, SES2_stop,subjInd, partition, ses1_start_file, ses2_start_file]...   
    = get_ses1_ses_2(ts_file_list, partition, concTS, SCAN_LEN, scan_name, sub_test)
    
    tests = struct(...
        'KET', 'sub-K',...
        'PTSD_COMP', 'sub-F',...
        'PTSD_CH', 'sub-D',...
        'PTSD_A', 'sub-E'...
        );
    
    tp2 = 'ses-2';
    if contains(sub_test,'PTSD_CH_AC')
        prefix = [tests.PTSD_CH; tests.PTSD_A];
    elseif contains(sub_test,'PTSD_ALL')
        prefix = [tests.PTSD_CH; tests.PTSD_A; tests.PTSD_COMP];    
    elseif ~contains(sub_test,'all')
        prefix = [tests.(sub_test)];
    end

    subjInd = [];
    partition_new = [];
    current_start = 1;
    ts_file_list_new = [];
    nsubjs = 0;
    SES2_start = 0;
    SES2_stop = 0;
    nsubjs_save = 0;
    for ts_file_i = 1:size(ts_file_list, 1)
        ts_file_name = ts_file_list(ts_file_i).name;
        SCAN_LEN_i = SCAN_LEN(ts_file_i);
        current_end = current_start + SCAN_LEN_i - 1;
        TS_i = concTS(current_start:current_end);  

        update = false;
        if contains(sub_test,'all')
            update = true;
        else
            for prefix_i = 1:size(prefix,1)
                if contains(ts_file_name, prefix(prefix_i,:))                        
                    update = true;
                end
            end
        end
        if update  
            if nsubjs == 0
                ses1_start_file = ts_file_i;
            end
            nsubjs = nsubjs+1; 
            ts_file_list_new = [ts_file_list_new; ts_file_list(ts_file_i)];
            if and(contains(ts_file_name, tp2), SES2_start == 0) %assumes all ses-1 come first
                SES1_stop = size(partition_new,1);
                SES2_start = SES1_stop + 1;
                nsubjs = 1;
                ses2_start_file = ts_file_i;
            end
            partition_new = [partition_new; partition(current_start:current_start+SCAN_LEN_i-1)];           
            subjInd = [subjInd, repelem(nsubjs:nsubjs, SCAN_LEN_i)]; 
            disp([num2str(SCAN_LEN_i), ' ', ts_file_name])
            
        end
        current_start = current_start + SCAN_LEN_i;
        
    end
    if SES2_stop == 0
        SES2_stop = size(partition_new,1);
    else
        nsubjs = nsubjs_save-1;
    end
    partition = partition_new;
    ts_file_list = ts_file_list_new;
    disp(size(ts_file_list))
end
