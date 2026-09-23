function [test, maxNumClusters, minNumClusters, nsplits, nparc, basedir, atlas, tokens, scanlab, TR, TR_path, PAIRED_RUNS] = setGlobalVar()
    %{ 
        test tokens:
        First part of the string will define raw data path and tokens to
        distinguish between runs
        Schaefer2018/AICHA/Lausanne - parcelation atlas to use    
    %}
    test = 'JOY_DOUBLE_Schaefer2018_tables_subs';

    minNumClusters = 2;
    maxNumClusters = 10;
    nsplits = 100;
    TR_path = '';
    PAIRED_RUNS = true;

    basedir = get_basedir(test); 

    if contains(test, 'Schaefer')
        disp('Schaefer')
        nparc = 454;
        atlas = 'E:\ptsd_ketamine\atlas\SchafferTian-Yeo.xlsx';  
        basedir = [basedir, 'Results_Schaefer2018_7Networks\'];
    elseif contains(test, 'Lausanne')
        nparc = 463;
        atlas = 'C:\Users\marko\ketamine\energy_landscape-main\energy_landscape-main\data\Lausanne_463_subnetworks.mat';
    else
        nparc = 384; %Number of ROIs - AICHA
        atlas = 'C:\Users\marko\ketamine\atlas\AICHA-Yeo.xlsx';
    end


    TR = -1;    
    if contains(test, 'PTSD_KET')
        tokens = {'sub-F*ses-1'; 'sub-F*ses-2';
                  'sub-E*ses-1'; 'sub-E*ses-2'; 
                  'sub-D*ses-1'; 'sub-D*ses-2';                   
                  'sub-K*ses-1'; 'sub-K*ses-2'};
        scanlab = {'COMP1'; 'CH1'; 'AC1'; 'KET'; 'COMP2'; 'CH2'; 'AC2'; 'Placebo'};
        TR_path = 'E:\ptsd_ketamine\fmri_scans\';
        test = 'PTSD_KET';
    elseif contains(test, 'PAIN_KET')
        tokens = {'sub-*ses-1'; 'sub-*ses-2'};
        scanlab = {'TP1', 'TP2'};
        TR = 2.5;
        test = 'PAIN_KET';
    elseif and(contains(test, 'JOY'), ~contains(test, 'DOUBLE'))
        tokens = {'sub-'};
        scanlab = {'all'};
        TR = 2;
        test = '';
        PAIRED_RUNS = false;
    elseif and(contains(test, 'JOY'), contains(test, 'DOUBLE'))
        tokens = {'sub-*ses-1'; 'sub-*ses-4'};
        scanlab = {'TP1', 'TP4'};
        TR = 2;
        test = 'TP1_TP4';
    end

