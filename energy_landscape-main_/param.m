
function [test, basedir, tokens, atlas, nparc] = param()

    experiment = 'PTSD_KET_Schaefer2018_cort_only'; %NOAM_KETAMINE_Lausanne;

    if contains(experiment, 'PTSD_KET')
        tokens = {'sub-F*ses-1'; 'sub-F*ses-2';
                  'sub-E*ses-1'; 'sub-E*ses-2'; 
                  'sub-D*ses-1'; 'sub-D*ses-2';                   
                  'sub-K*ses-1'; 'sub-K*ses-2'};
        test = 'PTSD_KET';
        basedir = 'E:\ptsd_ketamine\';
    elseif contains(experiment, 'PTSD_ONLY')
        tokens = {
          'sub-F*ses-1'; 'sub-F*ses-2';
          'sub-E*ses-1'; 'sub-E*ses-2'; 
          'sub-D*ses-1'; 'sub-D*ses-2';                   
          };
        test = 'PTSD_ONLY';
        basedir = 'E:\ptsd_ketamine\';
    elseif contains(experiment, 'KET_ONLY')
        tokens = {'sub-K*ses-1'; 'sub-K*ses-2'};
        test = 'KET_ONLY';
        basedir = 'E:\ptsd_ketamine\';
    end

    if contains(experiment, 'Schaefer2018')
        basedir = [basedir, 'Results_Schaefer2018_7Networks\'];
        atlas = 'E:\ptsd_ketamine\atlas\SchafferTian-Yeo.xlsx';
        nparc = 454;
    elseif contains(experiment, 'Lausanne')
        basedir = [basedir, 'Results_Lausanne\Raw'];
        atlas = 'C:\Users\marko\ketamine\energy_landscape-main\data\Lausanne_463_subnetworks.mat';
        nparc = 463;
    end
end
