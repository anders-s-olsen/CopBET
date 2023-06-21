clear
addpath('/mrdata/np2/p3/entropy/CopBET/functions')

atlasnames = {'HarvardOxford105','aal90','yeo17','Schaefer1000','Shen218',...
    'Craddock181','Lausanne462','Smith20','yeo7','SchaeferTian232'};

for atl = 1:numel(atlasnames)
    for test = 1:2
        if test==1
            load(['/mrdata/np2/p3/entropy/LSDdata/ROIdata/',atlasnames{atl},'/sub-001_ses-LSD_task-rest_run-01_bold'])
            data = V_roi;
            run_functions2(['/mrdata/np2/p3/entropy/LSDdata/sub-001/ses-LSD/func/sub-001_ses-LSD_task-rest_run-01_bold.nii.gz']);
            run_functions1(data);
        elseif test==2
            load(['/mrdata/np2/p3/entropy/LSDdata/ROIdata/',atlasnames{atl},'/sub-001_ses-LSD_task-rest_run-03_bold'])
            tbl.in{1} = data; %previous
            tbl.in{2} = V_roi;%current
            tbl2.in{1} = ['/mrdata/np2/p3/entropy/LSDdata/sub-001/ses-LSD/func/sub-001_ses-LSD_task-rest_run-01_bold.nii.gz'];
            tbl2.in{2} = ['/mrdata/np2/p3/entropy/LSDdata/sub-001/ses-LSD/func/sub-001_ses-LSD_task-rest_run-03_bold.nii.gz'];
            run_functions2(tbl);
            run_functions1(tbl);
        end
    end
end


function run_functions1(in)
% out = CopBET_DCC_entropy(in);
% output_check(out)

% out = CopBET_degree_distribution_entropy(in);
% output_check(out)

% out = CopBET_diversity_coefficient(in);
% output_check(out)

% out = CopBET_geodesic_entropy(in);
% output_check(out)

out = CopBET_LEiDA_transition_entropy(in,3);
output_check(out)

out = CopBET_metastate_series_complexity(in);
output_check(out)

out = CopBET_motif_connectivity_entropy(in);
output_check(out)

out = CopBET_sample_entropy(in);
output_check(out)

out = CopBET_temporal_entropy(in);
output_check(out)

out = CopBET_time_series_complexity(in);
output_check(out)

out = CopBET_von_Neumann_entropy(in);
output_check(out)
end

function run_functions2(in)
% smith_atlas = niftiread('/mrdata/np2/p3/entropy/ROIs/smith20_thresholded_2mm.nii');
% smith_atlas = logical(smith_atlas(:,:,:,[1,2,5,6,9,11,12,14,15]+1));
% out = CopBET_intranetwork_synchrony(in,smith_atlas);
% output_check(out)
end

function output_check(out)
for ses = 1:height(out)
    if any(isnan(out.entropy{ses}))
        warning('NaN produced')
    end
    if any(isinf(out.entropy{ses}))
        warning('inf produced')
    end
    if any(out.entropy{ses}==0)
        warning('0 produced')
    end
end
end