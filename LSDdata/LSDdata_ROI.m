clear
ROIpath = [pwd,'/Atlases'];
atlaslocs = {'CONN_atlas_2mm.nii',...
    'AAL90_2mm.nii',...
    'Yeo17_liberalmask_2mm.nii',...
    'Yeo7_2mm.nii',...
    'Schaefer_1000_2mm',...
    'Shen268_2mm.nii',...
    'Craddock200_2mm.nii',...
    'Lausanne463_2mm.nii',...
    'smith20_thresholded_3mm.nii',...
    'Schaefer_Tian_232_2mm.nii'};
atlasnames = {'HarvardOxford105','AAL90','yeo17','yeo7','Schaefer1000','Shen268',...
    'Craddock200','Lausanne463','Smith20','SchaeferTian232'};

subjects = dir('LSDdata/sub-*');
conditions = {'ses-LSD','ses-PLCB'};



for i = 1:numel(subjects)
    for j = 1:numel(conditions)
        for run = [1,3] %2 is music
            
            V = niftiread([subjects(i).folder,'/',subjects(i).name,'/',...
                conditions{j},'/func/',subjects(i).name,'_',conditions{j},...
                '_task-rest_run-0',num2str(run),'_bold.nii.gz']);
            V_sz = size(V);
            
            
            if V_sz(4)~=217
                warning('Wrong number of volumes')
            end
            
            for atlas = 1:numel(atlaslocs)
                ROIatlas = niftiread([ROIpath,atlaslocs{atlas}]);
                ROIsz = size(ROIatlas);
                
                if  ~isequal(V_sz(1:3),ROIsz(1:3))
                    error('wrong size')
                end
                
                
                if ndims(ROIatlas)==3
                    ROIs = unique(ROIatlas(ROIatlas>0));
                    num_rois = numel(ROIs); 
                    V_roi = nan(V_sz(4),num_rois);
                    for roi = 1:num_rois
                        for t = 1:V_sz(4)
                            tmp = V(:,:,:,t);
                            tmp2 = tmp(ROIatlas==ROIs(roi));
                            V_roi(t,roi) = nanmean(tmp2(tmp2~=0),'all');
                        end
                    end
                elseif ndims(ROIatlas)==4
                    num_rois = size(ROIatlas,4);
                    V_roi = nan(V_sz(4),num_rois);
                    for roi = 1:num_rois
                        for t = 1:V_sz(4)
                            tmp = V(:,:,:,t);
                            tmp2 = tmp(ROIatlas(:,:,:,roi));
                            V_roi(t,roi) = nanmean(tmp2(tmp2~=0),'all');
                        end
                    end
                end
                
                mkdir([pwd,'/LSDdata/ROIdata/',atlasnames{atlas}])
                save([pwd,'/ROIdata/',atlasnames{atlas},...
                    '/',subjects(i).name,'_',conditions{j},...
                    '_task-rest_run-0',num2str(run),'_bold.mat'],'V_roi')
                
                
            end
            disp(['Done with ',subjects(i).name,'_',conditions{j},...
                    '_task-rest_run-0',num2str(run),'_bold.mat'])
        end
    end
end

%%