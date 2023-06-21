%%
function [tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,ts_ROI2ROI)

if nargin==0||isempty(atlas)
%     error('An atlas needs to be specified')
    atlas = 'yeo7';
else
    %%%%%%%%%%%%%% function-wide settings and import stuff
    possible_atlases = {'AAL90','Craddock200','HarvardOxford105',...
        'Lausanne463','Schaefer1000','Shen268','yeo17','smith20','SchaeferTian232'};
    
    if all(~strcmpi(atlas,possible_atlases))
        disp(possible_atlases)
        error('Please input a different atlas name. Possible options above')
        
    end
    
end

subs = dir([pwd,'/LSDdata/sub-*']);
if isempty(subs)
    error('Please make sure you''re standing in the right directory')
end
nummats = numel(dir([pwd,'/LSDdata/ROIdata/AAL90/*.mat']));
conditions = {'ses-PLCB','ses-LSD'};

opts.subjects = {subs(:).name};

tblvarnames = [{'data','rp','subject', 'condition','session','num_vols','entropy'}];
tblvartypes = [{'cell','cell','cell','cell','double','double','cell'}];
tbl  = table('Size',[nummats,numel(tblvartypes)],...
    'VariableNames', tblvarnames,...
    'VariableTypes',tblvartypes);

tblcount = 1;

for sub = 1:numel(opts.subjects) %loop through subjects
    for cond = 1:numel(conditions)
        for ses = [1,3]
            % Fill tbl
            if strcmp(ts_ROI2ROI,'denoised_volumes')
                tbl.data{tblcount} = [subs(sub).folder,'/',subs(sub).name,'/',...
                conditions{cond},'/func/',subs(sub).name,'_',conditions{cond},...
                '_task-rest_run-0',num2str(ses),'_bold.nii.gz'];
            else
            load([pwd,'/LSDdata/ROIdata/',...
                atlas,'/',subs(sub).name,'_',conditions{cond},...
                '_task-rest_run-0',num2str(ses),'_bold']);
            tbl.data{tblcount} = V_roi;
            tbl.num_vols(tblcount) = size(tbl.data{tblcount},1);
            end
            tbl.subject{tblcount} = subs(sub).name;
            tbl.condition{tblcount} = conditions{cond};
            tbl.session(tblcount) = ses;
            
            tblcount = tblcount + 1;
        end
    end
end

% check?
if ~strcmp(ts_ROI2ROI,'denoised_volumes')
for h = 1:height(tbl)
    for h2 = h+1:height(tbl)
        if size(tbl.data{h},1)==size(tbl.data{h2},1)
            if norm(tbl.data{h}-tbl.data{h2},'fro')<1
                error(['equality problem for ',num2str(h),'-',num2str(h2)])
            end
        end
    end
end
end
data=[];
opts = [];
end

