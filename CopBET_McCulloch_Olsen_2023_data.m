%%
function [tbl,data,opts] = CopBET_McCulloch_Olsen_2023_data(atlas,ts_ROI2ROI,extravar,extravartype)

if nargin==0||isempty(atlas)
    % locate and load data sheet (grab latest)
    datafile = dir(['/mrdata/np2/p3/entropy/critical_files/data_*_Z.mat']);
    [~,latest_datasheet] = max(datetime({datafile.date}));
    load([datafile(latest_datasheet).folder,'/',datafile(latest_datasheet).name]);
    disp(['Loaded ',datafile(latest_datasheet).name])
    atlas = [];
else
    %%%%%%%%%%%%%% function-wide settings and import stuff
    possible_atlases = {'aal90','Craddock181','HarvardOxford105',...
        'Lausanne462','Schaefer1000','Shen218','yeo17','smith20','SchaeferTian232'};
    
    if all(~strcmpi(atlas,possible_atlases))
        disp(possible_atlases)
        error('Please input a different atlas name. Possible options above')
        
    end
    
    % locate and load data sheet (grab latest if more than 1)
    if nargin<2
        error('Please specify either ''ts'' or ''ROI2ROI'' as 2nd argument')
    end
    if strcmp(ts_ROI2ROI,'ts')
        datafile = dir(['/mrdata/np2/p3/entropy/critical_files/data_*',atlas,'_ts.mat']);
    elseif strcmp(ts_ROI2ROI,'ROI2ROI')
        datafile = dir(['/mrdata/np2/p3/entropy/critical_files/data_*',atlas,'_Z.mat']);
    else
        error('Please specify data sheet as either ''ts'' or ''ROI2ROI''')
    end
    if numel(datafile)>1
        [~,latest_datasheet] = max(datetime({datafile.date}));
        load([datafile(latest_datasheet).folder,'/',datafile(latest_datasheet).name]);
        disp(['Processing ',datafile(latest_datasheet).name])
    else
        load([datafile.folder,'/',datafile.name]);
        disp(['Processing ',datafile.name])
    end
end

data = SDI_PPL_rsfmri_tbl;

if nargin<3
    extravar = {};extravartype = {};
elseif nargin==3
    error('Please also specify variable type')
end

if ~iscell(extravar)||~iscell(extravartype)
    error('Please specify extra variables as cell')
end


% locate correct table column name for data
if any(~cellfun(@isempty,regexp(data.Properties.VariableNames,'ts_')))
    opts.datacol = data.Properties.VariableNames{~cellfun(@isempty,regexp(data.Properties.VariableNames,'ts_'))};
elseif any(~cellfun(@isempty,regexp(data.Properties.VariableNames,'ROI2ROI_')))
    opts.datacol = data.Properties.VariableNames{~cellfun(@isempty,regexp(data.Properties.VariableNames,'ROI2ROI_'))};
end
if strcmp(ts_ROI2ROI,'denoised_volumes')
    opts.datacol = data.Properties.VariableNames{~cellfun(@isempty,regexp(data.Properties.VariableNames,'rsfmri_location'))};
end

opts.subjects = unique(data.CIMBI_ID);
opts.include = strcmp(data.include,'yes');


tblvarnames = [{'data','rp','CIMBI', 'ses','time_since_admin', 'PPL','Occ', 'SDI', 'num_vols','mr','entropy'},extravar];
tblvartypes = [{'cell','cell','double','double','double','double','double','double','double','cell','cell'},extravartype];
tbl  = table('Size',[sum(opts.include),numel(tblvartypes)],...
    'VariableNames', tblvarnames,...
    'VariableTypes',tblvartypes);

opts.nrois = size(data.(opts.datacol){find(opts.include,1)},2);
opts.tssize = 0;
tblcount = 1;
for sub = 1:numel(opts.subjects) %loop through subjects
    
    %fmri indices for this subject
    f_sub = find(data.CIMBI_ID==opts.subjects(sub)&opts.include);
    
    [~,sortidx] = sort((data.Clock_time(f_sub)));
    f_sub_sort = f_sub(sortidx);
    
    opts.sescount = 1;
    for ses = 1:numel(f_sub_sort)
        % Fill tbl
        tbl = allfunc_filltbl(tbl,data,f_sub_sort(ses),tblcount,opts);
        tbl.ses(tblcount) = opts.sescount;
        tbl.sesidx(tblcount) = f_sub_sort(ses);
        if regexp(opts.datacol,'ts_','once')
            tbl.num_vols(tblcount) = size(data.(opts.datacol){f_sub_sort(ses)},1);
            opts.ts_size = opts.tssize + tbl.num_vols(height(tbl));
            
        end
        
        opts.sescount = opts.sescount + 1;
        tblcount = tblcount + 1;
    end
end

mr45idx = find(strcmp(tbl.mr,'mr45'));
mr001idx = find(strcmp(tbl.mr,'mr001'));
tbl = vertcat(tbl(mr45idx,:),tbl(mr001idx,:));


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

pa_out = find(~cellfun(@isempty,regexp(data.rsfmri_location(tbl.sesidx),'_pa','once')));
tbl(pa_out,:) = [];
% data(pa_out,:) = [];


function tbl = allfunc_filltbl(tbl,data,idx,count,opts)
mrID = data.rsfmri_raw_ID{idx};
[~,scan]=fileparts(data.rsfmri_location{idx});

tbl.time_since_admin(count) = data.time_since_admin(idx);
tbl.CIMBI(count) = data.CIMBI_ID(idx);
tbl.PPL(count) = 1020*data.psi_conc(idx);
tbl.Occ(count) = 0.766*tbl.PPL(count)./(1.95+tbl.PPL(count));
tbl.SDI(count) = data.SDI_score(idx);
if ischar(data.(opts.datacol){idx})
    
    tbl.data{count} = ['/mnt/mregdata2/Entropy_denoised_images/denoised_images/denoised',mrID,'_',scan,'.nii'];
    if ~isempty(regexp(data.rsfmri_location{idx},'/mr45','once'))
        tbl.mr{count} = 'mr45';
    elseif ~isempty(regexp(data.rsfmri_location{idx},'/mr001'))
        tbl.mr{count} = 'mr001';
    end
    rp = readtable(['/mrdata/np2/p3/entropy/data/rp_all/rp',mrID,'_',scan,'.txt']);
    tbl.rp{count} = table2array(rp);
else
    if ~isempty(regexp(data.rsfmri_location{idx},'/mr45','once'))
        tbl.mr{count} = 'mr45';
        tbl.data{count} = data.(opts.datacol){idx};
    elseif ~isempty(regexp(data.rsfmri_location{idx},'/mr001'))
        tbl.mr{count} = 'mr001';
        tbl.data{count} = resample(data.(opts.datacol){idx},2,5); %resample
    end
    sensible_data_check(tbl.data{count});
end



