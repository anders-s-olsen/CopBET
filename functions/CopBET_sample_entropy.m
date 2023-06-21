% out = CopBET_degree_distribution_entropy(in,keepdata,parallel)
%
% Copenhagen Brain Entropy Toolbox: Temporal entropy
% Evaluates temporal entropy as in Luppi et al., 2021. Sliding window
% connectivity matrices are constructed, The Louvain community detection
% algorithm is run for each matrix, and the module degree z-score and
% participation coefficient is evaluated. Based on these (concatenated
% across all windows), a K=2 kmeans algorithm is run on the cartographic
% profile. This presumably generates an integrative and segregated state,
% of which the entropy of the activity profile is evalated. All the above
% is performed for each subject, including a k-means with 500 replications,
% so this code takes forver to run. 
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before clustering, e.g.,
%   different subjects or scan sessions.
%   TR: TR for constructing tapered windows
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%   parallel: Whether to run temporal entropy in parallel
%
%
% Neurobiology Research Unit, 2023
% Please cite McCulloch, Olsen et al., 2023: "Navigating Chaos in
% Psychedelic Neuroimaging: A Rigorous Empirical Evaluation of the Entropic
% Brain Hypothesis" if you use CopBET in your studies. Please read the
% paper to get a notion of our recommendations regarding the use of the
% specific methodologies in the toolbox.

% ASO 9/3-2023

% potential tests:
% Viol mentioned some tests in her email

function out = CopBET_sample_entropy(in,atlas,keepdata,parallel,NRUspecific,doitornot)

if nargin<3
    keepdata = true;
    parallel = true;
    NRUspecific = false;
elseif nargin < 4
    parallel = true;
    NRUspecific = false;
elseif nargin<5
    NRUspecific = false;
elseif nargin<2
    error('Please specify both input data and an atlas (3D)')
end
if keepdata
    if any(strcmp(in.Properties.VariableNames,'entropy'))
        warning('Overwriting entropy column in datatbl')
    end
end

if ~istable(in)
    if isstr(in)
        % convert matrix to table with one entry
        tbl = table;
        tbl.in{1} = in;
        in = tbl;
    else
        error(['Please specify the input data as either a matrix (nxp, n>1)', ...
            'or a table of matrices tbl where the FIRST column contains the data',...
            'with a matrix for each row'])
    end
end

if parallel
    numworkers = 10;
else
    numworkers = 0;
end

brainVox = find(atlas(:)>0);
imgSize = size(atlas);
num_rois = numel(unique(atlas(atlas>0)));

r = 0.3;
m = 2;
scale = 1:5;

%load data
% parfor(ses = 1:height(in),numworkers)
for ses = 1:height(in)
    disp(['Working on entropy calculations for session: ',num2str(ses)])
    % 4D file
    if ~doitornot
        for a = 1:numel(scale)
        V = niftiread(['/mrdata/np2/p3/entropy/critical_files/Lebedev16_output_v3/',...
            num2str(in.sesidx(ses)),'_scale',num2str(a)]);
        V = double(V(:));
        
        for roi = 1:num_rois
            tmp = V(atlas(:)==roi);
            entropy{ses}(a,roi) = mean(tmp(tmp~=0));
        end
        end
    else
    path = in{ses,1}{1};
    image_4D = double(niftiread(path)); %4D series
    if NRUspecific
    if ~isempty(regexp(path,'mr001'))
        image_4D = NRUspecific_downsamplemr001data(image_4D);
    end
    end
    
    image_4D = image_4D - mean(image_4D,4);
    
    imsz = size(image_4D);
    
    if any(imsz(1:3)~=imgSize)
        error('Wrong image size')
    end
    
    
    % calculate sample entropy for all voxels
    for a = 1:numel(scale)
        MSE = nan(imgSize(1)*imgSize(2)*imgSize(3),1);
        parfor (vox = 1:length(brainVox),numworkers)
%         for vox = 1:length(brainVox)
            [row,col,sl] = ind2sub(imgSize,brainVox(vox));
            ts = squeeze(image_4D(row,col,sl,:));
            
            r_val = r*std(double(ts));
            tmp = sample_entropy(m,r_val,ts,scale(a));
            MSE(vox) = tmp(1);
            
            %         if ismember(vox,[10000:10000:numel(brainVox)])
            %             disp(['Done with ',num2str(vox),' of ',num2str(numel(brainVox))])%,' in ',num2str(toc),' seconds'
            %         end
            
        end
        MSE2 = nan(imgSize);
        for vox = 1:numel(brainVox)
            [row,col,sl] = ind2sub(imgSize,brainVox(vox));
            MSE2(row,col,sl) = MSE(vox);
        end
        
        for roi = 1:num_rois
            tmp = MSE2(atlas==roi);
            entropy{ses}(scale,roi) = mean(tmp(tmp~=0));
        end
        
    end
    
    end
    
end
if keepdata
    out = in;
    out.entropy = entropy';
else
    out = table;
    out.entropy = entropy';
end
end

