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

function out = CopBET_motif_connectivity_entropy(in,rp,atlas,TR,keepdata,parallel,NRUspecific)

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
        if isempty(rp)
            error('need to specify motion correction time series as well')
        end
        tbl.rp{1} = rp;
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

window_lengths = 15:1:150; %seconds

%load data
parfor(ses = 1:height(in),numworkers)
% for ses = 1:height(in)
    disp(['Working on entropy calculations for session: ',num2str(ses)])
    % 4D file
    path = in{ses,1}{1};
    rp = in{ses,2}{1};
    
    ROImeans = load_data(path,atlas,NRUspecific);
    rp = fwd_calc(rp);
    
    state_hist1 = state_dis(ROImeans,rp,TR,window_lengths);
    entropy{ses} = calc_shannon_ent(state_hist1);
    
end
if keepdata
    out = in;
    out.entropy = entropy';
else
    out = table;
    out.entropy = entropy';
end
end

%% functions
function ROImeans = load_data(path,atlas,NRUspecific)

image_4D = double(niftiread(path)); %4D series
if NRUspecific
    if ~isempty(regexp(path,'denoisedn'))
        image_4D = NRUspecific_downsamplemr001data(image_4D);
    end
end
datasz = size(image_4D);
num_rois = numel(unique(atlas(atlas>0)));
ROImeans = zeros(num_rois,datasz(4));

if any(datasz(1:3)~=size(atlas))
    error('wrong atlas size')
end

for vol = 1:datasz(4)
    v = image_4D(:,:,:,vol);
    for ROI = 1:num_rois
        ROImeans(ROI,vol) = mean(v(atlas==ROI));
    end
end

if any(isnan(ROImeans(:)))
    error(['Problem with nan for data '])
end

ROImeans = single(ROImeans);
ROImeans = (ROImeans - mean(ROImeans,2))./std(ROImeans,[],2);
end

function [word_count1] = state_dis(data,rp,TR,window_lengths)

possible_graphs = dec2bin(0:2^6-1)' - '0';
ROIpairs = [1,2;1,3;1,4;2,3;2,4;3,4];

eff_window_length = round(window_lengths/TR);

word_count1 = zeros(numel(window_lengths),size(possible_graphs,2));
% word_count2 = zeros(numel(window_lengths),size(possible_graphs,2));
for wl = 1:numel(window_lengths)
    
    if wl>1
        if eff_window_length(wl)==eff_window_length(wl-1)
            word_count1(wl,:) = word_count1(wl-1,:);
%             word_count2(wl,:) = word_count2(wl-1,:);
            continue
        end
    end
    
    window_starts = [0,eff_window_length(wl):eff_window_length(wl):size(data,2)]+1;
    
    % throw last window away
    if size(data,2)-window_starts(end) < eff_window_length(wl)
        window_starts(end)=[];
    end
    
    
    for w = 1:numel(window_starts)
        
        datatmp = data(:,window_starts(w):window_starts(w)+eff_window_length(wl)-1);
        rptmp = abs(rp(window_starts(w):window_starts(w)+eff_window_length(wl)-1,:));
        
        poss = 1:4;
        pval1 = nan(1,size(ROIpairs,1));
        for pairs = 1:size(ROIpairs,1)
            z = [datatmp(poss(poss~=ROIpairs(pairs,1)&poss~=ROIpairs(pairs,2)),:);rptmp']';
            [RHO1,pval1(pairs)] = partialcorr(datatmp(ROIpairs(pairs,1),:)',datatmp(ROIpairs(pairs,2),:)',z,'type','pearson','tail','both');
        end
        if any(isnan(pval1))
            error('NaN produced')
        end
        
        graph1 = pval1 < 0.05/6;
        
        [val,word1]=min(sum(abs(possible_graphs-graph1')));
        if val ~=0
            keyboard
        end
        
        word_count1(wl,word1) = word_count1(wl,word1)+1;
    end
    
    
end
end

function ent = calc_shannon_ent(word_count)
ent = nan(size(word_count,1),1);
for wl = 1:size(word_count,1)
    state_prob =word_count(wl,:)./sum(word_count(wl,:));
    ent(wl) = sum(state_prob(state_prob>0).*log2(1./state_prob(state_prob>0)));
end

end

function [fwd] = fwd_calc(rp)
radius = 50; %mm
ts = rp;
temp = rp(:,4:6)*180/pi; %convert to degrees
temp = (2*radius*pi/360)*temp;
ts(:,4:6) = temp;
dts = diff(ts);
fwd = [0;sum(abs(dts),2)];

fwd = (fwd - mean(fwd))/std(fwd); %standardize
end
