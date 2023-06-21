% out = CopBET_temporal_entropy(in,TR,keepdata,parallel)
%
% Copenhagen Brain Entropy Toolbox: Temporal entropy.
% Evaluates temporal entropy as in Luppi et al., 2021. Sliding window
% connectivity matrices are constructed, The Louvain community detection
% algorithm is run for each matrix, and the module degree z-score and
% participation coefficient is evaluated. Based on these (concatenated
% across all windows), a K=2 kmeans algorithm is run on the cartographic
% profile. This presumably generates an integrative and segregated state,
% of which the entropy of the activity profile is evalated. All the above
% is performed for each subject, including a k-means with 500 replications,
% so this code takes forever to run. 
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before clustering, e.g.,
%   different subjects or scan sessions.
%   TR: TR for constructing tapered windows
% name-value pairs:
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
% Check that the window looks sensible (qualitatively)
% check for nans all over
% Check that the k-means output looks sensible
% check for within-reasonable-range entropy values

function out = CopBET_temporal_entropy(in,TR,varargin)

if nargin<2
    error('Please specify the TR')
end

[out,numworkers,in] = CopBET_function_init(in,varargin);

window = construct_tapered_window(TR);

for ses = 1:height(in)
    disp(['Working on entropy calculations for session: ',num2str(ses)])
    tmp = in{ses,1}{1};
    tmp = tmp - mean(tmp);
    entropy(ses) = Luppi21_entropy(tmp,window,numworkers);
    
end
out.entropy = entropy';
end
%% functions

function window = construct_tapered_window(TR)

time_axis = 0:TR:22*TR; %This would start at TR if window size == 22TRs
if mod(numel(time_axis),2)==0
    time_axis(1) = [];
end

% construct rectangular and gaussian windows and convolve them
rect_win = ones(numel(time_axis),1);
gauss_win = gausswin(numel(time_axis),3*TR);%FWHM 3 TRs

window = conv(rect_win,gauss_win,'same');
window = window/sum(window); %make elements sum to 1

end

function entropy = Luppi21_entropy(data,window,numworkers)
% Much of the code below was given to us by Andrea Luppi


nNodes = size(data,2);

% Construct all window indices based on the input window
window_size_oneside = (numel(window)-1)/2;
end_point = size(data,1)-window_size_oneside;
window_midpoints = window_size_oneside+1:end_point;

% Preallocate matrices
dFCmat = nan(nNodes,nNodes,numel(window_midpoints));
ci = nan(nNodes,numel(window_midpoints));
q = nan(numel(window_midpoints),1);
WT = nan(nNodes,numel(window_midpoints));
BT = nan(nNodes,numel(window_midpoints));


% Loop through all windows (sliding across the data)
parfor (i = 1:numel(window_midpoints),numworkers)
% for i=1
    ci_iter = nan(nNodes, 1);
    q_iter = nan(100, 1);
    full_window = zeros(1,size(data,1));
    full_window(window_midpoints(i)-window_size_oneside:window_midpoints(i)+window_size_oneside) = window;
    
    data_win = data.*full_window';
    data_win(sum(data_win,2)==0,:) = [];

    dFCmat = corr(data_win);
    
    for iteration = 1:100
        
        ci_iter(:, 1)  = 1:nNodes;                   % initial community affiliations
        Q0 = -1;
        q_iter(iteration) = 0;            % initialize modularity values
        while q_iter(iteration)-Q0>1e-5;           % while modularity increases
            Q0 = q_iter(iteration);                % perform community detection

            %use improved version of louvain algorithm from BCT
            [ci_iter(:, 1),q_iter(iteration)] = community_louvain(dFCmat, [], ci_iter(:,1), 'negative_asym');
        end
        
        %update the timepoint-specific values based on whether current
        %iteration is better
        if iteration == 1
            q(i, 1) = q_iter(iteration);
            ci(:, i) = ci_iter(:, 1);
        else
            if q_iter > q(i, 1)
                q(i, 1) = q_iter(iteration);
                ci(:, i) = ci_iter(:, 1);
            end
        end    
        
    end
    
    WT(:,i) = module_degree_zscore(dFCmat,ci(:,i),0);
    BT(:,i) = participation_coef_sign(dFCmat,ci(:,i));
    
%     disp(['Done with time point ',num2str(i),' of ',num2str(numel(window_midpoints))])
    
end

% Step 4: 2-dimensional Cartographic Profile (CP)

xbins = [0:0.01:1.0]; ybins = [5:-.1:-5]; % 100 x 100 2d histogram
CP = zeros(size(xbins,2),size(ybins,2),numel(window_midpoints));
xNumBins = numel(xbins); yNumBins = numel(ybins);

for t = 1:numel(window_midpoints)
  Xi = round(interp1(xbins, 1:xNumBins, BT(:,t), 'linear', 'extrap') );
  Yi = round(interp1(ybins, 1:yNumBins, WT(:,t), 'linear', 'extrap') );
  Xi = max( min(Xi,xNumBins), 1);
  Yi = max( min(Yi,yNumBins), 1);
  CP(:,:,t) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
end

if numworkers>0
    idx=kmeans(reshape(CP,xNumBins * yNumBins,numel(window_midpoints))',2,'Distance','correlation',...
    'Replicates',500,'Display','off','Options',statset('UseParallel',1));
else
    idx=kmeans(reshape(CP,xNumBins * yNumBins,numel(window_midpoints))',2,'Distance','correlation',...
    'Replicates',500,'Display','off');
    
end
% This below gives the same as Luppi's version
[Prob,bins] = histcounts(idx,[0.5,1.5,2.5],'Normalization','probability');

entropy = -sum(Prob.*log2(Prob));

end

%% Functions from the 2013 BCT toolbox


function Z=module_degree_zscore(W,Ci,flag)
%MODULE_DEGREE_ZSCORE       Within-module degree z-score
%
%   Z=module_degree_zscore(W,Ci,flag);
%
%   The within-module degree z-score is a within-module version of degree
%   centrality.
%
%   Inputs:     W,      binary/weighted, directed/undirected connection matrix
%               Ci,     community affiliation vector
%               flag,   0, undirected graph (default)
%                       1, directed graph: out-degree
%                       2, directed graph: in-degree
%                       3, directed graph: out-degree and in-degree
%
%   Output:     Z,      within-module degree z-score.
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   Mika Rubinov, UNSW, 2008-2010

if ~exist('flag','var')
    flag=0;
end

switch flag
    case 0; % no action required
    case 1; % no action required
    case 2; W=W.';
    case 3; W=W+W.';
end

n=length(W);                        %number of vertices
Z=zeros(n,1);
for i=1:max(Ci)
    Koi=sum(W(Ci==i,Ci==i),2);
    Z(Ci==i)=(Koi-mean(Koi))./std(Koi);
end

Z(isnan(Z))=0;
end

function [Ppos Pneg]=participation_coef_sign(W,Ci)
%PARTICIPATION_COEF_SIGN     Participation coefficient
%
%   [Ppos Pneg] = participation_coef_sign(W,Ci);
%
%   Participation coefficient is a measure of diversity of intermodular
%   connections of individual nodes.
%
%   Inputs:     W,      undirected connection matrix with positive and
%                       negative weights
%
%               Ci,     community affiliation vector
%
%   Output:     Ppos,   participation coefficient from positive weights
%
%               Pneg,   participation coefficient from negative weights
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   2011, Mika Rubinov, UNSW

%   Modification History:
%   Mar 2011: Original
%   Sep 2012: Fixed treatment of nodes with no negative strength
%             (thanks to Alex Fornito and Martin Monti)


n=length(W);                                %number of vertices

Ppos = pcoef( W.*(W>0));
Pneg = pcoef(-W.*(W<0));

    function P=pcoef(W_)
        S   = sum(W_,2);                    %strength
        Gc  = (W_~=0)*diag(Ci);             %neighbor community affiliation
        Sc2 = zeros(n,1);                   %community-specific neighbors

        for i = 1:max(Ci);
            Sc2 = Sc2 + (sum(W_.*(Gc==i),2).^2);
        end

        P = ones(n,1) - Sc2./(S.^2);
        P(isnan(P)) = 0;
        P(~P) = 0;                            %p_ind=0 if no (out)neighbors
    end
end




