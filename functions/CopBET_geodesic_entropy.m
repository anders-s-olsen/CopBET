% out = CopBET_geodesic_entropy(in,varargin)
%
% Copenhagen Brain Entropy Toolbox: Geodesic entropy
% Calculates Geodesic entropy as in Viol et al., 2019.For each session, the
% correlation coefficient matrix is established. Then, the matrix is binarised
% according a large range of absolute value thresholds and the mean degree
% is calculated for each threshold. Then, for a range of desired mean
% degrees, the path length distribution is computed for each ROI. Potential 
% infinite path lengths and the diagonal are disregarded. The entropy of
% the path length distribution is then evaluated for each ROI. The
% characteristic path length is returned as the mean of ROI entropies. 
% 
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells), e.g., different subjects or scan sessions.
%   
%   varargin (name-value pairs):
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%   parallel: Whether to run temporal entropy in parallel (true)
%
%
% Neurobiology Research Unit, 2023
% Please cite McCulloch, Olsen et al., 2023: "Navigating Chaos in
% Psychedelic Neuroimaging: A Rigorous Empirical Evaluation of the Entropic
% Brain Hypothesis" if you use CopBET in your studies. Please read the
% paper to get a notion of our recommendations regarding the use of the
% specific methodologies in the toolbox.

% ASO March-April 2023

function out = CopBET_geodesic_entropy(in,varargin)

[out,numworkers,in] = CopBET_function_init(in,varargin);

% Select ranges and initialize variables
opts.thresholds = 0.01:0.001:0.99; %threshold values
opts.chosen_degrees = 1:100; %degrees across which to test
opts.bins = 0.5:100.5; %histogram bins

%load data
parfor(ses = 1:height(in),numworkers)
    disp(['Working on entropy calculations for session: ',num2str(ses)])
    ts = in{ses,1}{1};
    R = corr(ts);
    sensible_data_check(R,['correlation matrix, session ',num2str(ses)]);
    entropy{ses} = calc_viol19_entropy(R,opts);
    sensible_data_check(entropy{ses}(24:39),['entropy, session ',num2str(ses)]);
end

out.entropy = entropy';
end

%% functions

% test graph
% adj_mat = [0,1,1,0,1,1,1,0,0;1,0,1,1,0,0,0,0,0;1,1,0,0,1,0,0,0,0;0,1,0,0,0,0,0,0,0;1,0,1,0,0,0,0,0,0;1,0,0,0,0,0,1,0,0;1,0,0,0,0,1,0,1,1;0,0,0,0,0,0,1,0,1;0,0,0,0,0,0,1,1,0]

function entropy = calc_viol19_entropy(Z,opts)
bin_matrix=cell(1,length(opts.thresholds));
degree = zeros(length(bin_matrix),1);
entropy = zeros(1,length(opts.chosen_degrees));
for tt=1:length(opts.thresholds) %Loop through threshold values
    
    %binarise the matrix, establish graph
    bin_matrix{tt} = abs(Z)>opts.thresholds(tt);
    
    %Calculate mean degree
    degree(tt) = mean(degrees_und(bin_matrix{tt}));
    
end

for i=1:length(opts.chosen_degrees)
    % extract the appropriate matrix
    [~,idx] = min(abs(degree - opts.chosen_degrees(i)));
    G = graph(bin_matrix{idx});
    
    %%%%%% calculate the geodesic distance for all nodes
    entropy_n = nan(1,length(Z));
    paths = distances(G); % can be 0 (self-connection), int, or inf
    
    % inf obscures the probability density estimation. Solution: disregard
    paths(isinf(paths)) = NaN;
    paths(logical(eye(size(paths)))) = NaN; % dont count 0 values in diagonal
    
    for n = 1:length(Z)
        [p_n,~] = histcounts(paths(n,~isnan(paths(n,:))),opts.bins,'Normalization','probability');
        
        if abs(sum(p_n)-1)>0.01
            error('Poor density estimation')
        end
        
        entropy_n(n) = -nansum(p_n.*log(p_n));
    end
    % if all path lengths are inf -->
    entropy_n(entropy_n==0) = NaN;
    % characteristic entropy
    entropy(i) = nanmean(entropy_n);
    
end
end