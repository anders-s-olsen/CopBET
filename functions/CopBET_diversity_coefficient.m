% out = CopBET_diversity_coefficient(in,varargin)
%
% Copenhagen Brain Entropy Toolbox: Diversity coefficient
% Evaluates diversity coefficient as in Lebedev et al., 2015. First, the
% correlation coefficient matrix is established for all sessions. Then, the
% modularity_finetune_und algorithm is run 1000 times on the average
% (across sessions) correlation matrix, where the modularity assignment
% corresponding to the highest modularity score is extracted. For each
% session, the diversity coefficient (a measure of the diversity of the
% out-network connectivity, i.e., the rois NOT part of the same module) is
% evaluated using diversity_coef_sign. Only Hpos is reported (see BCT
% function for reference). 
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before modularity, e.g.,
%   different subjects or scan sessions.
%   
%   varargin (name-value pairs):
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%
%
% Neurobiology Research Unit, 2023
% Please cite McCulloch, Olsen et al., 2023: "Navigating Chaos in
% Psychedelic Neuroimaging: A Rigorous Empirical Evaluation of the Entropic
% Brain Hypothesis" if you use CopBET in your studies. Please read the
% paper to get a notion of our recommendations regarding the use of the
% specific methodologies in the toolbox.

% ASO March-April 2023

function out = CopBET_diversity_coefficient(in,varargin)

[out,~,in] = CopBET_function_init(in,varargin);

nrois = size(in{1,1}{1},2);
R_all = nan(nrois,nrois,height(in));
R_sub = cell(1,height(in));
for ses = 1:height(in)
    ts = in{ses,1}{1};
    R = corr(ts);
    sensible_data_check(R,['correlation matrix, session ',num2str(ses)]);
    R_all(:,:,ses) = R;
    R_sub{ses} = R;
end

optimal_assignment = optimal_modularity(R_all);

% Calculate diversity coefficient for every connectivity matrix
% (Lebedev reports Hpos)
for ses = 1:height(in)
    [Div_coef,~] = diversity_coef_sign(R_sub{ses}, optimal_assignment);
    entropy{ses} = Div_coef;
    sensible_data_check(Div_coef,['diversity coefficient, session ',num2str(ses)]);
end

out.entropy = entropy';


end

%% functions
% Unsure whether this corresponds to the actual algorithm used by Lebedev.
% However, we know that the assignment were performed on an average
% correlation matrix from placebo scans only. At some point on Alexander's
% github there is a modularity_louvain_und_sign algorithm, suggesting he
% used a binary matrix, but this was not explained in the paper. 
%
% We also have no way of knowing the general range of Hpos values as the
% figure in the paper is of too low resolution. 
function optimal_assignment = optimal_modularity(Z_all)

% Run modularity algorithm on an average connectivity matrix
wei_matrix_avg = mean(Z_all,3);

iterations = 1:1000;
mod_communities = cell(length(iterations), 1);
mod_score = nan(1,length(iterations));

num_rois  = size(wei_matrix_avg,1);
for xx=1:length(iterations)
    [M1,Q1]=modularity_finetune_und(wei_matrix_avg,1:num_rois,1);
    mod_communities{xx} = M1;
    mod_score(xx) = Q1;
end

% Extract best instance
[~,idx]=max(mod_score);
optimal_assignment = mod_communities{idx};
end

