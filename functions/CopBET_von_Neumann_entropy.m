% out = CopBET_von_Neumann_:entropy(in,keepdata)
%
% Copenhagen Brain Entropy Toolbox: Von-Neumann entropy.
% This function evaluates Von-Neumann entropy as introduced in Felippe et
% al., 2021. The Von-Neumann entropy is a Shannon entropy of a correlation
% matrix scaled by the number of regions. 
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before clustering, e.g.,
%   different subjects or scan sessions.
% name-value pairs:
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%
% Neurobiology Research Unit, 2023
% Please cite McCulloch, Olsen et al., 2023: "Navigating Chaos in
% Psychedelic Neuroimaging: A Rigorous Empirical Evaluation of the Entropic
% Brain Hypothesis" if you use CopBET in your studies. Please read the
% paper to get a notion of our recommendations regarding the use of the
% specific methodologies in the toolbox.

% ASO 9/3-2023

% potential tests:
% Check for nans
% check that entropy values are sensible...

function out = CopBET_von_Neumann_entropy(in,varargin)
[out,numworkers,in] = CopBET_function_init(in,varargin);

%load data
entropy = nan(height(in),1);
for ses = 1:height(in)
    
    R = corrcoef(in{ses,1}{1});
    
    rho = R./size(R,1);
    
%     entropy(ses) = -trace(rho*logm(rho));
    
    eigvals = eig(rho);
    
    % sometimes the eigenvalues are machine-precision negative...
    eigvals(eigvals<0)=eps;
    entropy(ses) = -sum(eigvals.*log(eigvals));
    
    % Normalize by maximum entropy
    entropy(ses) = entropy(ses)./log(size(R,1));
    
end

out.entropy = entropy;
end