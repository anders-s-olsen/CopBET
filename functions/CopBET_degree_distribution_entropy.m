% out = CopBET_degree_distribution_entropy(in,varargin)
%
% Copenhagen Brain Entropy Toolbox: Degree distribution entropy
% Evaluates temporal entropy as in Viol et al., 2017. For each session, the
% correlation coefficient matrix is established, where coefficients with a
% non-significant p-value are set to zero. Then, the matrix is binarised
% according a large range of absolute value thresholds and the mean degree
% is calculated for each threshold. Then, for a range of desired mean
% degrees, the degree distribution of the corresponding binarised
% matrix is calculated, and finally, the Shannon entropy is returned. 
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

% Copyright (C) 2023 Anders Stevnhoved Olsen & Drummond E-Wen McCulloch
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

function out = CopBET_degree_distribution_entropy(in,varargin)

[out,numworkers,in] = CopBET_function_init(in,varargin);

% Select ranges and initialize variables
opts.thresh=0.01:0.001:0.99; %threshold values
opts.chosen_degrees = 1:100; %degrees across which to test
opts.bins = 0.5:100.5; %histogram bins

%load data
parfor(ses = 1:height(in),numworkers)
    disp(['Working on entropy calculations for session: ',num2str(ses)])
    ts = in{ses,1}{1};
    [R,PVAL] = corr(ts);
    sensible_data_check(R,['correlation matrix, session ',num2str(ses)]);
    R(PVAL>0.05) = 0;
    entropy{ses} = calc_viol17_entropy(R,opts);
    sensible_data_check(entropy{ses},['entropy, session ',num2str(ses)]);
    
end
out.entropy = entropy';
end

%% functions

% test graph
% graph = [0,1,1,0,1,1,1,0,0;1,0,1,1,0,0,0,0,0;1,1,0,0,1,0,0,0,0;0,1,0,0,0,0,0,0,0;1,0,1,0,0,0,0,0,0;1,0,0,0,0,0,1,0,0;1,0,0,0,0,1,0,1,1;0,0,0,0,0,0,1,0,1;0,0,0,0,0,0,1,1,0]

function entropy_out = calc_viol17_entropy(Z,opts)
bin_matrix=cell(1,length(opts.thresh));
meandegree = zeros(length(bin_matrix),1);
entropy_out = zeros(length(opts.chosen_degrees),1);
for tt=1:length(opts.thresh) %Loop through threshold values
    
    %binarise the matrix
    bin_matrix{tt}=abs(Z)>opts.thresh(tt);
    
    %Calculate mean degree
    meandegree(tt) = mean(degrees_und(bin_matrix{tt}));
    
end

%Extract the matrix with each appropriate degree
for i=1:length(opts.chosen_degrees)
    [~,idx] = min(abs(meandegree - opts.chosen_degrees(i)));
    [Prob,~] = histcounts(degrees_und(bin_matrix{idx}),opts.bins,'Normalization','probability');
    %calculate the Shannon entropy of degree distribution
    entropy_out(i) = nansum(-Prob.*log(Prob));
end
end
