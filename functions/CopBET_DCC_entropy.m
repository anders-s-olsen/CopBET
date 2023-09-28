% out = CopBET_DCC_entropy(in,varargin)
%
% Copenhagen Brain Entropy Toolbox: Dynamic conditional correlation entropy
% Evaluates DCC entropy as in Barrett et al., 2020. The DCC toolbox is run
% for each session. This produces a matrix of correlation coeficients for
% each volume. Then, for each roi-roi pair, the probability distribution 
% over the correlation coefficients across volumes is established, from 
% which the Shannon entropy is evaluated. As a warning, this function takes
% several days to run, which is the reason it also has an option to
% actually do the computation or to read them from a path
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells)g, e.g., different subjects or scan sessions.
%   compute: If not interested in running the very long computations, a
%   path can be specified to the location of saved DCC matrices
%   
% name-value pairs:
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%   compute: whether to actually compute these values or retrieve those
%   stored in memory from previous computations
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
%%
function out = CopBET_DCC_entropy(in,compute,varargin)

if nargin<2
    error('Please specify whether to actually run these computations')
end

[out,numworkers,in] = CopBET_function_init(in,varargin);

%load data
parfor(ses = 1:height(in),numworkers)
% for ses = 1:height(in)
    
    disp(['Working on DCC entropy calculations for session: ',num2str(ses)])
    ts = in{ses,1}{1};
    if compute~=true
        Ct2 = load([compute,'/DCC',...
        num2str(in.sesidx(ses)),'.mat']);
        Ct2 = Ct2.data;
    else
    Ct2 = DCC(ts-mean(ts));
    end
    [var{ses},entropy{ses}] = Barrett_analysis_no_correction(Ct2);
    
end
out.entropy = entropy';
out.var = var;
end

function [variance_out,entropy_out] = Barrett_analysis_no_correction(Ct2)



num_rois = size(Ct2,1);
variance_out = var(Ct2,[],3);
entropy_out = zeros(num_rois, num_rois);

for i = 1:num_rois
    for ii = i:num_rois
        edge_hist = histcounts(Ct2(i,ii,:),'Normalization','probability');
        entropy_out(i,ii) = nansum(-edge_hist.*log(edge_hist));
    end
end
entropy_out = entropy_out + entropy_out';


%%%%%%%% Checks
% sensible_data_check(variance_out);
entropy_check = entropy_out;
entropy_check(entropy_check==diag(diag(entropy_check)))=1;
sensible_data_check(entropy_check,'entropy matrix');

% check Ct2 is symmetric:
for i = 1:size(Ct2,3)
if norm(Ct2(:,:,i)-Ct2(:,:,i)')>1e-10
    error('Ct2 matrix not symmetric')
end
end

% check if Ct2 has values outside of correlation coefficient range
if any(Ct2(:)<-1)||any(Ct2(:)>1)
    warning('Not correlation coefficient range')
end

end



