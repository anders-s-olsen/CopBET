% out = CopBET_metastate_series_complexity(in,keepdata,parallel)
%
% Copenhagen Brain Entropy Toolbox: Metastate series complexity.
% Evaluates Lempel-Ziv complexity (the LZ76 exhaustive algorithm) for a 
% binary series of metastates. The metastates are created by running
% K-means with K=4 states on the fMRI data using the correlation distance
% measure. Subsequently, since the four state centroids are typically
% diametrical opposites of one of the others, the four states are grouped
% to two 'metastates'. The activation sequence of these metastates then
% forms the input to the LZ76 algorithm.
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before clustering, e.g.,
%   different subjects or scan sessions.
% name-value pairs:
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%   parallel: Whether to run 200 replicates of k-means in parallel [true]
%
%
% Neurobiology Research Unit, 2023
% Please cite McCulloch, Olsen et al., 2023: "Navigating Chaos in
% Psychedelic Neuroimaging: A Rigorous Empirical Evaluation of the Entropic
% Brain Hypothesis" if you use CopBET in your studies. Please read the
% paper to get a notion of our recommendations regarding the use of the
% specific methodologies in the toolbox.

% ASO 9/3-2023

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

function out = CopBET_metastate_series_complexity(in,varargin)

[out,numworkers,in] = CopBET_function_init(in,varargin);

% Do LEiDA and concatenate data for clustering
disp('Concatenating data')
datasizes = nan(height(in),1);
for ses = 1:height(in)
    datasizes(ses) = size(in{ses,1}{1},1);
end
data_all = nan(sum(datasizes),size(in{1,1}{1},2));
c = 1;
for ses = 1:height(in)
    tmp = in{ses,1}{1};
    tmp = tmp-mean(tmp);
    data_all(c:c+datasizes(ses)-1,:) = tmp;
    c = c+datasizes(ses);
end

if any(isnan(data_all(:)))
    error('Error, not all atlas regions present in the data or some input data are nan')
end

disp('running kmeans and LZ calculations')
entropy = run_singleton_clustering(datasizes,data_all,numworkers);


out.entropy = entropy;



end

%% functions
function entropy = run_singleton_clustering(datasizes,ts_all,numworkers)
nreps = 200; % how many times to repeat clustering. will choose lowest error solution
distanceMethod = 'correlation';
maxI = 1000; % how many times you allow kmeans to try to converge

numClusters = 4;
if numworkers>0
[partition,centroids] = kmeans(ts_all,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI,...
    'Display','final','Options',statset('UseParallel',1));
else
    [partition,centroids] = kmeans(ts_all,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI,...
    'Display','final');
end

% group states
possible_states = 1:numClusters;
cencorr = corr(centroids');
for k = 1:numClusters
    [~,idx(k)] = min(cencorr(:,k));
end
for k = 1:numClusters
    if idx(idx(k))~=k
        error('Wrong metastate grouping')
    end
end

partition1 = [1,idx(1)];
partition2 = setdiff(possible_states,partition1);

newPartition = false(size(partition));
newPartition(ismember(partition,partition2)) = true;

tpts_traversed = 0;
entropy = nan(length(datasizes),1);
for h = 1:length(datasizes)
    
    T = tpts_traversed+1:tpts_traversed+datasizes(h);
    entropy(h) = calc_lz_complexity(newPartition(T),'exhaustive',1);
    
    tpts_traversed = tpts_traversed + datasizes(h);
%     disp(['Calculating LZ for #',num2str(h),' of ',num2str(height(tbl))])
end


end
