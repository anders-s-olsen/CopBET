% out = CopBET_sample_entropy(in,atlas,compute,varargin)
%
% Copenhagen Brain Entropy Toolbox: Sample entropy
% Evaluates sample entropy as in Lebedev et al., 2016. Sample entropy is
% very hard to explain, so perhaps it is best to read the paper
%
% Input:
%   in: char with the path to the input, denoised voxel-wise time series or
%   a table where the first column contains
%   chars (in cells), e.g., different subjects or scan sessions.
%   atlas: atlas matrix of values
%   compute: something, currently required to be true
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

function out = CopBET_sample_entropy(in,atlas,compute,varargin)

if nargin<3
    error('please specify atlas and whether to do computations or not')
end
if compute~=true
    error('this option is not currently implemented, please set compute to true')
end

[out,numworkers,in,NRUspecific] = CopBET_function_init(in,varargin);

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
    if NRUspecific
        if compute~=true
            for a = 1:numel(scale)
                V = niftiread([compute,'/',num2str(in.sesidx(ses)),'_scale',num2str(a)]);
                V = double(V(:));
                
                for roi = 1:num_rois
                    tmp = V(atlas(:)==roi);
                    entropy{ses}(a,roi) = mean(tmp(tmp~=0));
                end
            end
        end
        continue
    end
    path = in{ses,1}{1};
    image_4D = double(niftiread(path)); %4D series
    if NRUspecific
        if ~isempty(regexp(path,'denoisedn'))
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
        MSE = nan(length(brainVox),1);
        parfor (vox = 1:length(brainVox),numworkers)
            %         for vox = 1:length(brainVox)
            [row,col,sl] = ind2sub(imgSize,brainVox(vox));
            ts = squeeze(image_4D(row,col,sl,:));
            
            r_val = r*std(double(ts));
            tmp = sample_entropy(m,r_val,ts,scale(a));
            MSE(vox) = tmp(1);
            
        end
        MSE2 = nan(imgSize);
        for vox = 1:numel(brainVox)
            [row,col,sl] = ind2sub(imgSize,brainVox(vox));
            MSE2(row,col,sl) = MSE(vox);
        end
        
        for roi = 1:num_rois
            tmp = MSE2(atlas==roi);
            entropy{ses}(a,roi) = mean(tmp(tmp~=0));
        end
        
    end
    
    
    
end
out.entropy = entropy';
end

