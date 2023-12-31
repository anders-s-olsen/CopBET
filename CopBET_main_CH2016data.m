% This script serves as an introduction to the usage of the Copenhagen
% brain entropy toolbox available here: 
%
% https://github.com/anders-s-olsen/CopBET/
%
% Using this toolbox, 12 different brain entropy metrics may be evaluated.
% Some functions require a matrix (TxN) as input while others require the
% path to (denoised) 4D volumes in MNI-152 (2mm) space. Both options may be
% used for several subjects in a table, where the first column contains the
% matrices or paths. The functions have an option to retain the input table
% information and thus return output entropy values in a new column in the
% same table. The functions also have the option to be parallelized. Some
% functions require the input to be a path, including an atlas in the same
% space as the input data. 
%
% This MATLAB script evaluates 12 different entropy metrics for the openly 
% available Carhart-Harris 2016 acute IV LSD dataset available here: 
%
% https://openneuro.org/datasets/ds003059. 
%
% The user has the option to run the functions on the full dataset not
% included in the repository (instructions below), example data composed of
% two scans from the Carhart-Harris 2016 data set included in the respotory
% or the users own data. 
% This script is setup to run the functions using the example data. If the
% user wishes to run the functions on the full dataset, please download the
% data from OpenNeuro (link above) and put the data in CopBET/LSDdata
% (i.e., so that the folders sub-* are subfolders to CopBET/LSDdata). Then,
% please run the script LSDdata/LSDdata_ROI to extract the relevant
% atlases. Finally, the data-load-function
% CopBET_CarhartHarris_2016_data may take 'full' as input instead of
% 'example'. The example data set has been shortened to meet github
% requirements for file size. 
%
% Please run the script from the top CopBET directory for the
% 'CopBET_CarhartHarris_2016_data' function to work properly

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

clear
addpath(genpath(pwd))

[CopBETtbl,~,~] = CopBET_CarhartHarris_2016_data([],'ts','example');
%% LEiDA transition entropy (<1 minute on example data)
atlas = 'AAL90';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

% specify number of clusters as second argument
tbl = CopBET_LEiDA_transition_entropy(tbl,3,'keepdata',true,'parallel',true);
% outputs one entropy value pr scan in the table column tbl.entropy

plot_boxplots_CH2016(tbl.entropy,tbl,'LEiDA transition entropy')

% save in mastertable
CopBETtbl.LEiDA_transition_entropy = tbl.entropy;

%% Temporal Entropy (<5 minutes on example data)
atlas = 'SchaeferTian232';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

% Specify TR as 2nd argument
tbl = CopBET_temporal_entropy(tbl,2,'keepdata',true,'parallel',true);
% outputs one entropy value pr scan in the table column tbl.entropy

plot_boxplots_CH2016(tbl.entropy,tbl,'Temporal entropy')

CopBETtbl.temporal_entropy = tbl.entropy;

%% Metastate series complexity (<1 minute on example data)
atlas = 'Lausanne463';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

tbl = CopBET_metastate_series_complexity(tbl,'keepdata',true,'parallel',true);
% outputs one entropy value pr scan in the table column tbl.entropy

plot_boxplots_CH2016(tbl.entropy,tbl,'Metastate series complexity')

CopBETtbl.metastate_series_complexity = tbl.entropy;

%% von Neumann entropy (instant)
atlas = 'HarvardOxford_cort_subcort';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

tbl = CopBET_von_Neumann_entropy(tbl,'keepdata',true);
% outputs one entropy value pr scan in the table column tbl.entropy

plot_boxplots_CH2016(tbl.entropy,tbl,'von Neumann entropy')

CopBETtbl.von_Neumann_entropy = tbl.entropy;

%% Intranetwork synchrony (<1 minute on example data)
clearvars entropy

% This one should be given path to the denoised voxel-wise time series and
% an atlas. Here we extract 9 ROIs from smith20
smith_atlas = niftiread('Atlases/smith20_thresholded_2mm.nii');
smith_atlas = logical(smith_atlas(:,:,:,[1,2,5,6,9,11,12,14,15]+1));
[tbl,data,opts] = CopBET_CarhartHarris_2016_data([],'denoised_volumes','example');

tbl = CopBET_intranetwork_synchrony(tbl,smith_atlas,'keepdata',true,'parallel',true);
% outputs a vector of entropy values pr scan, each value corresponds to one
% network in the smith atlas

networks = {'motor','auditory','visual1','DMN','DAN','RFP','LFP','Salience','Visual2'};
% convert from table to matrix
for ses = 1:height(tbl)
    entropy(ses,:) = tbl.entropy{ses};
end

for net = 1:9
    CopBETtbl.(['Intranetwork_synchrony_',networks{net}]) = entropy(:,net);
    plot_boxplots_CH2016(entropy(:,net),tbl,['Intranetwork synchrony: ',networks{net}])
end

%% Dynamic conditional correlation (DCC) entropy (Several hours pr scan)
% This one doesn't work with the shortened scans. Please download the full
% dataset. 
clearvars entropy
atlas = 'Shen268';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','full');

% This one takes several days to run. Second argument is whether to actually run the script or to use saved previous outputs 
tbl = CopBET_DCC_entropy(tbl,true,'keepdata',true,'parallel',true);
% outputs a vector of entropy values pr scan, each value corresponds to one
% atlas edge/roi-to-roi connection (i.e., 268*267/2 unique values)

% unwrap roi-to-roi DCC edges to network-to-network
Shen268 = niftiread('Atlases/Shen268_2mm.nii');
Shen268labels = readtable('Atlases/shen_268_parcellation_networklabels_1.csv');

atlasnames = {'Medial_frontal','Frontoparietal','Deafult_mode','Subcortical_cerebellum',...
    'Motor','Visual1','Visual2','Visual_association'};

for network1 = 1:8
    for network2 = 1:8
        network_idx1 = find(Shen268labels.Network==network1);
        network_idx2 = find(Shen268labels.Network==network2);
        entropy = nan(height(tbl),1);
        for ses = 1:height(tbl)
            tmp = tbl.entropy{ses}(network_idx1,network_idx2);
            tmp = tmp(tmp~=0);
            entropy(ses) = mean(tmp);
        end
        CopBETtbl.(['DCC_entropy_',atlasnames{network1},'_',atlasnames{network2}]) = entropy;
        plot_boxplots_CH2016(entropy,tbl,['DCC entropy: ',atlasnames{network1},'_',atlasnames{network2}])
    end
end

%% Degree distribution entropy (instant)
clearvars entropy
atlas = 'HarvardOxford_cort_subcort';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

tbl = CopBET_degree_distribution_entropy(tbl,'keepdata',true,'parallel',true);
% outputs a vector of entropy values pr scan, each value corresponds to one
% integer "mean degree" at which the degree distribution entropy is evaluated

for degree = 1:100
    for ses = 1:height(tbl)
        entropy(ses,degree) = tbl.entropy{ses}(degree);
    end
    CopBETtbl.(['Degree_distribution_entropy_degree',num2str(degree)]) = entropy(:,degree);
end

% plot entropy for an example mean degree
degree_to_plot = 27;
plot_boxplots_CH2016(entropy(:,degree_to_plot),tbl,['Degree distribution entropy, degree ',num2str(degree_to_plot)])

%% Viol 2019 script (instant)
clearvars entropy
atlas = 'HarvardOxford_cort_subcort';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

tbl = CopBET_geodesic_entropy(tbl,'keepdata',true,'parallel',true);
% outputs a vector of entropy values pr scan, each value corresponds to one
% integer "mean degree" at which the degree distribution entropy is evaluated

for degree = 1:100
    for ses = 1:height(tbl)
        entropy(ses,degree) = tbl.entropy{ses}(degree);
    end
    CopBETtbl.(['Geodesic_entropy_degree',num2str(degree)]) = entropy(:,degree);
end

% plot entropy for an example mean degree
degree_to_plot = 27;
plot_boxplots_CH2016(entropy(:,degree_to_plot),tbl,['Geodesic entropy, degree ',num2str(degree_to_plot)])

%% Diversity coefficient (<1 minute on example data)
clearvars entropy
atlas = 'Craddock200';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');

tbl = CopBET_diversity_coefficient(tbl,'keepdata',true,'parallel',true);
% outputs a vector of entropy values pr scan, each value corresponds to one
% the diversity coefficient for a specific ROI in the (Craddock) atlas

for rois = 1:200
    for ses = 1:height(tbl)
        entropy(ses,rois) = tbl.entropy{ses}(rois);
    end
    CopBETtbl.(['Diversity_coefficient_roi',num2str(rois)]) = entropy(:,rois);
end

% plot entropy for an example ROI
ROI_to_plot = 111;
plot_boxplots_CH2016(entropy(:,ROI_to_plot),tbl,['Diversity Coefficient, ROI: ',num2str(ROI_to_plot)])

%% Sample entropy (several hours pr scan)
clearvars entropy

% This one requires an atlas to be specified. Here we use Yeo-17
atlas = niftiread('Atlases/Yeo17_liberal_2mm.nii');
[tbl,data,opts] = CopBET_CarhartHarris_2016_data([],'denoised_volumes','example');

tbl = CopBET_sample_entropy(tbl,atlas,true,'keepdata',true,'parallel',true);
% outputs a matrix of entropy values pr scan, each matrix contains values
% corresponding to a specific multi-scale sample entropy scale (1 to 5)
% along the rows, and a ROI in the (Yeo17) atlas along the columns. 

for scale = 1:5
    for session = 1:height(tbl)
        for ROI = 1:17
            entropy(session,scale,ROI) = tbl.entropy{session}(scale,ROI);
        end
    end
    for ROI = 1:17
        CopBETtbl.(['Sample_entropy_scale',num2str(scale),'_ROI',num2str(ROI)]) = entropy(:,scale,ROI);
    end
end
scale_to_plot = 2;
ROI_to_plot = 7;
plot_boxplots_CH2016(entropy(:,scale_to_plot,ROI_to_plot),tbl,['Sample entropy, scale: ',num2str(scale_to_plot),', ROI: ',num2str(ROI_to_plot)])

%% Motif connectivity entropy (<1 minute on example data)
clearvars entropy

% This one requires an atlas. Here we specify 4 ROIs in one atlas:
ROI1 = logical(niftiread('Atlases/r2Tag_-34_-22_-16_LHip_5mm_sphere.nii'));
ROI2 = logical(niftiread('Atlases/r2Tag_26_-22_-16_RHip_5mm_sphere.nii'));
ROI3 = logical(niftiread('Atlases/r2Tag_-2_22_28_LACC_5mm_sphere.nii'));
ROI4 = logical(niftiread('Atlases/r2Tag_4_34_18_RACC_5mm_sphere.nii'));
atlas = ROI1+2*ROI2+3*ROI3+4*ROI4;

[tbl,data,opts] = CopBET_CarhartHarris_2016_data([],'denoised_volumes','example'); 

% This one also needs motion correction time series in the second column of
% 'tbl'. For the CH2016, we don't have these available. For the sake of
% example, here we make some random data:
for i = 1:height(tbl)
    tbl.rp{i} = randn(217,6);
end

tbl = CopBET_motif_connectivity_entropy(tbl,atlas,2,'keepdata',true,'parallel',true);
% outputs a vector of entropy values, where each value corresponds to a
% window-length used to evaluate a windowed connectivity matrix inside the
% function.

window_lengths = 15:150;
for wl = 1:136
    for session = 1:height(tbl)
        entropy(session,wl) = tbl.entropy{session}(wl);
    end
    CopBETtbl.(['Motif_conenctivity_entropy_windowsize',num2str(window_lengths(wl))]) = entropy(:,wl);
end

% plot entropy for an example window length
wl_to_plot = 4;
plot_boxplots_CH2016(entropy(:,wl_to_plot),tbl,['Motif connectivity entropy, window size: ',num2str(wl_to_plot)])

%% Varley script, temporal LZ78 (<5 minutes on example data)
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');
tbl = CopBET_time_series_complexity(tbl,'LZ78temporal',true,true);

CopBETtbl.Time_series_complexity_temporal = tbl.entropy;
% outputs one entropy value pr scan in the table column tbl.entropy

plot_boxplots_CH2016(tbl.entropy,tbl,'Time series complexity, temporal')

%% Varley script, spatial LZ78 (<5 minutes on example data)
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts','example');
tbl = CopBET_time_series_complexity(tbl,'LZ78spatial',true,true);

CopBETtbl.Time_series_complexity_spatial = tbl.entropy;
% outputs one entropy value pr scan in the table column tbl.entropy

plot_boxplots_CH2016(tbl.entropy,tbl,'Time series complexity, spatial')

