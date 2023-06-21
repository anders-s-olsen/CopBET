%
% This MATLAB script evaluates 12 different entropy metrics for the openly 
% available Carhart-Harris 2016 acute IV LSD dataset available here: 
%
% https://openneuro.org/datasets/ds003059. 
%
% Please run the script LSDdata/LSDdata_ROI to extract the relevant atlases
%
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
% Please run the script from the top CopBET directory for the
% 'CopBET_CarhartHarris_2016_data' function to properly load data

clear
addpath(genpath('/mrdata/np2/p3/entropy/CopBET'))
LSDsessions = [1,3];

[CopBETtbl,~,~] = CopBET_CarhartHarris_2016_data([],'ts');
%% LEiDA transition entropy

atlas = 'aal90';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');

% specify number of clusters
tbl = CopBET_LEiDA_transition_entropy(tbl,3,'keepdata',true,'parallel',true);

%%% Figure
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Kringelbach20, p=',num2str(p)]})
end

CopBETtbl.LEiDA_transition_entropy = tbl.entropy;

%% Temporal Entropy
atlas = 'SchaeferTian232';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');

% Specify TR as 2nd argument
tbl = CopBET_temporal_entropy(tbl,2,'keepdata',true,'parallel',true);

%%% Figure
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Luppi21, p=',num2str(p)]})
end

CopBETtbl.temporal_entropy = tbl.entropy;

%% metastate series complexity
atlas = 'Lausanne463';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');
tbl = CopBET_metastate_series_complexity(tbl,'keepdata',true,'parallel',true);

%%% Figure
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Singleton22, p=',num2str(p)]})
end

CopBETtbl.metastate_series_complexity = tbl.entropy;

%% von Neumann entropy
atlas = 'HarvardOxford105';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');
tbl = CopBET_von_Neumann_entropy(tbl,'keepdata',true);

%%% Figure
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Felippe21, p=',num2str(p)]})
end

CopBETtbl.von_Neumann_entropy = tbl.entropy;

%% Intranetwork synchrony
clearvars entropy

% This one should be given path to the denoised voxel-wise time series and
% an atlas. Here we extract 9 ROIs from smith20
smith_atlas = niftiread('Atlases/smith20_thresholded_2mm.nii');
smith_atlas = logical(smith_atlas(:,:,:,[1,2,5,6,9,11,12,14,15]+1));

[tbl,data,opts] = CopBET_CarhartHarris_2016_data([],'denoised_volumes');

tbl = CopBET_intranetwork_synchrony(tbl,smith_atlas,'keepdata',true,'parallel',true);
networks = {'motor','auditory','visual1','DMN','DAN','RFP','LFP','Salience','Visual2'};

% convert from table to matrix
for ses = 1:height(tbl)
    entropy(ses,:) = tbl.entropy{ses};
end

for net = 1:9
    CopBETtbl.(['Intranetwork_synchrony_',networks{net}]) = entropy(:,net);
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(entropy(tbl.session==LSDsessions(LSDses),net),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),net),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),net)]','k-')
        [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),net),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),net));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['CH14, ',networks{net},', p=',num2str(p)]})
    end
    
end

%% Dynamic conditional correlation (DCC) entropy
clearvars entropy
atlas = 'Shen268';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');

% This one takes several days to run. 
tbl = CopBET_DCC_entropy(tbl,true,'keepdata',true,'parallel',true);

% unwrap roi-to-roi DCC edges to network-to-network
Shen218 = niftiread('/mrdata/np2/p3/entropy/ROIs/shen_2mm_268_parcellation_wo_cerebellum.nii');
Shen268labels = readtable('/mrdata/np2/p3/entropy/ROIs/shen_268_parcellation_networklabels_1.csv');

atlasnames = {'Medial_frontal','Frontoparietal','Deafult_mode','Subcortical_cerebellum',...
    'Motor','Visual1','Visual2','Visual_association'};

for network1 = 1:8
    for network2 = 1:8
        network_idx1 = find(Shen218labels==network1);
        network_idx2 = find(Shen218labels==network2);
        entropy = nan(height(tbl),1);
        for ses = 1:height(tbl)
            tmp = tbl.entropy{ses}(network_idx1,network_idx2);
            tmp = tmp(tmp~=0);
            entropy(ses) = mean(tmp);
            
        end
        CopBETtbl.(['DCC_entropy_',atlasnames{network1},'_',atlasnames{network2}]) = entropy;
        
        figure,hold on,
        for LSDses = 1:2
            subplot(1,2,LSDses)
            boxplot(entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
            plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
            [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
            title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Barrett20, ',atlasnames{network1},'_',atlasnames{network2},', p=',num2str(p)]})
        end
        
    end
end

%% Degree distribution entropy
clearvars entropy
atlas = 'HarvardOxford105';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');
tbl = CopBET_degree_distribution_entropy(tbl,'keepdata',true,'parallel',true);

for degree = 1:100
    for ses = 1:height(tbl)
        entropy(ses,degree) = tbl.entropy{ses}(degree);
    end
    CopBETtbl.(['Degree_distribution_entropy_degree',num2str(degree)]) = entropy(:,degree);
end

degree_toplot = 27;
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(entropy(tbl.session==LSDsessions(LSDses),degree_toplot),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree_toplot),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree_toplot)]','k-')
    [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree_toplot),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree_toplot));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Viol17, p=',num2str(p)]})
end

%% Viol 2019 script
clearvars entropy
atlas = 'HarvardOxford105';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');
tbl = CopBET_geodesic_entropy(tbl,'keepdata',true,'parallel',true);

for degree = 1:100
    for ses = 1:height(tbl)
        entropy(ses,degree) = tbl.entropy{ses}(degree);
    end
    CopBETtbl.(['Geodesic_entropy_degree',num2str(degree)]) = entropy(:,degree);
end

degree_toplot = 27;
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(entropy(tbl.session==LSDsessions(LSDses),degree_toplot),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree_toplot),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree_toplot)]','k-')
    [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree_toplot),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree_toplot));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Viol19, p=',num2str(p)]})
end

%% Diversity coefficient
clearvars entropy
atlas = 'Craddock200';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');

tbl = CopBET_diversity_coefficient(tbl,'keepdata',true,'parallel',true);


for rois = 1:200
    for ses = 1:height(tbl)
        entropy(ses,rois) = tbl.entropy{ses}(rois);
    end
    CopBETtbl.(['Lebedev15_mrspecific_roi',num2str(rois)]) = entropy(:,rois);

    % you can choose to plot all 200 rois...
%     figure,hold on,
%     for LSDses = 1:2
%         subplot(1,2,LSDses)
%         boxplot(entropy(tbl.session==LSDsessions(LSDses),rois),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
%         plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),rois),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),rois)]','k-')
%         [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),rois),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),rois));
%         title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Lebedev15, p=',num2str(p)]})
%     end

end
%% Sample entropy
clearvars entropy

% This one requires an atlas to be specified. Here we use Yeo-17
atlas = niftiread('Atlases/Yeo17_liberal_2mm.nii');
[tbl,data,opts] = CopBET_CarhartHarris_2016_data([],'denoised_volumes');
tbl = CopBET_sample_entropy(tbl,atlas,true,'keepdata',true,'parallel',true);


for scale = 1:5
    for session = 1:height(tbl)
        for ROI = 1:17
            entropy(session,scale,ROI) = tbl.entropy{session}(scale,ROI);
        end
    end
    for ROI = 1:17
        CopBETtbl.(['Lebedev16_scale',num2str(scale),'_ROI',num2str(ROI)]) = entropy(:,scale,ROI);
        
        if strcmp(dataset,'LSD')
            figure,hold on,
            for LSDses = 1:2
                subplot(1,2,LSDses)
                boxplot(entropy(tbl.session==LSDsessions(LSDses),scale,ROI),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
                plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),scale,ROI),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),scale,ROI)]','k-')
                [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),scale,ROI),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)),scale,ROI);
                title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Lebedev16, p=',num2str(p)]})
            end
        end
        
    end
end

%% Tagliazucchi script
clearvars entropy

% This one requires an atlas. Here we specify 4 ROIs in one atlas:
ROI1 = logical(niftiread('Atlases/r2Tag_-34_-22_-16_LHip_5mm_sphere.nii'));
ROI2 = logical(niftiread('Atlases/r2Tag_26_-22_-16_RHip_5mm_sphere.nii'));
ROI3 = logical(niftiread('Atlases/r2Tag_-2_22_28_LACC_5mm_sphere.nii'));
ROI4 = logical(niftiread('Atlases/r2Tag_4_34_18_RACC_5mm_sphere.nii'));
atlas = ROI1+2*ROI2+3*ROI3+4*ROI4;
[tbl,data,opts] = CopBET_CarhartHarris_2016_data([],'denoised_volumes');
tbl = CopBET_motif_connectivity_entropy(tbl,atlas,2,'keepdata',true,'parallel',true);

window_lengths = 15:150;
for wl = 1:136
    for session = 1:height(tbl)
        entropy(session,wl) = tbl.entropy{session}(wl);
    end
    CopBETtbl.(['Tagliazucchi14_windowsize',num2str(window_lengths(wl))]) = entropy(:,wl);
    
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(entropy(tbl.session==LSDsessions(LSDses),wl),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),wl),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),wl)]','k-')
        [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),wl),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),wl));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Tagliazucchi14, p=',num2str(p)]})
    end
    
end

%% Varley script, temporal LZ78
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');
tbl = CopBET_time_series_complexity(tbl,'LZ78temporal',true,true);

CopBETtbl.LZ78temporal = tbl.entropy;

figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Varley20, LZ78temporal, p=',num2str(p)]})
end

%% Varley script, spatial LZ78
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,'ts');
tbl = CopBET_time_series_complexity(tbl,'LZ78spatial',true,true);

CopBETtbl.LZ78spatial = tbl.entropy;

figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Varley20, LZ78spatial, p=',num2str(p)]})
end

