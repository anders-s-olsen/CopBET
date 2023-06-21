%%
clear
addpath(genpath('/mrdata/np2/p3/entropy/CopBET'))
dataset = 'NRU';
LSDsessions = [1,3];

load('/mrdata/np2/p3/entropy/code/Analyses/mastertbl6.mat')

[CopBETtbl,~,~] = CopBET_load_data(dataset,[],'ts');
% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Kringelbach script
% LEiDA_transition_entropy
atlas = 'aal90';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
if strcmp(dataset,'NRU')
    mr45idx = find(strcmp(tbl.mr,'mr45'));
    mr001idx = find(strcmp(tbl.mr,'mr001'));
    tbl45 = CopBET_LEiDA_transition_entropy(tbl(mr45idx,:),2,3,true,true);
    tbl001 = CopBET_LEiDA_transition_entropy(tbl(mr001idx,:),2);
    tbl = vertcat(tbl45,tbl001);
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),tbl.entropy(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),tbl.entropy(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Kringelbach20_mrspecific(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Kringelbach20_mrspecific(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
elseif strcmp(dataset,'LSD')
    tbl = CopBET_LEiDA_transition_entropy(tbl,2,3,true,true);
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Kringelbach20, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Kringelbach2020'],'-dpng','-r300')
end

CopBETtbl.Kringelbach20_mrspecific = tbl.entropy;



%% Luppi script
% Temporal_entropy, takes forever because of 500 kmeans replicates PER
% SUBJECT
atlas = 'SchaeferTian232';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_temporal_entropy(tbl,2,true,true);

CopBETtbl.Luppi21 = tbl.entropy;

if strcmp(dataset,'LSD')
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Luppi21, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Luppi2021'],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),tbl.entropy(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),tbl.entropy(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Luppi21(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Luppi21(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
end
%% Singleton script
% metastate series complexity
atlas = 'Lausanne462';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
if strcmp(dataset,'NRU')
    mr45idx = find(strcmp(tbl.mr,'mr45'));
    mr001idx = find(strcmp(tbl.mr,'mr001'));
    tbl45 = CopBET_metastate_series_complexity(tbl(mr45idx,:),true,true);
    tbl001 = CopBET_metastate_series_complexity(tbl(mr001idx,:),true,true);
    tbl = vertcat(tbl45,tbl001);
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),tbl.entropy(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),tbl.entropy(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Singleton21_mrspecific(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Singleton21_mrspecific(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
elseif strcmp(dataset,'LSD')
    tbl = CopBET_metastate_series_complexity(tbl,true,true);
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Singleton22, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Singleton2022'],'-dpng','-r300')
end

CopBETtbl.Singleton21_mrspecific = tbl.entropy;


%% Felippe script
% von Neumann entropy
atlas = 'HarvardOxford105';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_von_Neumann_entropy(tbl,true);

CopBETtbl.Felippe22 = tbl.entropy;
if strcmp(dataset,'LSD')
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Felippe21, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Felippe2021'],'-dpng','-r300')
elseif strcmp(dataset, 'NRU')
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),tbl.entropy(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),tbl.entropy(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Felippe22(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Felippe22(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
end
%% Carhart Harris 2014 script (only works on bohr)
clearvars entropy
% intranetwork_synchrony, smith20
% This one should be given path to the denoised voxel-wise time series and
% an atlas
% specific for this function: load smith atlas
smith_atlas = niftiread('/mrdata/np2/p3/entropy/ROIs/smith20_thresholded_2mm.nii');
smith_atlas = logical(smith_atlas(:,:,:,[1,2,5,6,9,11,12,14,15]+1));
[tbl,data,opts] = CopBET_load_data(dataset,[],'denoised_volumes');
tbl = CopBET_intranetwork_synchrony(tbl,smith_atlas,true,true,true);
networks = {'motor','auditory','visual1','DMN','DAN','RFP','LFP','Salience','Visual2'};

for ses = 1:height(tbl)
    entropy(ses,:) = tbl.entropy{ses};
end

if strcmp(dataset,'LSD')
    for net = 1:9
        figure,hold on,
        for LSDses = 1:2
            subplot(1,2,LSDses)
            boxplot(entropy(tbl.session==LSDsessions(LSDses),net),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
            plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),net),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),net)]','k-')
            [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),net),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),net));
            title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['CH14, ',networks{net},', p=',num2str(p)]})
        end
        print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/CH2014',networks{net}],'-dpng','-r300')
    end
elseif strcmp(dataset,'NRU')
    
    for net = 1:9
        CopBETtbl.(['CarhartHarris14_',networks{net}]) = entropy(:,net);
        figure,
        subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),entropy(strcmp(tbl.mr,'mr45'),net),'k.'),lsline
        subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),entropy(strcmp(tbl.mr,'mr001'),net),'k.'),lsline
        subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['CarhartHarris14_',networks{net}])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
        subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['CarhartHarris14_',networks{net}])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    end
    
    save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
end
%% Barrett script
clearvars entropy
atlas = 'Shen218';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
% tbl = CopBET_DCC_entropy(tbl,true,false);
tbl = CopBET_DCC_entropy(tbl,true,true);

Shen218 = niftiread('/mrdata/np2/p3/entropy/ROIs/shen_2mm_268_parcellation_wo_cerebellum.nii');
Shen268labels = readtable('/mrdata/np2/p3/entropy/ROIs/shen_268_parcellation_networklabels_1.csv');
un = unique(Shen218);un(1)=[];
Shen218labels = Shen268labels.Network(un); %remove cerebellar regions
atlasnames = {'Medial_frontal','Frontoparietal','Deafult_mode','Subcortical_cerebellum',...
    'Motor','Visual1','Visual2','Visual_association'};
for network1 = 1:8
    for network2 = 1:8
        network_idx1 = find(Shen218labels==network1);
        network_idx2 = find(Shen218labels==network2);
        entropy = nan(height(tbl),1);
        for ses = 1:height(tbl)
            tmp = tbl.entropy{ses}(network_idx1,network_idx2);
            %             tmp2 = tbl.entropy{ses}(network_idx2,network_idx1);
            
            tmp = tmp(tmp~=0);
            %             tmp2 = tmp2(tmp2~=0);
            %             entropy(ses) = mean([tmp;tmp2]);
            entropy(ses) = mean(tmp);
            
        end
        CopBETtbl.(['Barrett20_',atlasnames{network1},'_',atlasnames{network2}]) = entropy;
        if strcmp(dataset,'NRU')
%             tbl.(['Barrett20_',atlasnames{network1},'_',atlasnames{network2}]) = entropy;
%             figure,
%             subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),entropy(strcmp(tbl.mr,'mr45')),'k.'),lsline
%             subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),entropy(strcmp(tbl.mr,'mr001')),'k.'),lsline
%             subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['Barrett20_',atlasnames{network1},'_',atlasnames{network2}])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
%             subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['Barrett20_',atlasnames{network1},'_',atlasnames{network2}])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
        elseif strcmp(dataset,'LSD')
            figure,hold on,
            for LSDses = 1:2
                subplot(1,2,LSDses)
                boxplot(entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
                plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
                [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
                title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Barrett20, ',atlasnames{network1},'_',atlasnames{network2},', p=',num2str(p)]})
            end
            print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Barrett20_',atlasnames{network1},'_',atlasnames{network2}],'-dpng','-r300')
        end
        
    end
end



% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Viol 2017 script
clearvars entropy
atlas = 'HarvardOxford105';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_degree_distribution_entropy(tbl,true,true);

for degree = 1:100
    for ses = 1:height(tbl)
        entropy(ses,degree) = tbl.entropy{ses}(degree);
    end
    CopBETtbl.(['Viol17_degree',num2str(degree)]) = entropy(:,degree);
end

if strcmp(dataset,'LSD')
    degree = 27;
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(entropy(tbl.session==LSDsessions(LSDses),degree),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree)]','k-')
        [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Viol17, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Viol2017_',num2str(degree)],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
    
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),ex(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),ex(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['Viol17_degree',num2str(d)])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['Viol17_degree',num2str(d)])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
end
%% Viol 2019 script
clearvars entropy
atlas = 'HarvardOxford105';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_geodesic_entropy(tbl,true,true);

for degree = 1:100
    for ses = 1:height(tbl)
        entropy(ses,degree) = tbl.entropy{ses}(degree);
    end
    CopBETtbl.(['Viol19_degree',num2str(degree)]) = entropy(:,degree);
end


if strcmp(dataset,'LSD')
    degree = 27;
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(entropy(tbl.session==LSDsessions(LSDses),degree),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree)]','k-')
        [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),degree),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),degree));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Viol19, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Viol2019_',num2str(degree)],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
    
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),ex(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),ex(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['Viol19_degree',num2str(d)])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['Viol19_degree',num2str(d)])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
end
%% Lebedev 2015 script
clearvars entropy
atlas = 'Craddock181';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
if strcmp(dataset,'NRU')
mr45idx = find(strcmp(tbl.mr,'mr45'));
mr001idx = find(strcmp(tbl.mr,'mr001'));
tbl45 = CopBET_diversity_coefficient(tbl(mr45idx,:),true,true);
tbl001 = CopBET_diversity_coefficient(tbl(mr001idx,:),true,true);
tbl = vertcat(tbl45,tbl001);
elseif strcmp(dataset,'LSD')
    tbl = CopBET_diversity_coefficient(tbl,true,true);
end

for rois = 1:181
    for ses = 1:height(tbl)
        entropy(ses,rois) = tbl.entropy{ses}(rois);
        
    end
    CopBETtbl.(['Lebedev15_mrspecific_roi',num2str(rois)]) = entropy(:,rois);

if strcmp(dataset,'LSD')
    degree = 27;
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(entropy(tbl.session==LSDsessions(LSDses),rois),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),rois),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),rois)]','k-')
        [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),rois),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)),rois);
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Lebedev15, p=',num2str(p)]})
        print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Lebedev2015_',num2str(rois)],'-dpng','-r300')
    end
elseif strcmp(dataset,'NRU')
    
    % figure,
    % subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),ex(strcmp(tbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),ex(strcmp(tbl.mr,'mr001')),'k.'),lsline
    % subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['Viol19_degree',num2str(d)])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
    % subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['Viol19_degree',num2str(d)])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
    % save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
end
end
%% Lebedev 2016 script
clearvars entropy
atlas = niftiread('/mrdata/np2/p3/entropy/ROIs/Yeo2011_17Networks_MNI152_xmm_liberal_atlas.nii');
[tbl,data,opts] = CopBET_load_data(dataset,[],'denoised_volumes');
tbl = CopBET_sample_entropy(tbl,atlas,true,true,true,true);


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
                print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Lebedev2016_',num2str(scale),'_',num2str(ROI)],'-dpng','-r300')
            end
        end
        
    end
end

% a = 1;
% roi = 1;
% for ses = 1:height(tbl)
%     ex(ses) = tbl.entropy{ses}(a,roi);
% end

% figure,
% subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),ex(strcmp(tbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),ex(strcmp(tbl.mr,'mr001')),'k.'),lsline
% subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['Lebedev16_scale',num2str(a),'_ROI',num2str(roi)])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['Lebedev16_scale',num2str(a),'_ROI',num2str(roi)])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Tagliazucchi script
clearvars entropy
ROI1 = logical(niftiread('/mrdata/np2/p3/entropy/ROIs/r2Tag_-34_-22_-16_LHip_5mm_sphere.nii'));
ROI2 = logical(niftiread('/mrdata/np2/p3/entropy/ROIs/r2Tag_26_-22_-16_RHip_5mm_sphere.nii'));
ROI3 = logical(niftiread('/mrdata/np2/p3/entropy/ROIs/r2Tag_-2_22_28_LACC_5mm_sphere.nii'));
ROI4 = logical(niftiread('/mrdata/np2/p3/entropy/ROIs/r2Tag_4_34_18_RACC_5mm_sphere.nii'));
atlas = ROI1+2*ROI2+3*ROI3+4*ROI4;
[tbl,data,opts] = CopBET_load_data(dataset,[],'denoised_volumes');
tbl = CopBET_motif_connectivity_entropy(tbl,[],atlas,2,true,true,true);

window_lengths = 15:150;
for wl = 1:136
    for session = 1:height(tbl)
        entropy(session,wl) = tbl.entropy{session}(wl);
    end
    CopBETtbl.(['Tagliazucchi14_windowsize',num2str(window_lengths(wl))]) = entropy(:,wl);
    
    if strcmp(dataset,'LSD')
        figure,hold on,
        for LSDses = 1:2
            subplot(1,2,LSDses)
            boxplot(entropy(tbl.session==LSDsessions(LSDses),wl),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
            plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),wl),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),wl)]','k-')
            [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses),wl),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses),wl));
            title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Tagliazucchi14, p=',num2str(p)]})
            print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Tagliazucchi2014_',num2str(scale),'_',num2str(ROI)],'-dpng','-r300')
        end
    elseif strcmp(dataset,'NRU')
        
        if ismember(wl,[1:10:136])
            figure,
            subplot(2,2,1),plot(tbl.PPL(strcmp(tbl.mr,'mr45')),CopBETtbl.(['Tagliazucchi14_windowsize',num2str(window_lengths(wl))])(strcmp(tbl.mr,'mr45')),'k.'),lsline
            subplot(2,2,2),plot(tbl.PPL(strcmp(tbl.mr,'mr001')),CopBETtbl.(['Tagliazucchi14_windowsize',num2str(window_lengths(wl))])(strcmp(tbl.mr,'mr001')),'k.'),lsline
            subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.(['Tagliazucchi14_windowsize',num2str(window_lengths(wl))])(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
            subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.(['Tagliazucchi14_windowsize',num2str(window_lengths(wl))])(strcmp(mastertbl.mr,'mr001')),'k.'),lsline
        end
    end
end
save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Varley script, 
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_time_series_complexity(tbl,'LZ76temporal',true,true);

CopBETtbl.LZ76temporal = tbl.entropy;

if strcmp(dataset,'LSD')
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Varley20, LZ76temporal, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Varley2020_LZ76temporal'],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
%     figure,
% subplot(2,2,1),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr45')),CopBETtbl.LZ76temporal(strcmp(CopBETtbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,2),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr001')),CopBETtbl.LZ76temporal(strcmp(CopBETtbl.mr,'mr001')),'k.'),lsline
% subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr001')),'k.'),lsline

end


% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Varley script, temporal LZ78
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_time_series_complexity(tbl,'LZ76spatial',true,true);

CopBETtbl.LZ76spatial = tbl.entropy;

if strcmp(dataset,'LSD')
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Varley20, LZ76spatial, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Varley2020_LZ76spatial'],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
%     figure,
% subplot(2,2,1),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr45')),CopBETtbl.LZ76spatial(strcmp(CopBETtbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,2),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr001')),CopBETtbl.LZ76spatial(strcmp(CopBETtbl.mr,'mr001')),'k.'),lsline
% subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr001')),'k.'),lsline

end
% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Varley script, temporal LZ78
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_time_series_complexity(tbl,'LZ78temporal',true,true);

CopBETtbl.LZ78temporal = tbl.entropy;

if strcmp(dataset,'LSD')
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Varley20, LZ78temporal, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Varley2020_LZ78temporal'],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
%     figure,
% subplot(2,2,1),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr45')),CopBETtbl.LZ78temporal(strcmp(CopBETtbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,2),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr001')),CopBETtbl.LZ78temporal(strcmp(CopBETtbl.mr,'mr001')),'k.'),lsline
% subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr001')),'k.'),lsline

end
% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')
%% Varley script, spatial LZ78
atlas = 'Schaefer1000';
[tbl,data,opts] = CopBET_load_data(dataset,atlas,'ts');
tbl = CopBET_time_series_complexity(tbl,'LZ78spatial',true,true);

CopBETtbl.LZ78spatial = tbl.entropy;

if strcmp(dataset,'LSD')
    figure,hold on,
    for LSDses = 1:2
        subplot(1,2,LSDses)
        boxplot(tbl.entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
        plot([ones(15,1),2*ones(15,1)]',[tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
        [h,p,ci,stats] = ttest(tbl.entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),tbl.entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
        title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],['Varley20, LZ78spatial, p=',num2str(p)]})
    end
    print(['/mrdata/np2/p3/entropy/CopBET/LSDdataresults/Varley2020_LZ78spatial'],'-dpng','-r300')
elseif strcmp(dataset,'NRU')
%     figure,
% subplot(2,2,1),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr45')),CopBETtbl.LZ78spatial(strcmp(CopBETtbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,2),plot(CopBETtbl.PPL(strcmp(CopBETtbl.mr,'mr001')),CopBETtbl.LZ78spatial(strcmp(CopBETtbl.mr,'mr001')),'k.'),lsline
% subplot(2,2,3),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr45')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr45')),'k.'),lsline
% subplot(2,2,4),plot(mastertbl.PPL(strcmp(mastertbl.mr,'mr001')),mastertbl.Varley20(strcmp(mastertbl.mr,'mr001')),'k.'),lsline

end
% save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl','CopBETtbl')

%% fix table
load('/mrdata/np2/p3/entropy/CopBET/CopBETtbl')
markforremove = false(width(CopBETtbl),1);

for col = 1:width(CopBETtbl)
    
    if ~isempty(regexp(CopBETtbl.Properties.VariableNames{col},'Viol17_degree'))
        if length(CopBETtbl.Properties.VariableNames{col})==14
            CopBETtbl.Properties.VariableNames{col} = ['Viol17_degree','0',CopBETtbl.Properties.VariableNames{col}(end)];
        end
        if str2double(CopBETtbl.Properties.VariableNames{col}(end-1:end))>48
            markforremove(col)=true;
        end
    elseif ~isempty(regexp(CopBETtbl.Properties.VariableNames{col},'Viol19_degree'))
        if length(CopBETtbl.Properties.VariableNames{col})==14
            CopBETtbl.Properties.VariableNames{col} = ['Viol19_degree','0',CopBETtbl.Properties.VariableNames{col}(end)];
        end
        if str2double(CopBETtbl.Properties.VariableNames{col}(end-1:end))>53
            markforremove(col)=true;
        end
    elseif ~isempty(regexp(CopBETtbl.Properties.VariableNames{col},'Lebedev15_mrspecific_roi'))
        if length(CopBETtbl.Properties.VariableNames{col})==25
            CopBETtbl.Properties.VariableNames{col} = ['Lebedev15_mrspecific_roi','00',CopBETtbl.Properties.VariableNames{col}(end)];
        elseif length(CopBETtbl.Properties.VariableNames{col})==26
            CopBETtbl.Properties.VariableNames{col} = ['Lebedev15_mrspecific_roi','0',CopBETtbl.Properties.VariableNames{col}(end-1:end)];
        end
    elseif ~isempty(regexp(CopBETtbl.Properties.VariableNames{col},'Tagliazucchi14_windowsize'))
        if length(CopBETtbl.Properties.VariableNames{col})==27
            CopBETtbl.Properties.VariableNames{col} = ['Tagliazucchi14_windowsize','0',CopBETtbl.Properties.VariableNames{col}(end-1:end)];
        end
    elseif ~isempty(regexp(CopBETtbl.Properties.VariableNames{col},'Lebedev16_scale'))
        if length(CopBETtbl.Properties.VariableNames{col})==21
            CopBETtbl.Properties.VariableNames{col} = ['Lebedev16_scale',CopBETtbl.Properties.VariableNames{col}(16),'_ROI0',CopBETtbl.Properties.VariableNames{col}(end)];
        end
    end
    
end

for i = 1:height(CopBETtbl) %without pa
    
    sesid = regexp(data.rsfmri_location{CopBETtbl.sesidx(i)},'rest','once');
    ses = data.rsfmri_location{CopBETtbl.sesidx(i)}(sesid:end);
    
    rp = table2array(readtable(['/mrdata/np2/p3/entropy/data/rp_all/rp',data.rsfmri_raw_ID{CopBETtbl.sesidx(i)},'_',ses]));
    
    rp = rp-mean(rp);
    FDtrans = sum(abs(diff(rp(:,1:3))),2);
    FDrad = 50*pi/180*sum(abs(diff(rp(:,4:6))),2);
    
    CopBETtbl.FD(i) = mean([0;sum([FDtrans,FDrad],2)]);
    
end

CopBETtbl(:,markforremove)=[];
CopBETtbl.data = [];
CopBETtbl.rp = [];
CopBETtbl.entropy = [];
CopBETtbl.num_vols = [];
CopBETtbl.sesidx = [];
save('/mrdata/np2/p3/entropy/CopBET/CopBETtbl_trimmed','CopBETtbl')
writetable(CopBETtbl,'/mrdata/np2/p3/entropy/CopBET/CopBETtbl_trimmed.csv')











