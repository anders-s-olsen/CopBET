function [tbl,data,opts] = CopBET_load_data(dataset,atlas,ts_ROI2ROI)

if strcmp(dataset,'LSD')
    [tbl,data,opts] = CopBET_CarhartHarris_2016_data(atlas,ts_ROI2ROI);
elseif strcmp(dataset,'NRU')
    [tbl,data,opts] = CopBET_McCulloch_Olsen_2023_data(atlas,ts_ROI2ROI);
end