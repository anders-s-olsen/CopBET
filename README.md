# Copenhagen Brain Entropy Toolbox

![License](https://img.shields.io/badge/License-MIT-blue.svg)

The Copenhagen Brain Entropy Toolbox is a MATLAB toolbox that provides a collection of functions for evaluating 12 different entropy metrics described in the paper "Navigating the chaos of psychedelic neuroimaging: A multi-metric evaluation of acute psilocybin effects on brain entropy" by Drummond McCulloch, Anders S Olsen et al. This toolbox allows researchers to analyze and quantify the entropy of fMRI-recorded brain activity using various methods.

## Installation

To use the Copenhagen Brain Entropy Toolbox, follow the steps below:

1. Clone this repository to your local machine:

   ```shell
   git clone https://github.com/your-username/copenhagen-brain-entropy-toolbox.git
   ```

2. Add the toolbox directory to your MATLAB path using the `addpath` function:

   ```matlab
   addpath(genpath('/path/to/copenhagen-brain-entropy-toolbox');)
   ```

## Usage

The toolbox provides a set of MATLAB functions that can be used to evaluate different entropy metrics on brain data. The metrics take either a matrix (time by regions) or a table in which the first column of the table contains matrices of size (time by regions). The output will always be a table with a column with the name "entropy". The functions have an option to keep the input table information in the output, thus adding a new entropy-column to the existing table. 
Note that two functions need voxel-wise input data and thus, the input data for these should be either a char or a table or chars of the path to the (denoised) 4D-volumes for the corresponding subject. These functions also require an atlas (in the same space!) to be specified. Here's an example of how to use the toolbox:

```matlab
% Load your brain data (i.e., time series)
data1 = load('brain_data1.mat'); %subject 1
data2 = load('brain_data2.mat'); %subject 2

% Structure the data in a table:
tbl = table;
tbl.data{1} = data1;
tbl.data{2} = data2;
tbl.subject = [1;2];

% Compute entropy using a specific metric, specifying to add entropy
% as a new column in the existing table and to parallelize computation.
tbl = CopBET_FUNCTION(tbl,atlas 'keepdata',true,'parallel',true);

```

Replace `'CopBET_FUNCTION'` with the name of the desired entropy metric:
* CopBET_DCC_entropy
* CopBET_degree_distribution_entropy
* CopBET_diversity_coefficient
* CopBET_geodesic_entropy
* CopBET_intranetwork_synchrony
* CopBET_LEiDA_transition_entropy
* CopBET_metastate_series_complexity
* CopBET_motif_connectivity_entropy
* CopBET_sample_entropy
* CopBET_temporal_entropy
* CopBET_time_series_complexity
* CopBET_von_Neumann_entropy

The functions that require voxel-wise input data and a specified atlas are `'CopBET_intranetwork_synchrony'`, `'CopBET_motif_connectivity_entropy'`, and `'CopBET_sample_entropy'`. Please refer to the paper for more details on available metrics and their usage. Note that while some of these functions return scalar entropy values per row in the input table, some return a vector or matrix of entropy values, which need to be further processed afterwards. For example, `CopBET_degree_distribution_entropy` returns an entropy value for each "mean degree". 

## Example script

We have posted a script `CopBET_main_CH2016data` showing examples of how to use the functions with the openly available acute IV LSD dataset https://openneuro.org/datasets/ds003059. In the script, the dataset is assumed to be present as a folder in the CopBET folder. 

## Atlases

We provided a folder with all the atlases that we used in our paper to replicate previous studies, in 2mm MNI-152 space. However, for the analyses in our paper we stripped away cerebellar ROIs due to inconsistent inclusion in the field of view. These atlases may be used in conjunction with the LSD dataset. 

## Contributing

Please email Anders S. Olsen (ansol@dtu.dk) if you would like to contribute to the toolbox or otherwise have comments and questions to the code not answered in the paper. 

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more information.

## References

If you use the Copenhagen Brain Entropy Toolbox in your research, please cite the following paper:

Drummond McCulloch, Anders S Olsen, et al. "Navigating the chaos of psychedelic neuroimaging: A multi-metric evaluation of acute psilocybin effects on brain entropy". *XX*, 2023.

Please replace "YYYY" with the publication year of the paper.
