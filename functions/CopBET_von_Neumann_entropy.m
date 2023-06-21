% out = CopBET_von_Neumann_entropy(in,keepdata)
%
% Copenhagen Brain Entropy Toolbox: Von Neumann entropy
% Evaluates XX
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before clustering, e.g.,
%   different subjects or scan sessions.
%   keepdata: Indicates whether the output table also should contain the
%   input data, i.e., by adding an extra column containing entropy values.
%   Defaults to true
%
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

function out = CopBET_von_Neumann_entropy(in,keepdata)

if nargin<2
    keepdata = true;
elseif nargin<1
   error('Please specify input data')
end
if keepdata
    if any(strcmp(in.Properties.VariableNames,'entropy'))
        warning('Overwriting entropy column in datatbl')
    end
end

if ~istable(in)
    if ismatrix(in)
        % convert matrix to table with one entry
        tbl = table;
        tbl.in{1} = in;
        in = tbl;
    else
        error(['Please specify the input data as either a matrix (nxp, n>1)', ...
            'or a table of matrices tbl where the FIRST column contains the data',...
            'with a matrix for each row'])
    end
end

%load data
entropy = nan(height(in),1);
for ses = 1:height(in)
    
    R = corrcoef(in{ses,1}{1});
%     R = corr(in{ses,1}{1});
    
    rho = R./size(R,1);
    
%     entropy(ses) = -trace(rho*logm(rho));
    
    eigvals = eig(rho);
    
    % sometimes the eigenvalues are machine-precision negative...
    eigvals(eigvals<0)=eps;
    entropy(ses) = -sum(eigvals.*log(eigvals));
    
    % Normalize by maximum entropy
    entropy(ses) = entropy(ses)./log(size(R,1));
    
end
if keepdata
    out = in;
    out.entropy = entropy;
else
    out = table;
    out.entropy = entropy;
end
end