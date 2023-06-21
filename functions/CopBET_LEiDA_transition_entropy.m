% out = CopBET_LEiDA_transition_entropy(in,K,varargin)
%
% Copenhagen Brain Entropy Toolbox: LEiDA transition entropy.
% Evaluates the leading eigenvectors of instantaneous phase coherence
% matrices, clusters the data using k-means with squared euclidean distance
% function, and evaluates the Markov transition rate entropy, all according
% to the specifications in Kringelbach et al., 2020: "Dynamic coupling of
% whole-brain neuronal and neurotransmitter systems" and the associated
% Github repository.
%
% Input:
%   in: a matrix (nxp,n>1) or a table where the first column contains
%   matrices (in cells) to be concatenated before clustering, e.g.,
%   different subjects or scan sessions.
%   K: number of clusters for k-means clustering. 
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

% potential tests:
% Check that cluster centroids make sense (qualitatively) and correct size
% Check unit norm of data points
% Check nans/infs

function out = CopBET_LEiDA_transition_entropy(in,K,varargin)

if nargin<2
    error('Please specify a number of components, K')
end

[out,numworkers,in] = CopBET_function_init(in,varargin);

% Do LEiDA and concatenate data for clustering
disp('Running LEiDA')
datasizes = nan(height(in),1);
for ses = 1:height(in)
    datasizes(ses) = size(in{ses,1}{1},1);
end
Lead_eigs_all = nan(sum(datasizes),size(in{1,1}{1},2));
c = 1;
for ses = 1:height(in)
    tmp = in{ses,1}{1};
    tmp = tmp - mean(tmp);
    Lead_eigs_all(c:c+datasizes(ses)-1,:) = LEiDA(tmp)';
    sensible_data_check(Lead_eigs_all(c:c+datasizes(ses)-1,:),['Leading eigenvectors, session',num2str(ses)]);
    c = c+datasizes(ses);
end


disp('Running k-means')
if numworkers>0
    a = gcp;
    if isempty(a)
        parpool(numworkers)
    end
    [IDX,C]=kmeans(Lead_eigs_all,K,'Distance','sqeuclidean',...
        'Replicates',200,'Display','final','Options',statset('UseParallel',1));
else
    [IDX,C]=kmeans(Lead_eigs_all,K,'Distance','sqeuclidean',...
        'Replicates',200,'Display','final');
end

entropy = tpm_entropy(IDX,datasizes);
sensible_data_check([entropy{:}],'entropy');
out.entropy = [entropy{:}]';



end

%% functions
function Leading_Eig = LEiDA(ts)
% More or less copy-pasted from https://github.com/decolab/pnas-neuromod/blob/master/LEiDA_PsiloData.m
% This function can be optimized heavily, however, we tried to keep it as
% much as the original to prevent differences from the original paper.



% fnq=1/(2*TR);                 % Nyquist frequency
% flp = 0.04;                    % lowpass frequency of filter (Hz)
% fhi = 0.07;                    % highpass
% Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
% k=2;                          % 2nd order butterworth filter
% [bfilt,afilt]=butter(k,Wn);   % construct the filter

[Tmax,N_areas] = size(ts);
Phase_BOLD=zeros(N_areas,Tmax); %OBS transposed
%%%%%%%%%%%%% from PNAS code
for seed=1:N_areas
    tstmp=detrend(ts(:,seed));
%     signal_filt =filtfilt(bfilt,afilt,tstmp);
%     Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    Phase_BOLD(seed,:) = angle(hilbert(tstmp));
    % Phase_BOLD = angle(hilbert(ts))';
end
Leading_Eig = nan(N_areas,Tmax);

for t=1:Tmax
    
    %Calculate the Instantaneous FC (BOLD Phase Synchrony)
    iFC=zeros(N_areas);
    for n=1:N_areas
        for p=1:N_areas
            iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
        end
    end
    
    % Get the leading eigenvector
    
    [V1,~]=eigs(iFC,1);
    % Make sure the largest component is negative
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end
    % Save V1 from all frames in all fMRI sessions in Leading eig
    Leading_Eig(:,t)=V1; %vertcat(V1,V2);
end
end

function entropy = tpm_entropy(IDX,datasizes)
% More or less copy-pasted from https://github.com/decolab/pnas-neuromod/blob/master/LEiDA_PsiloData.m

K = numel(unique(IDX(IDX>0)));

tpts_traversed = 1;
% entropy = nan(length(datasizes),1);
for ses = 1:length(datasizes)
    T = tpts_traversed:tpts_traversed+datasizes(ses)-1;
    
    Ctime=IDX(T);
    
    PTRANSITION=zeros(K,K);
    for c1=1:K
        for c2=1:K
            sumatr=0;
            for t=1:length(Ctime)-1
                if Ctime(t)==c1 && Ctime(t+1)==c2
                    sumatr=sumatr+1;
                end
            end
            if ~isempty((find(Ctime(1:length(Ctime)-1)==c1)))
                PTRANSITION(c1,c2)=sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
            end
        end
    end
    
    entropy{ses} = EntropyMarkov(PTRANSITION);
    
    tpts_traversed = tpts_traversed + datasizes(ses);
end
end

function H = EntropyMarkov(P)
number_states=size(P,1);
[V,D] = eig( P' );    % eigendecomp of P' (P=transition matrix)
[~,ii]= max(diag(D)) ;
st = V(:,ii);
p1 = abs(st)/sum(abs(st)); % p(i) (stationary!)

% Remark: you could also get an estimate of p(i) directly from the observations
% p(i) = sum( state == i )/T;

% Markov entropy:
Hi = zeros(1,number_states);
for row=1:number_states
    Hi(row) = -p1(row)*sum( (P(row,:)+eps).*log2(P(row,:)+eps) );
end

H = sum(Hi)/log2(number_states);
end





