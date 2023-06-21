function s = sample_entropy(m,r,u,a)

% Function to compute multiscale sample entropy
% compute sample entropy, when a = 1
% Inputs:
% - m : window length (scalar)
% - r : regularity (scalar)
% - u : time series (a 1D col vector)
% - a : skipping parameter

% Check inputs for consistency
u = double(u(:)');
N = length(u);
if m+1 > N
    error('window size cannot be larger than series length!');
end

% for multiscale entropy
if a > 1
   n = fix(N/a);
   uu = zeros(1, n);
   for i = 1:n
       uu(i) = sum(u(((i-1)*a+1):(i*a)))/a;
   end
else
    n = N;
    uu = u;
end

% initialization
nNoMatch = 0;
cor = zeros(1,2);
dataMat = zeros(m+1,n-m);
for i = 1:m+1
    dataMat(i,:) = uu(i:n-m+i-1);
end

% calculation
for d = m:m+1
    count = zeros(1,n-m);
    tempMat = dataMat(1:d,:);
    for i = 1:n-d
        dist = max(abs(tempMat(:,i+1:n-m) - repmat(tempMat(:,i),1,n-m-i)));
        D = (dist < r);
        count(i) = sum(D)/(n-m-1);
    end
    cor(d-m+1) = sum(count)/(n-m);
end
    
if cor(2) ~= 0 & cor(1) ~= 0
    sampEn = -log(cor(2)/cor(1));
else
    sampEn = 0;
    nNoMatch = nNoMatch + 1;
end

s = [sampEn nNoMatch];