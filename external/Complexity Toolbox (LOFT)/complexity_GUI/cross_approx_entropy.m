function s = cross_approx_entropy(m,r,u,v)

% Function to compute Cross-ApEn
% Inputs:
% - m : window length (scalar)
% - r : regularity (scalar)
% - a : first time series (a 1D row vector)
% - b : second time series (a 1D row vector)
%
% Example:
% crossApen = cross_approx_entropy(5,1,a,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check inputs for consistency %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in case these are not row vectors, convert them to row vectors
u = double(u(:)');
v = double(v(:)');
nNoMatch = 0;
if length(u) ~= length(v)
    error('time series a and b must be the same length!');
end

N = length(u);

if m+1 > N
    error('window size cannot be larger than N!')
end

% r = round(0.25*std(u));     % fixed at 25% of standard deviation of time series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute PHI(m) and PHI(m+1) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute PHI(m)
% Pre-compute yj vectors for faster execution
y = zeros(N-m+1,m);
for i = 0:m-1
    y(:,i+1)= v(1+i:N-m+1+i);
end

for i = 1:N-m+1
    % Compute C_i^m(r)(a||b)
    xi = u(i:i+m-1);
    % Count number of j <= N-m+1 such that d[xi,yj]<=r
    count = 0;
    for j = 1:m
        diffXY(:,j) = abs(y(:,j)-xi(1,j));
    end
    d = max(diffXY,[],2);   % max for each row
    count = length(find(d<=r));        
    Cm(i) = count/(N-m+1);
end
ind = find(Cm>0);
if (length(ind)>0)
    phi_m = mean(log(Cm(ind)));
else
    phi_m = 0;
    nNoMatch = nNoMatch+1;
end
clear diffXY d ind

% Compute PHI(m+1)
% Pre-compute yj vectors for faster execution
y1 = zeros(N-(m+1)+1,m+1);
for i = 0:m
    y1(:,i+1)= v(1+i:N-(m+1)+1+i);
end

for i = 1:N-(m+1)+1
    % Compute C_i^m(r)(a||b)
    xi = u(i:i+(m+1)-1);
    % Count number of j <= N-m+1 such that d[xi,yj]<=r
    count = 0;
    for j = 1:m+1
        diffXY(:,j) = abs(y1(:,j)-xi(1,j));
    end
    d = max(diffXY,[],2);   % max for each row
    count = length(find(d<=r));
    Cm1(i) = count/(N-(m+1)+1);
end
ind = find(Cm1>0);
if length(ind)> 0
    phi_m1 = mean(log(Cm1(ind)));
else
    phi_m1 = 0;
    nNoMatch = nNoMatch+1;
end

% Cross-ApEn calculation = log(mean(Cm)./mean(Cm1));
crossApen = phi_m-phi_m1;

s = [crossApen nNoMatch];