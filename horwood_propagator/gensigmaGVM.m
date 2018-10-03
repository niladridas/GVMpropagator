function [cansigma,sigma,wt,A] = gensigmaGVM(varargin)
% GENSIGMAGVM generate sigma points from a given GVM distribution 
%  Number of samples are 2n+3 where n is the length of the first argument
%  mu
%   ... 

%% AUTHOR    : Niladri Das 
%% $DATE     : 27-Sep-2018 12:22:50 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : gensigmaGVM.m 
%
% 
mu = varargin{1};
P = varargin{2};
alpha = varargin{3};
beta = varargin{4};
gam = varargin{5};
kap = varargin{6};
%
n = length(mu); % Size of the state space excluding the circular domain
%% Generate sigma points from Canonical GVM
xi = sqrt(3);
I2 = besselj(2,kap);
I1 = besselj(1,kap);
I0 = besselj(0,kap);
B2k = 1-(I2/I0);
B1k = 1-(I1/I0);
eta = acos(0.5*B2k/B1k-1);
cansigma = zeros(n+1,2*n+3);
cansigma(end,2) = eta;
cansigma(end,3) = -eta;
for i = 4:2:(2*n+3)
    cansigma(i/2-1,i) = xi;
    cansigma(i/2-1,i+1) = -xi;
end
%% Generate weights of the sigma points
wt = ones(2*n+3,1);
wt(4:end,1) = (1/6)*wt(4:end,1);
wt(2,1) = B1k^2/(4*B1k-B2k);
wt(3,1) = wt(3,1);
wt(1,1) = 1-2*wt(2,1)-10*(1/6);
%% Generate Sigma points for input GVM from canonical GVM sigma ponits
A = chol(P,'lower');
sigma = zeros(n+1,2*n+3);
sigma(1:n,:) = A*cansigma(1:n,:)+repmat(mu,1,2*n+3);
tmp1 = 0.5*sum((gam*cansigma(1:n,:)).*(gam*cansigma(1:n,:)),1);
tmp2 = beta'*cansigma(1:n,:)+alpha*ones(1,2*n+3);
sigma(end,:) = tmp1 + tmp2 +cansigma(end,:);
% ===== EOF ====== [gensigmaGVM.m] ======  
