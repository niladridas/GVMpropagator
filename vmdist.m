function vm = vmdist(theta,alpha,kappa)
% VMDIST ... 
%  
%   ... 
% Support of VM is any interval of 2*pi 
%
%% AUTHOR    : Niladri Das 
%% $DATE     : 25-Sep-2018 14:26:05 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : gvmdist.m 

% disp(' !!!  You must enter code into this file < gvmdist.m > !!!')
num = exp(-2*kappa*(sin(0.5*(theta-alpha))^2));
den = 2*pi*exp(-kappa)*besselj(0,kappa);
vm = num/den;
end
% ===== EOF ====== [vmdist.m] ======  
