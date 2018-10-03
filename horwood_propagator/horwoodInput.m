%% AUTHOR    : Niladri Das 
%% $DATE     : 27-Sep-2018 11:58:30 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : horwoodInput.m 
% ipmu = [7136.635;0;0;0;0]; % in Kilometers
ipmu = [7136.6358*10^3;0;0;0;0]; % in Kilometers
ipP = diag([20^2*10^6;10^-6;10^-6;10^-6;10^-6]); % in (Kilometer)^2;
ipalpha = 0;
ipbeta = zeros(5,1);
ipGamma = zeros(5,5);
ipkappa = 3.282806*10^7;
% ===== EOF ====== [horwoodInput.m] ======  
