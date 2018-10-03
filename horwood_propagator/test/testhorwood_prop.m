% TESTHORWOOD_PROP ... 
%  
%   ... 
% In Horwoods propagator the Kappa parameter is assumed to be constant
%% AUTHOR    : Niladri Das 
%% $DATE     : 03-Oct-2018 11:33:38 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : testhorwood_prop.m 
clc;clear;close;
startGVMequi.ipmu = [7136.6358*10^3;0;0;0;0]; % in Kilometers
startGVMequi.ipP = diag([20^2*10^6;10^-6;10^-6;10^-6;10^-6]); % in (Kilometer)^2;
startGVMequi.ipalpha = 0;
startGVMequi.ipbeta = zeros(5,1);
startGVMequi.ipGamma = zeros(5,5);
startGVMequi.ipkappa = 3.282806*10^7;
addpath(genpath('main_propagator'))
initialize_propagator;
%% REFERENCE TIME (this is not the start time for the simulation)
Mjd_UTC = 57844; % Reference time.
AuxParam.Mjd_UTC = Mjd_UTC; 
acceleration_function = @(t,mee_states) accel_mee(t,mee_states,eopdata,...
    AuxParam,MJD_J2000,GM_Earth,Arcs,Cnm,Snm);
%%
start_mjdutc = 57844.0270833333; % in MJDUTC
Tstart = (start_mjdutc-AuxParam.Mjd_UTC)*86400; % in Seconds w.r.t reference MJD-UTC
D1 = (1/(96*2)); % Days from start_mjdutc
end_mjdutc = start_mjdutc + D1;% in MJDUTC
Tend = (end_mjdutc-AuxParam.Mjd_UTC)*86400; % in Seconds w.r.t reference MJD-UTC
% In horwoods work the time is varied from 05*100 minutes to 8*100 minutes,
% where 100 minutes is one orbital period.
% The integrator needs time in seconds.
%
[endGVM] = horwood_prop(Tstart,Tend,startGVMequi,acceleration_function,options,GM_Earth);
% ===== EOF ====== [testhorwood_prop.m] ======  
