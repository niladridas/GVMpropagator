% HORWOODPAPERMAIN ... 
%  
%   ... 

%% AUTHOR    : Niladri Das 
%% $DATE     : 27-Sep-2018 12:55:32 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : horwoodpapermain.m 
tic
clc;clear;close;
horwoodInput;
[cansigma,sigma,wt,A] = gensigmaGVM(ipmu,ipP,ipalpha,ipbeta,ipGamma,ipkappa);
sigmamee =  equi2mequi(sigma);
% Sigma points are all in equinoctial orbital element
% Our integrator takes modified equinoctial orbital elements.
% One unit of orbital period is assumed to be 100 minutes
addpath(genpath('main_propagator'))
initialize_propagator;
%% REFERENCE TIME (this is not the start time for the simulation)
Mjd_UTC = 57844; % Reference time.
%
AuxParam.Mjd_UTC = Mjd_UTC; 
acceleration_function = @(t,mee_states) accel_mee(t,mee_states,eopdata,...
    AuxParam,MJD_J2000,GM_Earth,Arcs,Cnm,Snm);
%% SIMULATION START TIME 
start_mjdutc = 57844.0270833333; % in MJDUTC
Tstart = (start_mjdutc-AuxParam.Mjd_UTC)*86400; % in Seconds w.r.t reference MJD-UTC
D1 = (1/(96*2)); % Days from start_mjdutc
end_mjdutc = start_mjdutc + D1;% in MJDUTC
Tend = (end_mjdutc-AuxParam.Mjd_UTC)*86400; % in Seconds w.r.t reference MJD-UTC
% In horwoods work the time is varied from 05*100 minutes to 8*100 minutes,
% where 100 minutes is one orbital period.
% The integrator needs time in seconds.
sigmameeprop = zeros(size(sigmamee,1),size(sigmamee,2));
% Propagated mee sigma points
earth_example; hold on;
reftime = AuxParam.Mjd_UTC;
for j = 1:size(sigmamee,2)
    [Tlist,X_states] = radau(acceleration_function,[Tstart Tend],sigmamee(:,j),options); 
    disp(j);
    sigmameeprop(:,j) = X_states(end,:)';
    X_statesECEF = zeros(size(X_states,1),3);
    for i = 1:size(X_states,1)
     [r_ECI,v] = mee2eci(GM_Earth,X_states(i,:)');
     MJD = (Tlist(i)/86400)+reftime;
     JD = MJD+2400000.5;
     r = ECItoECEF(JD,r_ECI');
     X_statesECEF(i,:) = [r(1),r(2),r(3)];%
    end
    p = plot3(X_statesECEF(:,1),X_statesECEF(:,2),X_statesECEF(:,3)); 
    p.Color = [.49 1 .63];
    p.LineWidth = 2;
    hold on;
%     keyboard
end
hold off;
%% Collect the propagated particles
%% Recover approximating Posterior GVM from the propagated sigma points
% Sample mean of the first 5 states
% Convert the propagated sigmameeprop to sigmaprop
sigmaprop =  mequi2equi(sigmameeprop);
save('sigmaprop.mat','sigmaprop');
save('wt.mat','wt');
% muapprox = bsxfun(@times,wt',sigmaprop(1:5,:));
covxequi = zeros(5,5);
for i = 1:size(sigmamee,2)
    covxequi = covxequi + wt(i)*sigmaprop(1:5,i)*sigmaprop(1:5,i)';
end
Atilde = chol(covxequi,'lower');
%% The mu and alpha of the approximating GVM is 
muapprox = sigmaprop(1:5,1);
alphaapprox = sigmaprop(6,1);
% All in equi frame and not in mequi frame
%% Calculate the beta of the approximating GVM
deltat = Tend - Tstart; % In seconds
a0 = muapprox(1,1);
n0 = sqrt(GM_Earth/a0^3);
temp2 = ipbeta+A'*[-1.5*n0/(a0*deltat);0;0;0;0];
betaapprox = Atilde\A*temp2;
%% Calculate the Gamma of the approximating GVM
temp3 = zeros(5,5);
temp3(1,1) = 15*n0/(4*a0^2*deltat);
temp4 = ipGamma+A'*temp3*A;
Gammaapprox = Atilde\A*temp4*A'*inv(Atilde)';
Gammaapproxvec = vec(Gammaapprox);
%% Recover Optimal Posterior GVM from the propagated sigma points
% This step calculates the alpha, beta and Gamma from the best available
% alpha, beta and Gamma calculated in the previous step
% We need all the sigmaprop samples in this optimization problem 
optfunction1 = @(x) alphabetaGamma(x,sigmaprop,cansigma,Atilde,muapprox,ipkappa);
% x = fminunc(optfunction1,[alphaapprox;betaapprox;Gammaapproxvec]);
x = fminunc(optfunction1,[alphaapprox;betaapprox;Gammaapproxvec]);
alphaoptimal = x(1,1);
betaoptimal = x(2:6,1);
Gammaoptimal = reshape(x(7:31,1),[5,5]);
toc;
% ===== EOF ====== [horwoodpapermain.m] ====== %
