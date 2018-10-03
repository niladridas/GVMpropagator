function [endGVM] = horwood_prop(Tstart,Tend,startGVMequi,acceleration_function,options,GM_Earth)
% HORWOOD_PROP ... 
%  
%   ...
% In Horwoods propagator the Kappa parameter is assumed to be constant

%% AUTHOR    : Niladri Das 
%% $DATE     : 03-Oct-2018 11:29:05 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : horwood_prop.m 
    [cansigma,sigma,wt,A] = gensigmaGVM(startGVMequi.ipmu,startGVMequi.ipP,...
        startGVMequi.ipalpha,startGVMequi.ipbeta,startGVMequi.ipGamma,startGVMequi.ipkappa);
    sigmamee =  equi2mequi(sigma);
    sigmameeprop = zeros(size(sigmamee,1),size(sigmamee,2));
    for j = 1:size(sigmamee,2)
        [~,X_states] = radau(acceleration_function,[Tstart Tend],sigmamee(:,j),options); 
        sigmameeprop(:,j) = X_states(end,:)';
    end
    %% Collect the propagated particles
    %% Recover approximating Posterior GVM from the propagated sigma points
    % Sample mean of the first 5 states
    % Convert the propagated sigmameeprop to sigmaprop
    sigmaprop =  mequi2equi(sigmameeprop);
    muxequi = bsxfun(@times,wt',sigmaprop(1:5,:));
    muxequi = sum(muxequi,2);
    covxequi = zeros(5,5);
    for i = 1:size(sigmamee,2)
        covxequi = covxequi + wt(i)*sigmaprop(1:5,i)*sigmaprop(1:5,i)';
    end
    Atilde = chol(covxequi,'lower');
    %% The mu and alpha of the approximating GVM is 
    muapprox = muxequi;%sigmaprop(1:5,1);
    alphaapprox = sigmaprop(6,1);
    % All in equi frame and not in mequi frame
    %% Calculate the beta of the approximating GVM
    deltat = Tend - Tstart; % In seconds
    a0 = muapprox(1,1);
    n0 = sqrt(GM_Earth/a0^3);
    temp2 = startGVMequi.ipbeta+A'*[-1.5*n0/(a0*deltat);0;0;0;0];
    betaapprox = Atilde\A*temp2;
    %% Calculate the Gamma of the approximating GVM
    temp3 = zeros(5,5);
    temp3(1,1) = 15*n0/(4*a0^2*deltat);
    temp4 = startGVMequi.ipGamma+A'*temp3*A;
    Gammaapprox = Atilde\A*temp4*A'*inv(Atilde)';
    Gammaapproxvec = vec(Gammaapprox);
    %% Recover Optimal Posterior GVM from the propagated sigma points
    % This step calculates the alpha, beta and Gamma from the best available
    % alpha, beta and Gamma calculated in the previous step
    % We need all the sigmaprop samples in this optimization problem 
    optfunction1 = @(x) alphabetaGamma(x,sigmaprop,cansigma,Atilde,muapprox,startGVMequi.ipkappa);
    % x = fminunc(optfunction1,[alphaapprox;betaapprox;Gammaapproxvec]);
    options = optimoptions('fminunc','Display','off');
    x = fminunc(optfunction1,[alphaapprox;betaapprox;Gammaapproxvec],options);
    endGVM.alphaoptimal = x(1,1);
    endGVM.betaoptimal = x(2:6,1);
    endGVM.Gammaoptimal = reshape(x(7:31,1),[5,5]);
    endGVM.muoptimal = muapprox;
    endGVM.Poptimal = covxequi;
    endGVM.kappaoptimal = startGVMequi.ipkappa;
end
% ===== EOF ====== [horwood_prop.m] ======  
