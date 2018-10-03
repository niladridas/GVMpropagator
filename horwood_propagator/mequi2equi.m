function [equi] = mequi2equi(mequi)
% MEQUI2EQUI: Converts modified equinoctial orbital elements to  equinoctial
% orbital elements.
% Input can be a single vector of size 6x1 or
% set of vectors as a matrix with size 6xN
%  
%   ... 

%% AUTHOR    : Niladri Das 
%% $DATE     : 2-Oct-2018 13:07:19 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : mequi2equi.m 
    N = size(mequi,2);
    n = size(mequi,1);
    if n ~= 6
        error('Check sample size');
    end
    equi = zeros(6,N);
    for i = 1:N
        equi(3,i) = mequi(2,i);  
        equi(2,i) = mequi(3,i) ;
        equi(5,i) = mequi(4,i);
        equi(4,i) = mequi(5,i) ;
        equi(6,i) = mequi(6,i) ;
        e2 = mequi(2,i)^2+mequi(3,i)^2;
        equi(1,i) = mequi(1,i)/(1-e2);
    end
end
% ===== EOF ====== [equi2mequi.m] ======  
