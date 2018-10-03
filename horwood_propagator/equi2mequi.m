function [mequi] = equi2mequi(equi)
% EQUI2MEQUI: Converts equinoctial orbital elements to modified equinoctial
% orbital elements.
% Input can be a single vector of size 6x1 or
% set of vectors as a matrix with size 6xN
%  
%   ... 

%% AUTHOR    : Niladri Das 
%% $DATE     : 27-Sep-2018 13:07:19 $ 
%% DEVELOPED : 9.3.0.713579 (R2017b) 
%% FILENAME  : equi2mequi.m 
    N = size(equi,2);
    n = size(equi,1);
    if n ~= 6
        error('Check sample size');
    end
    mequi = zeros(6,N);
    for i = 1:N
        mequi(2,i) =  equi(3,i);
        mequi(3,i) =  equi(2,i);
        mequi(4,i) =  equi(5,i);
        mequi(5,i) =  equi(4,i);
        mequi(6,i) =  equi(6,i);
        e2 = equi(2,i)^2+equi(3,i)^2;
        mequi(1,i) = equi(1,i)*(1-e2);
    end
end
% ===== EOF ====== [equi2mequi.m] ======  
