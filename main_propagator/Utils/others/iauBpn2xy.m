%   - - - - - - - - - -
%    i a u B p n 2 x y
%   - - - - - - - - - -
% 
%   Extract from the bias-precession-nutation matrix the X,Y coordinates
%   of the Celestial Intermediate Pole.
% 
%   This function is part of the International Astronomical Union's
%   SOFA (Standards Of Fundamental Astronomy) software collection.
% 
%   Status:  support function.
% 
%   Given:
%      rbpn              celestial-to-true matrix (Note 1)
% 
%   Returned:
%      x,y               Celestial Intermediate Pole (Note 2)
% 
%   Notes:
% 
%   1) The matrix rbpn transforms vectors from GCRS to true equator (and
%      CIO or equinox) of date, and therefore the Celestial Intermediate
%      Pole unit vector is the bottom row of the matrix.
% 
%   2) The arguments x,y are components of the Celestial Intermediate
%      Pole unit vector in the Geocentric Celestial Reference System.
% 
%   Reference:
%      "Expressions for the Celestial Intermediate Pole and Celestial
%      Ephemeris Origin consistent with the IAU 2000A precession-
%      nutation model", Astron.Astrophys. 400, 1145-1154
%      (2003)
% 
%      n.b. The celestial ephemeris origin (CEO) was renamed "celestial
%           intermediate origin" (CIO) by IAU 2006 Resolution 2.
% 
%   This revision:  2010 January 18
% 
%   SOFA release 2012-03-01
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = iauBpn2xy(rbpn)

% Extract the X,Y coordinates.
x = rbpn(3,1);
y = rbpn(3,2);

