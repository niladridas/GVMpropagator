%  - - - - - - - -
%   i a u E o r s
%  - - - - - - - -
%
%  Equation of the origins, given the classical NPB matrix and the
%  quantity s.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rnpb          classical nutation x precession x bias matrix
%     s             the quantity s (the CIO locator)
%
%  Returned (function value):
%                   the equation of the origins in radians.
%
%  Notes:
%  1)  The equation of the origins is the distance between the true
%      equinox and the celestial intermediate origin and, equivalently,
%      the difference between Earth rotation angle and Greenwich
%      apparent sidereal time (ERA-GST).  It comprises the precession
%      (since J2000.0) in right ascension plus the equation of the
%      equinoxes (including the small correction terms).
%
%  2)  The algorithm is from Wallace & Capitaine (2006).
%
% References:
%     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
%     Wallace, P. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
%
%  This revision:  2008 May 26
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eo = iauEors(rnpb, s)

% Evaluate Wallace & Capitaine (2006) expression (16).
x = rnpb(3,1);
ax = x / (1 + rnpb(3,3));
xs = 1 - ax * x;
ys = -ax * rnpb(3,2);
zs = -x;
p = rnpb(1,1) * xs + rnpb(1,2) * ys + rnpb(1,3) * zs;
q = rnpb(2,1) * xs + rnpb(2,2) * ys + rnpb(2,3) * zs;

if ((p ~= 0) || (q ~= 0))
    eo = s - atan2(q, p);
else
    eo = s;
end

