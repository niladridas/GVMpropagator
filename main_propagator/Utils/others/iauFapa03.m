%  - - - - - - - - - -
%   i a u F a p a 0 3
%  - - - - - - - - - -
%
%  Fundamental argument, IERS Conventions (2003):
%  general accumulated precession in longitude.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     t         TDB, Julian centuries since J2000.0 (Note 1)
%
%  Returned (function value):
%               general precession in longitude, radians (Note 2)
%
%  Notes:
%
%  1) Though t is strictly TDB, it is usually more convenient to use
%     TT, which makes no significant difference.
%
%  2) The expression used is as adopted in IERS Conventions (2003).  It
%     is taken from Kinoshita & Souchay (1990) and comes originally
%     from Lieske et al. (1977).
%
%  References:
%     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
%     48, 187
%     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
%     Astron.Astrophys. 58, 1-16
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2009 December 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = iauFapa03(t)

% General accumulated precession in longitude.
a = (0.024381750 + 0.00000538691 * t) * t;

