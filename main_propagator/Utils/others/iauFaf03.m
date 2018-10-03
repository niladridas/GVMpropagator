%  - - - - - - - - -
%   i a u F a f 0 3
%  - - - - - - - - -
%
%  Fundamental argument, IERS Conventions (2003):
%  mean longitude of the Moon minus mean longitude of the ascending
%  node.
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
%               F, radians (Note 2)
%
%  Notes:
%  1) Though t is strictly TDB, it is usually more convenient to use
%     TT, which makes no significant difference.
%
%  2) The expression used is as adopted in IERS Conventions (2003) and
%     is from Simon et al. (1994).
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
%     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
%
%  This revision:  2009 December 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = iauFaf03(t)

global const

% Mean longitude of the Moon minus that of the ascending node
% (IERS Conventions 2003).
a = mod( 335779.526232 + ...
        t * ( 1739527262.8478 + ...
        t * (       - 12.7512 + ...
        t * (        - 0.001037 + ...
        t * (          0.00000417 ) ) ) ), const.TURNAS ) * const.DAS2R;

