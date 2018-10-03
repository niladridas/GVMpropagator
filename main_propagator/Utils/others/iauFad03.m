%  - - - - - - - - -
%   i a u F a d 0 3
%  - - - - - - - - -
%
%  Fundamental argument, IERS Conventions (2003):
%  mean elongation of the Moon from the Sun.
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
%               D, radians (Note 2)
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
function a = iauFad03(t)

global const

% Mean elongation of the Moon from the Sun (IERS Conventions 2003).
a = mod(    1072260.703692 + ...
        t * ( 1602961601.2090 + ...
        t * (        - 6.3706 + ...
        t * (          0.006593 + ...
        t * (        - 0.00003169 ) ) ) ), const.TURNAS ) * const.DAS2R;

