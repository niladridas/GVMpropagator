%  - - - - - - - - - -
%   i a u P o m 0 0
%  - - - - - - - - - -
%
%  Form the matrix of polar motion for a given date, IAU 2000.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     xp,yp        coordinates of the pole (radians, Note 1)
%     sp           the TIO locator s' (radians, Note 2)
%
%  Returned:
%     rpom         polar-motion matrix (Note 3)
%
%  Notes:
%  1) The arguments xp and yp are the coordinates (in radians) of the
%     Celestial Intermediate Pole with respect to the International
%     Terrestrial Reference System (see IERS Conventions 2003),
%     measured along the meridians to 0 and 90 deg west respectively.
%
%  2) The argument sp is the TIO locator s', in radians, which
%     positions the Terrestrial Intermediate Origin on the equator.  It
%     is obtained from polar motion observations by numerical
%     integration, and so is in essence unpredictable.  However, it is
%     dominated by a secular drift of about 47 microarcseconds per
%     century, and so can be taken into account by using s' = -47*t,
%     where t is centuries since J2000.0.  The function iauSp00
%     implements this approximation.
%
%  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
%     that it is the final rotation when computing the pointing
%     direction to a celestial source.
%
%  Called:
%     iauIr        initialize r-matrix to identity
%     iauRz        rotate around Z-axis
%     iauRy        rotate around Y-axis
%     iauRx        rotate around X-axis
%
%  Reference:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rpom = iauPom00(xp, yp, sp)

% Construct the matrix.
rpom = zeros(3);
rpom = iauIr(rpom);
rpom = iauRz(sp, rpom);
rpom = iauRy(-xp, rpom);
rpom = iauRx(-yp, rpom);

