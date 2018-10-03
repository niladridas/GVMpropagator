function [function_val] = alphabetaGamma(x,sigmaprop,cansigma,Atilde,muapprox,ipkappa)
    alphaapproxopt = x(1,1);
    betaapproxopt = x(2:6,1);
    Gammaapproxopt = reshape(x(7:31,1),[5,5]);
    Nsamples = size(sigmaprop,2);
    r = 0;
    for i = 1:Nsamples
        temp5 =  cansigma(1:5,i)'*cansigma(1:5,i);
        temp61 = Atilde\(sigmaprop(1:5,i)-muapprox);
        temp62 = temp61'*temp61;
        temp71 = sigmaprop(6,i)-alphaapproxopt-betaapproxopt'*temp61 - 0.5*temp61'*Gammaapproxopt*temp61;
        temp72 = 4*ipkappa*(sin(0.5*cansigma(6,i))^2-sin(0.5*temp71)^2);
        r = r + (temp5 - temp62 + temp72)^2;
    end
    function_val = r;
end