function [aoutOverainCrit] = calculateTripleStabilityCriteria(min,mout,eout,i)
% Equation 1 from Toonen+2020: https://arxiv.org/pdf/2004.07848.pdf
% Stability criteria from Mardling & Aarseth (1999, 2001)
% See also Vigna-Gomez+2020b
% qout = mout/min = m3/(m1+m2)
% eout: eccentricity outer orbit
% i: inclination

qout = mout./min;
aoutOverainCrit = (2.8./(1-eout)).*(1-(0.3*i/pi)).*((1+qout).*(1+eout)./sqrt(1-eout)).^(2.0/5);

end