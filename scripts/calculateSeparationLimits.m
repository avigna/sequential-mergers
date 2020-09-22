function [amaxRsol] = calculateSeparationLimits(M1,M2)
% Function based on eq. 5.10 from Peters 1964
% [m1]=[m2]=Msol
% amax = Rsol
% e=0

% Constants in cgs
Gcgs = 6.6743*10^-8;
ccgs = 2.99792458*10^10;
HubbleTimeSeconds = 4.55*10^17;
cmToRsol = 1/(6.957*10^10);

% Convert mass variables to grams 
m1g = M1.*1.989*10^33;
m2g = M2.*1.989*10^33;

amax = ((256./5).*(Gcgs^3).*(ccgs^-5).*m1g.*m2g.*(m1g+m2g).*HubbleTimeSeconds).^(1/4);
amaxRsol = amax*cmToRsol;

end