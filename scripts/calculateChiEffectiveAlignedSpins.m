function [chiEffectiveAligned] = calculateChiEffectiveAlignedSpins(M1,chi1,M2,chi2)
% Effective spin calculation 
% m1,2: masses
% a1,2: dimensionless effective spins aligned with the orbit
% Following Vigna-Gomez+2020b (eq. 5):
% M1 = M_{BBH,in}, chi1 = \chi_{BBH,in}
% M2 = M_{BH,3}, chi2 = \chi_{BH,3}
chiEffectiveAligned = ((M1.*chi1)+(M2.*chi2))./(M1+M2);

end