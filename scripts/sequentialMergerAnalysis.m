function sequentialMergerAnalysis(R,M1,M2,M3,sBBH,separationInner,separationOuter,printValues)
% Script for sequential merger analysis from Vigna-Gomez+2020b
% 
% [R] = Rsol
% [M] = Msol
% R = 9;  
% M1 = 40;
% M2 = M1;
% M3 = 140;
% frad = 0.05;
% spinBBH = 0.68;
factorL2 = 1.32;    % From Marchant+2016

Minner = M1+M2;
% Estimate fraction of radiated mass from GWs
if M1>=M2
	qin = M2/M1;
else
    qin = M1/M2;
end
tempLoader = importdata('../data/Mrad_fraction_chi0.dat');
q_frac = tempLoader.data(:,1);
Mrad_frac = tempLoader.data(:,2);
[minVal,minIdx] = min(abs(qin-q_frac));
frad = Mrad_frac(minIdx);
MBBHin = Minner*(1-frad);
finalTotalMass = MBBHin+M3;

if MBBHin/M3 <= 1
    finalMassRatio =MBBHin/M3;
else
    finalMassRatio = M3/MBBHin;
end

chirpMass = calculateChirpMass(MBBHin,M3);
chiEffective = calculateChiEffectiveAlignedSpins(MBBHin,sBBH,M3,0);

%
aoutOverainCrit = calculateTripleStabilityCriteria(Minner,M3,0.0,0.0);
rocheRadius = max(calculateRocheRadius(M1,M2),calculateRocheRadius(M2,M1));
outerLagrangianPoint = factorL2*rocheRadius;

ain_minimum = R/outerLagrangianPoint;
separationInnerMinimum = round(ain_minimum);
aouter_minimum = aoutOverainCrit.*separationInnerMinimum;

if separationInner > 0
    aouter = separationInner.*aoutOverainCrit;
else
    aouter = -1;
end


maxSeparationPristine = calculateSeparationLimits(Minner,M3);
maxSeparation = calculateSeparationLimits(MBBHin,M3);

[efinal,aFactorChange] = calculateBlaauwKick(Minner,MBBHin,M3);
minOuterSeparation = aFactorChange*aouter_minimum;

if separationOuter > 0
    maxOuterSeparation = aFactorChange*separationOuter;
else
	maxOuterSeparation = -1;
end
    
if printValues
    Minner
    qin
    frad
    MBBHin
    M3
    finalMassRatio
    chirpMass
    chiEffective    
    R
    rocheRadius
    outerLagrangianPoint
    aoutOverainCrit
    separationInnerMinimum
    aouter_minimum
    separationInner    
    separationOuter
    aouter
    efinal
    minOuterSeparation
    maxOuterSeparation
    maxSeparation

end

end