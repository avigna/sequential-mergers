function calculateGravitationalWaveKick(q)
% Calculate gravitational-wave kick from non-spinning binary black-hole merger
% Based on Gonzalez+2007

etaKick = @(q) q/(1+q)^2;

A = 1.2*10^4;
B = -0.93;

vKick = @(eta) (1.2*10^4)*eta*eta*sqrt(1-(4*eta))*(1+(-0.93*eta));

vKick(etaKick(q))
end