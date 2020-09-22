function [rocheRadius] = calculateRocheRadius(M1,M2)
% Roche radius following Eggleton 1983.
% rocheRadius as seen by star 1, e.g. donor star.
q = M1/M2;
rocheRadius = 0.49./(0.6+(log(1+q.^(1.0./3))./(q.^(2.0./3))));

end