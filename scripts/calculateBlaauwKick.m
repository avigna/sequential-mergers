function [efinal,aFactorChange] = calculateBlaauwKick(Minitial,Mfinal,Mcompanion)
% Following Postnov & Yungelson 2014, page 32
% Blaauw kick from spherically symmetric mass loss
% Following Vigna-Gomez+2020b (eqs. 2 and 3):
% Minitial    = M_{BH,1} + M_{BH,2}
% Mfinal      = (1-f_{rad})(M_{BH,1} + M_{BH,2}) = M_{BBH,in}
% Mcompanion  = M_{BH,3}
chi = (Minitial+Mcompanion)./(Mfinal+Mcompanion);
aFactorChange = 1./(2-chi);
efinal = (Minitial-Mfinal)./(Mfinal+Mcompanion);

end