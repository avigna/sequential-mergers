function figure2_totalMassDistribution_COMPAS
tic;

%
minBHmass = 2.5;
maxBHmass = 300;
lowGap = 43;
uppGap = 124;
massRatioLimit = 0.9;

%
M=importdata('../data/Mrad_fraction_chi0.dat');
massRatioRad = M.data(:,1);
fRad = M.data(:,2);

%
CHE = importdata('../data/totalMassMergingBBHsCHE.mat');
nonCHE = importdata('../data/totalMassMergingBBHsNonCHE.mat');

% Calculating post-merger mass
% CHE channel
[minVal, idx] = min(abs(CHE.massRatioMergingBBHsCHE - massRatioRad));
totalMass_CHE = CHE.totalMassMergingBBHsCHE;
massRatioRad_CHE = massRatioRad(idx)';
fRad_CHE = fRad(idx)';
totalMassBBHin_CHE = totalMass_CHE.*(1.-fRad_CHE);

clear minVal
clear idx
% Non-CHE channel
[minVal, idx] = min(abs(nonCHE.massRatioMergingBBHsNonCHE - massRatioRad));
totalMass_nonCHE = nonCHE.totalMassMergingBBHsNonCHE;
massRatioRad_nonCHE = massRatioRad(idx)';
fRad_nonCHE = fRad(idx)';
totalMassBBHin_nonCHE = totalMass_nonCHE.*(1.-fRad_nonCHE);

%
totalMass = [totalMassBBHin_CHE, totalMassBBHin_nonCHE];
massRatio = [CHE.massRatioMergingBBHsCHE, nonCHE.massRatioMergingBBHsNonCHE];

totalMassCHE = totalMassBBHin_CHE(find(CHE.massRatioMergingBBHsCHE >= massRatioLimit));
totalMassNonCHE = totalMassBBHin_nonCHE(find(nonCHE.massRatioMergingBBHsNonCHE >= massRatioLimit));

%
xBins = linspace(0,90,19);
[Ntotal,edgesTotal] = histcounts(totalMass,xBins);
edgesTotal = edgesTotal(2:end) - (edgesTotal(2)-edgesTotal(1))/2;

[NCHEafter,edgesCHEafter] = histcounts(totalMassCHE,xBins);
[NNonCHEafter,edgesNonCHEafter] = histcounts(totalMassNonCHE,xBins);

Y = zeros(size(Ntotal));
% Plot
sz=18;
fs=18;
lw=2;
alphaNum0=0.8;
alphaNum1=0.4;

stringX = '$\rm \textit{M}_{BH,3}\ [M_{\odot}]$';
stringY = '$\rm \textit{M}_{BBH,in}\ [M_{\odot}]$';

color=lines(7);

clf
% Position = [left bottom width height]
% axes('Position',[0.1300    0.800    0.7750    0.15])
axes('Position',[0.1300    0.2    0.7750    0.3])
hold on

plot(edgesTotal,NNonCHEafter+NCHEafter,'k','Linewidth',lw)

h0=fill([edgesTotal flip(edgesTotal)], [Y flip(NNonCHEafter)],'k','EdgeColor','none');
set(h0,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area       

h1=fill([edgesTotal flip(edgesTotal)], [NNonCHEafter flip(NNonCHEafter+NCHEafter)],'k','EdgeColor','none');
set(h1,'FaceAlpha',alphaNum1)% if size(y0,1)==2 %plot shaded area       

ylabel('\rm Counts','Interpreter','Latex','FontSize',fs)
legend( '$\rm{Post-merger}: 0.9 \le q_{\rm{in}}$',...
        '$\rm{Standard}: 0.9 \le q_{\rm{in}}$',...
        '$\rm{Compact}: 0.9 \le q_{\rm{in}}$',...
        'Box','Off',...
        'Interpreter','Latex',...
        'Location','NorthWest')
ylim([0 1000])
xlim([2*minBHmass 2*maxBHmass])
xticks([5 10 43 86 124+minBHmass 248 600]);
% xticklabels([]);
% yticks([0 0.025 0.05]);

set(gca, 'XScale', 'log')
set(gca,'Fontsize',fs)
box on
hold off

print(gcf,'../figures/figure2_totalMassDistribution_COMPAS.png','-dpng','-r400');

toc;
end