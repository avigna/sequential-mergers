function figure2_massGap
tic;

%
minBHmass = 2.5;
maxBHmass = 300;
lowTotalMassStevenson = 80;
lowGap = 43;
uppGap = 124;
massRatioLimit = 0.9;

% Loading LIGO posteriors
% GW170729
GW170729_90 = importdata('../data/GW170729_mirroredM1M2_90percent_contour.dat');
GW170729_90_m1 = GW170729_90(:,1);
GW170729_90_m2 = GW170729_90(:,2);

% GW190521
GW190521_90 = importdata('../data/GW190521_mirroredM1M2_90percent_contour.dat');
GW190521_90_m1 = GW190521_90(:,1);
GW190521_90_m2 = GW190521_90(:,2);
% Creating histogram
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
totalMass_temp = [totalMassBBHin_CHE, totalMassBBHin_nonCHE];
massRatio_temp = [CHE.massRatioMergingBBHsCHE, nonCHE.massRatioMergingBBHsNonCHE];

totalMass = totalMass_temp(massRatio_temp >= massRatioLimit);

xBins = linspace(0,90,19);
[Ntotal,edgesTotal] = histcounts(totalMass,xBins);
NtotalNorm = Ntotal./max(Ntotal);
edgesTotal = edgesTotal(2:end) - (edgesTotal(2)-edgesTotal(1))/2;

% 
massTripleCompanion = logspace(log10(minBHmass),log10(maxBHmass),2000);
totalMassInnerBBH = logspace(log10(2*minBHmass),log10(2*maxBHmass),2000);
[X,Y] = meshgrid(massTripleCompanion,totalMassInnerBBH);
totalMassTripleMerger = X+Y;

minAlpha = 0.2;
colorTotalMassTripleMerger = zeros(size(totalMassTripleMerger));

for num=1:length(edgesTotal)-1;
colorTotalMassTripleMerger(find(Y >= edgesTotal(num) & Y < edgesTotal(num+1) & totalMassTripleMerger<=2*uppGap & totalMassTripleMerger>=lowTotalMassStevenson)) = NtotalNorm(num)+minAlpha;
end
% colorTotalMassTripleMerger(find(totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=1;
% colorTotalMassTripleMerger(find(Y<uppGap & X<lowGap & totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=minAlpha;
% colorTotalMassTripleMerger(find(totalMassTripleMerger>80 & totalMassTripleMerger<=100))=medAlpha;
colorTotalMassTripleMerger(find(Y>uppGap & Y<2*uppGap & totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=minAlpha;
% colorTotalMassTripleMerger(find(X>uppGap & X<2*uppGap & Y<lowGap & totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=minAlpha;

% Blank inside the mass gap
colorTotalMassTripleMerger(find(Y>lowGap*2 & Y<uppGap+minBHmass))=0.0;
colorTotalMassTripleMerger(find(X>lowGap & X<uppGap))=0.0;

regionColor = ones(size(totalMassTripleMerger));
regionColor(find(X>Y))=2;
regionColor(find(X<lowGap & Y>uppGap))=3;
regionColor(find(X<Y & X>uppGap))=5;

% Visually delimitate regions
colorTotalMassTripleMerger(find(totalMassTripleMerger>99 & totalMassTripleMerger<101))=0;
% colorTotalMassTripleMerger(find(Y>lowGap-0.5 & Y<lowGap+0.5 & totalMassTripleMerger > uppGap))=0;


% Plot
sz=10;
fs=18;
lw=2.0;
alphaNum0=0.7;
alphaNum1=0.2;

stringX = '$\rm \textit{M}_{BH,3}\ [M_{\odot}]$';
stringY = '$\rm \textit{M}_{BBH,in}\ [M_{\odot}]$';

% Based on Du Buisson + 2020
% BH masses 43 < M/Msol < 124
N = 100;
xPISNGap = linspace(lowGap,uppGap+0.4,N);              % Based on Du Buisson+2020
Xfull = linspace(minBHmass,maxBHmass,N);
Y1 = repmat([2*minBHmass;2*maxBHmass],1,N);                   % Based on Stevenson+2019
Y2 = repmat([2*lowGap;uppGap+minBHmass],1,N);

X2=massTripleCompanion;
Z2=ones(size(X2));
totalMass2 = Z2.*2*uppGap;
totalMass3 = Z2.*lowGap;
maxGapFarmer=56;
totalMassFarmer = Z2.*maxGapFarmer;

color=lines(7);
% Plot
numTop = 3;

clf
hold on

s=surf(Y,X,regionColor,'HandleVisibility','Off');
set(s, 'EdgeColor', 'none');
set(s, 'AlphaData', colorTotalMassTripleMerger, 'AlphaDataMapping', 'none');
set(s, 'FaceAlpha', 'flat');

indexBot = find(X2 <= lowGap);
indexTop = find(X2 >= uppGap);
plot3(totalMass2(indexBot),X2(indexBot),Z2(indexBot),'--k','Linewidth',lw,'HandleVisibility','Off')
plot3(totalMass2(indexTop),X2(indexTop),Z2(indexTop),'--k','Linewidth',lw)

% plot3(totalMass3(indexBot),X2(indexBot),2.*Z2(indexBot),'-.k','Linewidth',lw)
% plot3(totalMass3(indexTop),X2(indexTop),2.*Z2(indexTop),'-.k','Linewidth',lw,'HandleVisibility','Off')

purpleColor = [0.4660    0.6740    0.1880];

plot3(totalMassFarmer(indexBot),X2(indexBot),2.*Z2(indexBot),'-.k','Linewidth',lw)
plot3(totalMassFarmer(indexTop),X2(indexTop),2.*Z2(indexTop),'-.k','Linewidth',lw,'HandleVisibility','Off')

plot3(GW170729_90_m1,GW170729_90_m2,numTop.*ones(size(GW170729_90_m2)),'Color',purpleColor,'Linewidth',lw)
plot3(GW190521_90_m1,GW190521_90_m2,numTop.*ones(size(GW190521_90_m2)),'--','Color',purpleColor,'Linewidth',lw)

h1=fill([Y1(1,:) flip(Y1(2,:))], [xPISNGap flip(xPISNGap)],'k','EdgeColor','none','HandleVisibility','Off');
set(h1,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area       

h2=fill([Y2(1,:) flip(Y2(2,:))],[Xfull flip(Xfull)],'k','EdgeColor','none','HandleVisibility','Off');
set(h2,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area   

text(10,140,numTop,'\bf IMRI','FontSize',fs,'Color','k','Interpreter','Latex')
text(65,140,numTop,'\bf A','FontSize',fs,'Color','k','Interpreter','Latex')
text(72,35,numTop,'\bf B','FontSize',fs,'Color','k','Interpreter','Latex')
text(63,18,numTop,'\bf U','FontSize',fs,'Color','k','Interpreter','Latex')
text(140,20,numTop,'\bf C','FontSize',fs,'Color','k','Interpreter','Latex')
ht=text(150,10,numTop,'\bf IMRI','FontSize',fs,'Color','k','Interpreter','Latex');
set(ht,'Rotation',270);

text(6,75,numTop,'\bf Mass gap','FontSize',fs,'Color','w','Interpreter','Latex')
hgap=text(105,10,numTop,'\bf Mass gap','FontSize',fs,'Color','w','Interpreter','Latex');
set(hgap,'Rotation',270);

text(300,10,numTop,'\bf CHE','FontSize',fs,'Color','k','Interpreter','Latex')
text(300,200,numTop,'\bf CHE','FontSize',fs,'Color','k','Interpreter','Latex')

plot3(49.5,33.5,numTop,'dk','MarkerFaceColor','k','MarkerSize',sz)

% blueString = '$\rm \textit{M}_{BH,out} < \textit{M}_{BBH,in}$';
% redString = '$\rm \textit{M}_{BH,out} > \textit{M}_{BBH,in}$';
% yellowString = '$\rm \textit{M}_{BBH,in} \gg \textit{M}_{BH,out}$';
% 
% text(10,5,numTop,blueString,'FontSize',fs,'Color','k','Interpreter','Latex')
% text(6,30,numTop,redString,'FontSize',fs,'Color','k','Interpreter','Latex')
% ht2=text(400,30,numTop,yellowString,'FontSize',fs,'Color','k','Interpreter','Latex')
% set(ht2,'Rotation',270);

% GW170729MBBHin = 49.5;
% GW170729MBBH3 = 33.5;
% plot3(GW170729MBBHin,GW170729MBBH3,numTop,'*k','MarkerSize',10)

legend( 'CHE limit',...
        'Single BH limit',...
        'GW170729',...
        'GW190521',...
        'Triple',...
        'Interpreter','Latex',...
        'Location','SouthWest',...
        'FontSize',fs,...
        'Box','Off');


xlabel(stringY,'Interpreter','Latex','FontSize',fs)
ylabel(stringX,'Interpreter','Latex','FontSize',fs)
ylim([minBHmass maxBHmass])
xlim([2*minBHmass 2*maxBHmass])
yticks([2.5 10 43 124 200 300]);
xticks([5 10 43 86 124+minBHmass 248 600]);

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'Fontsize',fs)

box on
hold off
colormap(lines(5))

print(gcf,'../figures/figure2_massGap.png','-dpng','-r400');

toc;
end