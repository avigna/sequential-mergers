function figure4_plotQuantityOfInterest
tic;

minBHmass = 2.5;
maxBHmass = 300;
lowGap = 43;
uppGap = 124;

massTripleCompanion = logspace(log10(minBHmass),log10(maxBHmass),2000);
totalMassInnerBBH = logspace(log10(2*minBHmass),log10(2*maxBHmass),2000);
[X,Y] = meshgrid(massTripleCompanion,totalMassInnerBBH);
totalMassTripleMerger = X+Y;
massRatioTripleMerger = X./Y;
massRatioTripleMergerInverse = Y./X;
massRatioTripleMerger(find(massRatioTripleMerger>1))=massRatioTripleMergerInverse(find(massRatioTripleMerger>1));

% Color PISN mass gap regions of interest
minAlpha = 0.3;
colorTotalMassTripleMerger = zeros(size(totalMassTripleMerger)).*minAlpha;
colorTotalMassTripleMerger(find(totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=1;
colorTotalMassTripleMerger(find(totalMassTripleMerger>80 & totalMassTripleMerger<=100))=1;
colorTotalMassTripleMerger(find(Y>uppGap & Y<2*uppGap & totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=1;
colorTotalMassTripleMerger(find(X>uppGap & X<2*uppGap & Y<lowGap & totalMassTripleMerger>100 & totalMassTripleMerger<2*uppGap))=1;
colorTotalMassTripleMerger(find(Y>lowGap*2 & Y<uppGap+minBHmass))=0.0;
colorTotalMassTripleMerger(find(X>lowGap & X<uppGap))=0.0;

% Calculate quantities of interests
chiEff = 0.68./(1+(X./Y));
amaxRsol = calculateSeparationLimits(X,Y);
aoutOverainCrit = calculateTripleStabilityCriteria(Y,X,0.0,0.0);
amaxRsol_cropped = amaxRsol;
amaxRsol_cropped(find(log10(amaxRsol)<1.5)) = 30;
amaxRsol_cropped(find(log10(amaxRsol)>log10(160))) = 150;

% Visually delimitate regions
colorTotalMassTripleMerger(find(totalMassTripleMerger>99 & totalMassTripleMerger<101))=0;
% colorTotalMassTripleMerger(find(Y>lowGap-0.5 & Y<lowGap+0.5 & totalMassTripleMerger > uppGap))=0;

% totalMassTripleMerger(jj)

% Choose what to plot
prompt = 'Choose to plot:\n (a) mass ratio \n (b) max outer separation \n (c) chi effective \n (d) triple critical stability criteria \n';
chooseWhatToPlot = input(prompt,'s');

if chooseWhatToPlot == 'a'
    quantityOfInterest = massRatioTripleMerger;
    numTop=1;
    plotString = '../figures/figure4_massRatio.png'
    cbarLimits = [0 1]
    cbarTicks = [0:0.2:1];
    cbarLabelString = '$q$'
elseif chooseWhatToPlot == 'b'
    quantityOfInterest = log10(amaxRsol_cropped);
    numTop=3;
    plotString = '../figures/figure4_maxSeparation.png'
    cbarLimits = log10([30 150])
    cbarTicks = log10([30, 50, 100, 130, 150]);
    cbarTickLabels = {'30','50','100','130','150'};
    cbarLabelString = '$\rm \textit{a}_{out,max}\ [R_{\odot}]$'
%     quantityOfInterest = (amaxRsol_cropped);
%     numTop=3;
%     plotString = '../Plots/figure4_maxSeparation.png'
%     cbarLimits = ([30 150])
%     cbarTicks = ([30, 50, 100, 130, 150]);
%     cbarTickLabels = {'30','50','100','130','150'};
%     cbarLabelString = '$\rm \textit{a}_{out,max}\ [R_{\odot}]$'
elseif chooseWhatToPlot == 'c'
    quantityOfInterest = chiEff;
    numTop=3;
    plotString = '../figures/figure4_chiEffective.png'
    cbarLimits = [0 0.7]
    cbarTicks = [0:0.1:0.7];
    cbarLabelString = '$\chi_{\rm{eff}}$';   
elseif chooseWhatToPlot == 'd'
    quantityOfInterest = aoutOverainCrit;
    numTop=20;
    plotString = '../figures/figure4_stabilityCriteria.png'
    cbarLimits = [2.8 14.5]
    cbarTicks = [3:14];
    cbarLabelString = '$\rm (\textit{a}_{out}/\textit{a}_{in})|_{\rm crit}$'
else
    warning('No reasonable choice to choose what to plot')
end

% Plot
fs=20;
lw=2;
alphaNum0=0.7;

stringX = '$\rm \textit{M}_{BH,3}\ [M_{\odot}]$';
stringY = '$\rm \textit{M}_{BBH,in}\ [M_{\odot}]$';

% Based on Du Buisson + 2020
% BH masses 43 < M/Msol < 124
N = 100;
xPISNGap = linspace(lowGap,uppGap+0.4,N);              % Based on Du Buisson+2020
Xfull = linspace(minBHmass,maxBHmass,N);
Y1 = repmat([2*minBHmass;2*maxBHmass],1,N);                   % Based on Stevenson+2019
Y2 = repmat([2*lowGap;uppGap+minBHmass+0.4],1,N);

X2=massTripleCompanion;
Z2=ones(size(X2));
totalMass2 = Z2.*2*uppGap;
maxGapFarmer=56;
totalMassFarmer = Z2.*maxGapFarmer;

clf
hold on
% Color plot
s=surf(Y,X,quantityOfInterest);
set(s, 'EdgeColor', 'none');
set(s, 'AlphaData', colorTotalMassTripleMerger, 'AlphaDataMapping', 'none');
set(s, 'FaceAlpha', 'flat');

% Plot vertical lines
indexBot = find(X2 <= lowGap);
indexTop = find(X2 >= uppGap);
plot3(totalMass2(indexBot),X2(indexBot),Z2(indexBot),'--k','Linewidth',lw)
plot3(totalMass2(indexTop),X2(indexTop),Z2(indexTop),'--k','Linewidth',lw,'HandleVisibility','Off')
plot3(totalMassFarmer(indexBot),X2(indexBot),10.*Z2(indexBot),'-.k','Linewidth',lw)
plot3(totalMassFarmer(indexTop),X2(indexTop),10.*Z2(indexTop),'-.k','Linewidth',lw,'HandleVisibility','Off')

% Shade regions
h1=fill([Y1(1,:) flip(Y1(2,:))], [xPISNGap flip(xPISNGap)],'k','EdgeColor','none');
set(h1,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area       
h2=fill([Y2(1,:) flip(Y2(2,:))],[Xfull flip(Xfull)],'k','EdgeColor','none');
set(h2,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area       
   
% Add text
text(10,140,numTop,'\bf IMRI','FontSize',fs,'Color','w','Interpreter','Latex')
text(65,140,numTop,'\bf A','FontSize',fs,'Color','w','Interpreter','Latex')
text(70,35,numTop,'\bf B','FontSize',fs,'Color','w','Interpreter','Latex')
text(62,18,numTop,'\bf U','FontSize',fs,'Color','w','Interpreter','Latex')
text(140,20,numTop,'\bf C','FontSize',fs,'Color','w','Interpreter','Latex')
ht=text(150,10,numTop,'\bf IMRI','FontSize',fs,'Color','w','Interpreter','Latex');
set(ht,'Rotation',270);

% text(6,75,numTop,'\bf Mass gap','FontSize',fs,'Color','w','Interpreter','Latex')
% hgap=text(105,12,numTop,'\bf Mass gap','FontSize',fs,'Color','w','Interpreter','Latex')
% set(hgap,'Rotation',270);

xlabel(stringY,'Interpreter','Latex','FontSize',fs)
ylabel(stringX,'Interpreter','Latex','FontSize',fs)
ylim([minBHmass maxBHmass])
xlim([2*minBHmass 2*maxBHmass])
yticks([2.5 10 43 124 200 300]);
xticks([5 10 43 86 124 248 600]);

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'Fontsize',fs)

cbar = colorbar;
cbar.Label.Interpreter = 'latex';
cbar.Limits = cbarLimits;
cbar.Ticks = (cbarTicks);
cbar.Label.String = cbarLabelString;
cbar.Label.FontSize = fs;

if chooseWhatToPlot == 'b'
    cbar.TickLabels = cbarTickLabels;
end

box on
hold off
colormap(parula(1000))

print(gcf,plotString,'-dpng','-r400');

toc;
end