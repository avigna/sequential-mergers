function figure1_massParameterSpace_heavyBBHs
tic;

% Loading LIGO posteriors
% GW170729
GW170729_90 = importdata('../data/GW170729_mirroredM1M2_90percent_contour.dat');
GW170729_90_m1 = GW170729_90(:,1);
GW170729_90_m2 = GW170729_90(:,2);

% GW190521
GW190521_90 = importdata('../data/GW190521_mirroredM1M2_90percent_contour.dat');
GW190521_90_m1 = GW190521_90(:,1);
GW190521_90_m2 = GW190521_90(:,2);

minBHmass = 2.5;
maxBHmass = 300;
lowGap = 43;
uppGap = 124;

massTripleCompanion = logspace(log10(minBHmass),log10(maxBHmass),2000);
totalMassInnerBBH = logspace(log10(2*minBHmass),log10(2*maxBHmass),2000);
[X,Y] = meshgrid(massTripleCompanion,totalMassInnerBBH);
totalMassTripleMerger = X+Y;

colorTotalMassTripleMerger = ones(size(totalMassTripleMerger));
colorTotalMassTripleMerger(find(Y>lowGap*2 & Y<uppGap+minBHmass))=0.0;
colorTotalMassTripleMerger(find(X>lowGap & X<uppGap))=0.0;

regionColor = ones(size(totalMassTripleMerger));
regionColor(find(X>Y))=2;
regionColor(find(X<lowGap & Y>uppGap))=3;
regionColor(find(X<Y & X>uppGap))=5;

YGWTC=linspace(2.5,42.5,5000);
ZGWTC=ones(size(YGWTC));

% Plot
sz=10;
fs=18;
lw=2.0;
lw2=3;
alphaNum0=0.7;
alphaNum1=0.2;
numTop = 3;

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


clf
hold on
numTop2 = 3;
s=surf(Y,X,regionColor,'HandleVisibility','Off');
set(s, 'EdgeColor', 'none');
set(s, 'AlphaData', colorTotalMassTripleMerger, 'AlphaDataMapping', 'none');
set(s, 'FaceAlpha', 'flat');

% indexBot = find(X2 <= lowGap);
% indexTop = find(X2 >= uppGap);
% plot3(totalMass2(indexBot),X2(indexBot),Z2(indexBot),'--k','Linewidth',lw,'HandleVisibility','Off')
% plot3(totalMass2(indexTop),X2(indexTop),Z2(indexTop),'--k','Linewidth',lw,'HandleVisibility','Off')

% plot3(GW170729_50_m1,GW170729_50_m2,numTop.*ones(size(GW170729_50_m2)),'w','Linewidth',lw,'HandleVisibility','Off')
plot3(GW170729_90_m1,GW170729_90_m2,numTop.*ones(size(GW170729_90_m2)),'w','Linewidth',lw)

% plot3(GW190521_50_m1,GW190521_50_m2,numTop.*ones(size(GW190521_50_m2)),'--w','Linewidth',lw,'HandleVisibility','Off')
plot3(GW190521_90_m1,GW190521_90_m2,numTop.*ones(size(GW190521_90_m2)),'--w','Linewidth',lw)

h1=fill([Y1(1,:) flip(Y1(2,:))], [xPISNGap flip(xPISNGap)],'k','EdgeColor','none','HandleVisibility','Off');
set(h1,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area       

h2=fill([Y2(1,:) flip(Y2(2,:))],[Xfull flip(Xfull)],'k','EdgeColor','none','HandleVisibility','Off');
set(h2,'FaceAlpha',alphaNum0)% if size(y0,1)==2 %plot shaded area       

blueString = '$\rm \textit{M}_{BH,3}< \textit{M}_{BBH,in}$';
redString = '$\rm \textit{M}_{BBH,in}\ll \textit{M}_{BH,3}$';
yellowString = '$\rm \textit{M}_{BH,3}\ll \textit{M}_{BBH,in}$';

text(7,3,numTop,blueString,'FontSize',fs,'Color','k','Interpreter','Latex')
text(7,200,numTop,redString,'FontSize',fs,'Color','k','Interpreter','Latex')
text(140,10,numTop,yellowString,'FontSize',fs,'Color','k','Interpreter','Latex')

text(6,75,numTop,'\bf Mass gap','FontSize',fs,'Color','w','Interpreter','Latex')
hgap=text(105,10,numTop,'\bf Mass gap','FontSize',fs,'Color','w','Interpreter','Latex');
set(hgap,'Rotation',270);

plot3(49.5,33.5,numTop,'d','Color','w','MarkerFaceColor','w','MarkerSize',sz)

lgd=legend( 'GW170729',...
        'GW190521',...
        'Triple',...
        'Interpreter','Latex',...
        'FontSize',fs,...
        'Position',[0.6632 0.6 0.2279 0.1662],...
        'TextColor','white',...
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

print(gcf,'../figures/figure1_massParameterSpace_heavyBBHs.png','-dpng','-r400');

toc;
end