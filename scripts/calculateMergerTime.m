function [tc] = calculateMergerTime(m1,m2,a)
% Function based on eq. 5.10 from Peters 1964
% [m1]=[m2]=Msol
% amax = Rsol
% e=0.0
Gcgs = 6.6743*10^-8;
ccgs = 2.99792458*10^10;
HubbleTimeSeconds = 4.55*10^17;
RsolTocm = (6.957*10^10);
cmToRsol = 1/RsolTocm;
MyearToSec = (3.154e+7)*(10^6);
secToMyear = 1/MyearToSec;

m1g = m1.*1.989*10^33;
m2g = m2.*1.989*10^33;
acm = a.*RsolTocm;

timeCoalescenceSeconds = (acm.^4)./((256./5).*(Gcgs^3).*(ccgs^-5).*m1g.*m2g.*(m1g+m2g));
tc=timeCoalescenceSeconds.*secToMyear;

end