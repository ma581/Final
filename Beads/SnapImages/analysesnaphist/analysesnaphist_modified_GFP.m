function analysesnaphist_modified

D = dir('*1.mat');  %Default '*1.mat'
prefixes = [];
presub = [];
for i = 1:length(D),
name = D(i).name;
L(i) = length(name);
p = findstr('-1.mat',name);
prefixes{i} = name(1:p+1);
presub{i} = name(1:p-1);
end;



for d = 1:length(prefixes)
        plothist2(prefixes{d},presub{d},'cfp','yfp')
end

function plothist2(input,functnm,xaxis,yaxis);
%Order of colors in structure = [cygrx]
load(input);
lardg = eval(functnm);
numcolors = length(lardg(1).mean);
means = reshape([lardg.mean],numcolors,[])';  %Grabs mean data

%Taking data near center of image
cen = reshape([lardg.Centroid],2,[])'; %Grabs centroid data (for filtering below)
dis = sqrt((cen(:,1)-672).^2+(cen(:,2)-512).^2); %d = Distance from center of image
cirIndex = dis<400; %Index those are within 400 pixels of center
meansflat = means(cirIndex,:);  %Filter data on those in center

%Gating the data based on rfp data
gate = 1.0; % A gate of 1 is no gate.
medianrfp = median(means(:,4));
rfpfiltIndex = meansflat(:,4)>(1-gate)*medianrfp...
    & meansflat(:,3)<(1+gate)*medianrfp;
meansflat = meansflat(rfpfiltIndex,:);

% %Filtering data based on rfp
% topbot = 0.20; %Percentage of top and bottom cutoff
% maxrfp = max(means(:,2));
% minrfp = min(means(:,2));
% diff = maxrfp-minrfp;
%
% %Remove anything in the top or bottom
% rfpfiltIndex = meansflat(:,2)>topbot*diff+minrfp...
%     & meansflat(:,2)<maxrfp-topbot*diff;
% meansflat = meansflat(rfpfiltIndex,:);

%Substracting the mode of background values for each cell
back = reshape([lardg.bval],numcolors,[])';
for i = 1:numcolors
    backmode(i) = mode(back(:,i));
    meansflatcorr(:,i) = meansflat(:,i)-backmode(i);
end

%----
%YFP
%---

yfpd = meansflatcorr(:,2);
if min(yfpd)<0
    yfpd = yfpd - min(yfpd);  %Crap code occurs when backmode is too high
end

logy = log10(yfpd);
rfpd = meansflatcorr(:,4);

yfpn = yfpd./rfpd;
ycv = std(yfpd)/mean(yfpd);
ymean = nearest(mean(yfpd))

%----
%GFP
%---
gfpd = meansflatcorr(:,3);
if min(gfpd)<0
    gfpd = gfpd - min(gfpd);  %Crap code occurs when backmode is too high
end


% Marta - modifying the code in order to extract CFP data (2013/06/12)

%----
%CFP
%---

cfpd = meansflatcorr(:,1);
if min(cfpd)<0
    cfpd = cfpd - min(cfpd);  %Crap code occurs when backmode is too high
end

logc = log10(cfpd);
rfpd = meansflatcorr(:,4);

cfpn = cfpd./rfpd;
ccv = std(cfpd)/mean(cfpd);
cmean = nearest(mean(cfpd));
% End of "modifying the code in order to extract CFP data" (2013/06/12)


% figure;
% subplot(2,1,1);
% numbins = 20;
% [bincounts, binpositions] = hist(yfpd,numbins);
% normcounts = bincounts./sum(bincounts)*100; % Sets y axis as percentage of cells
% bar(binpositions, normcounts,1);

% %Begin Gammafit Code
% [parm,parmci] = gamfit(yfpd);
% muhat = parm(1);
% sigmahat = parm(2);
% hold on
% binwidth = binpositions(2) - binpositions(1);
% histarea = binwidth*sum(normcounts);
% x = linspace(binpositions(1),binpositions(end),100);
% y = gampdf(x,muhat,sigmahat);
% plot(x,histarea*y,'r','LineWidth',2)
% %End Gammafit Code

% set(gca,'FontSize',16)
% title(input,'Fontsize',16)
% xlabel('YFP','Fontsize',16);
% ylabel('% Total cells', 'Fontsize',16);
% a = axis;
% % a(1,2) = 2000;
% axis([0 a(1,2) a(1,3) a(1,4)]);
% text(a(1,2)-0.55*a(1,2), a(1,4)-0.1*a(1,4), ['CV ',...
%     sprintf('%6.2f',ycv)],'Fontsize',14);
% text(a(1,2)-0.55*a(1,2), a(1,4)-0.2*a(1,4), ['Mean ',...
%     num2str(ymean)],'Fontsize',14);
% 
% text(a(1,2)-0.55*a(1,2), a(1,4)-0.3*a(1,4), ['Numcell ',...
%     num2str(length(yfpd))],'Fontsize',14); %Checks number of cells

% % SUBPLOT(3,1,2), hist(yfpn,50)
% % title(input,'Fontsize',14)
% % xlabel('YFP normalized to RFP');
% % a = axis;
% % axis([0 a(1,2) a(1,3) a(1,4)]);
% % text(0, a(1,4)-0.1*a(1,4), ['CV = ',...
% %     num2str(std(yfpn)/mean(yfpn))],'Fontsize',12);
% % text(0, a(1,4)-0.2*a(1,4), ['Mean = ', num2str(mean(yfpn))],'Fontsize',12);

% subplot(2,1,2), hist(rfpd,50)
% xlabel('RFP');
% a = axis;
%  axis([0 a(1,2) a(1,3) a(1,4)]);
% text(0, a(1,4)-0.1*a(1,4), ['CV = ',...
%     num2str(std(rfpd)/mean(rfpd))],'Fontsize',12);
% text(0, a(1,4)-0.2*a(1,4), ['Mean = ', num2str(mean(rfpd))],'Fontsize',12);

% Marta - again
figure;
subplot(2,2,1)
numbins = 20;
[bincountsc, binpositionsc] = hist(cfpd,numbins);
normcountsc = bincountsc./sum(bincountsc)*100; % Sets y axis as percentage of cells
bar(binpositionsc, normcountsc,1);

%Begin Gammafit Code
[parmc,parmcic] = gamfit(cfpd);
muhatc = parmc(1);
sigmahatc = parmc(2);
hold on
binwidthc = binpositionsc(2) - binpositionsc(1);
histareac = binwidthc*sum(normcountsc);
xc = linspace(binpositionsc(1),binpositionsc(end),100);
yc = gampdf(xc,muhatc,sigmahatc);
plot(xc,histareac*yc,'r','LineWidth',2)
%End Gammafit Code

set(gca,'FontSize',16)
title(input,'Fontsize',16)
xlabel('CFP','Fontsize',16);
ylabel('% Total cells', 'Fontsize',16);
a = axis;
% a(1,2) = 2000;
axis([0 a(1,2) a(1,3) a(1,4)]);
text(a(1,2)-0.55*a(1,2), a(1,4)-0.1*a(1,4), ['CV ',...
    sprintf('%6.2f',ccv)],'Fontsize',14);
text(a(1,2)-0.55*a(1,2), a(1,4)-0.2*a(1,4), ['Mean ',...
    num2str(cmean)],'Fontsize',14);

text(a(1,2)-0.55*a(1,2), a(1,4)-0.3*a(1,4), ['Numcell ',...
    num2str(length(cfpd))],'Fontsize',14); %Checks number of cells


subplot(2,2,2)
numbins = 20;
[bincounts, binpositions] = hist(yfpd,numbins);
normcounts = bincounts./sum(bincounts)*100; % Sets y axis as percentage of cells
bar(binpositions, normcounts,1);

%Begin Gammafit Code
[parm,parmci] = gamfit(yfpd);
muhat = parm(1);
sigmahat = parm(2);
hold on
binwidth = binpositions(2) - binpositions(1);
histarea = binwidth*sum(normcounts);
x = linspace(binpositions(1),binpositions(end),100);
y = gampdf(x,muhat,sigmahat);
plot(x,histarea*y,'r','LineWidth',2)
%End Gammafit Code

set(gca,'FontSize',16)
title(input,'Fontsize',16)
xlabel('YFP','Fontsize',16);
ylabel('% Total cells', 'Fontsize',16);
a = axis;
% a(1,2) = 2000;
axis([0 a(1,2) a(1,3) a(1,4)]);
text(a(1,2)-0.55*a(1,2), a(1,4)-0.1*a(1,4), ['CV ',...
    sprintf('%6.2f',ycv)],'Fontsize',14);
text(a(1,2)-0.55*a(1,2), a(1,4)-0.2*a(1,4), ['Mean ',...
    num2str(ymean)],'Fontsize',14);

text(a(1,2)-0.55*a(1,2), a(1,4)-0.3*a(1,4), ['Numcell ',...
    num2str(length(yfpd))],'Fontsize',14); %Checks number of cells

subplot(2,2,3), hist(rfpd,50)
xlabel('RFP');
a = axis;
 axis([0 a(1,2) a(1,3) a(1,4)]);
text(0, a(1,4)-0.1*a(1,4), ['CV = ',...
    num2str(std(rfpd)/mean(rfpd))],'Fontsize',12);
text(0, a(1,4)-0.2*a(1,4), ['Mean = ', num2str(mean(rfpd))],'Fontsize',12);

subplot(2,2,4), hist(gfpd,50)
xlabel('GFP');
a = axis;
 axis([0 a(1,2) a(1,3) a(1,4)]);
text(0, a(1,4)-0.1*a(1,4), ['CV = ',...
    num2str(std(gfpd)/mean(gfpd))],'Fontsize',12);
text(0, a(1,4)-0.2*a(1,4), ['Mean = ', num2str(mean(gfpd))],'Fontsize',12);


% figure
% subplot(2,2,1:2)
% text(.5,.5,'subplot(2,2,1:2)',...
%     'FontSize',14,'HorizontalAlignment','center')
% subplot(2,2,3)
% text(.5,.5,'subplot(2,2,3)',...
%     'FontSize',14,'HorizontalAlignment','center')
% subplot(2,2,4)
% text(.5,.5,'subplot(2,2,4)',...
%     'FontSize',14,'HorizontalAlignment','center')
% Marta - end

print('-djpeg',[input,'jpg'])

%%%%%%%%%%% Below for gamma fits.
% figure
%
% [parm,parmci] = gamfit(yfpd)
% [parm2,parmci2] = gamfit(rfpd)
% muhat = parm(1)
% sigmahat = parm(2)
%  numbins = 50;
%  means = yfpd
% %  hist(means,numbins)
%  hold on
%  [bincounts,binpositions] = hist(means,numbins);
%  binwidth = binpositions(2) - binpositions(1);
%  histarea = binwidth*sum(bincounts);
%  x = binpositions(1):0.001:binpositions(end);
% y = gampdf(x,muhat,sigmahat);
% plot(x,histarea*y,'r','LineWidth',2)





