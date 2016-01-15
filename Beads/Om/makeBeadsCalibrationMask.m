
load('beadsCalib100.mat')
sortedBeads = sort(reshape(beadsCalib,1,1392*1040));
top25 = sortedBeads(nearest(end*0.95):end);
figure;
[C,h] = contourf(beadsCalib/mean(top25));set(gca,'YDir','reverse');
figure;
ContourLevel = 0.9
[C,h] = contourf(beadsCalib/mean(top25),[ContourLevel,ContourLevel]);set(gca,'YDir','reverse');
Cnan = [C(1,~isnan(C(1,2:end)));C(2,~isnan(C(1,2:end)))];
figure;
hold on;

ExcludeBorder = 10;
CnanB = Cnan(:,Cnan(1,:)>ExcludeBorder);
CnanB = CnanB(:,CnanB(2,:)>ExcludeBorder);

plot(CnanB(1,2:end),CnanB(2,2:end),'.');
convHull = convhull(CnanB(1,2:end),CnanB(2,2:end));
convHullnTess = convhulln(CnanB(:,2:end)');

plot(CnanB(1,convHull+1),CnanB(2,convHull+1),'g-o')
plot(CnanB(1,convHullnTess(:,1)+1),CnanB(2,convHullnTess(:,1)+1),'r.-')

set(gca,'YDir','reverse');
xlim([1,1392]);
ylim([1,1040]);
hold off;
contmask = poly2mask(CnanB(1,convHull+1),CnanB(2,convHull+1),1040,1392);
figure;
imshow(contmask);
