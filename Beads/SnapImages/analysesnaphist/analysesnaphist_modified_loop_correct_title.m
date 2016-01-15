function analysesnaphist_modified_loop_correct_title(varargin)
global savedata;
D = dir('*1.mat');  %Default '*1.mat'
prefixes = [];
presub = [];

%-------------------------------------------------------------------------
%Gets filenames
%--------------------------------------------------------------------------
for i = 1:length(D),
name = D(i).name;
L(i) = length(name);
p = findstr('-1.mat',name);
prefixes{i} = name(1:p+1);
presub{i} = name(1:p-1);
end;

%-------------------------------------------------------------------------
%Reading input variables; so far just wheter to export txt files
%--------------------------------------------------------------------------
savedata=0; %default setting

%Checking input
if length(varargin)>0
    for i = 1:length(varargin),
        theparam = lower(varargin{i});
        switch(theparam),
            case 'exportdata',
                savedata = 1;
                disp('Data is exported');
        end
    end
else
    disp('Data is not exported')
end

%-------------------------------------------------------------------------
%Gets number of colour channels to analyse
%--------------------------------------------------------------------------
load(prefixes{1});
lardg = eval(presub{1});
numcolors = length(lardg(1).mean);
disp(num2str(numcolors));

%-------------------------------------------------------------------------
%Gets used colors
%--------------------------------------------------------------------------
colorsuffixes={'*-c.tif','*-y.tif','*-g.tif','*-r.tif'};
colornames={'CFP','YFP','GFP','RFP'};
usedcolors={};
for i=1:length(colorsuffixes)
    if length(dir(colorsuffixes{i})) >0;
        usedcolength=length(usedcolors);
        usedcolors{usedcolength+1}=colornames{i};
        
    end
end

%-------------------------------------------------------------------------
%Generating Graphs
%--------------------------------------------------------------------------

for d = 1:length(prefixes)
        plothist2(prefixes{d},presub{d},colornames,usedcolors)
end

%--------------------------------------------------------------------------
%Defining plothist2 function 
%--------------------------------------------------------------------------

function plothist2(input,functnm,colornames,usedcolors);
%--------------------------------------------------------------------------
%Setting the default values
%--------------------------------------------------------------------------
%global savedata;




%--------------------------------------------------------------------------
%Loading inputdata
%Order of colors in structure = [cygrx]
%--------------------------------------------------------------------------
load(input);
input2=input;
lardg = eval(functnm);
numcolors = length(lardg(1).mean); %gets number of colors which is equivalent to the number of collums in lardg(1)
means = reshape([lardg.mean],numcolors,[])';  % rehapes lardge into a table with two columns and Grabs mean data

%--------------------------------------------------------------------------
%Taking only cells which are 400 pixel far away from the center
%--------------------------------------------------------------------------
%cen = reshape([lardg.Centroid],2,[])'; %reshapes centeroid matrix into a table with two colums (x and y)Grabs centroid data (for filtering below)
%dis = sqrt((cen(:,1)-696).^2+(cen(:,2)-520).^2); %d = Distance from center of image
%cirIndex = dis<400; %Index those are within 400 pixels of center (gives out a logical matrix)
%meansflat = means(cirIndex,:);  %Filter data on those in center (only intensities (rfp cfp values) of celles close to the middle of the image (up to 400 pixels from the center) surive.
meansflat = means;
%--------------------------------------------------------------------------
%Gating RFP Q: What does it do besides calc the median rfp
%--------------------------------------------------------------------------

gate = 1.0; % A gate of 1 is no gate.
medianrfp = median(means(:,2));
rfpfiltIndex = meansflat(:,2)>(1-gate)*medianrfp...
    & meansflat(:,2)<(1+gate)*medianrfp;
meansflat = meansflat(rfpfiltIndex,:);

%--------------------------------------------------------------------------
%Substracting the mode of background values for each cell
%--------------------------------------------------------------------------
back = reshape([lardg.bval],numcolors,[])'; % arranges background (bval) into a two collum array
for i = 1:numcolors
    backmode(i) = mode(back(:,i)); % gets the most frequent value out of back
    meansflatcorr(:,i) = meansflat(:,i)-backmode(i); %subtracts the background from the meansflat (rfp and cfp values)
end

%--------------------------------------------------------------------------
%Choosing plot positions depending on #of colors
%--------------------------------------------------------------------------
plotpos={{[1 2] [3 4]},{[1,2],3,4},{1,2,3,4}};
numcolors;
if numcolors == 2
usedplotpos=plotpos{1,1};
  elseif numcolors == 3
usedplotpos=plotpos{1,2};
  elseif numcolors == 4
usedplotpos=plotpos{1,3};
  else
disp('Too much imput!')
end

%--------------------------------------------------------------------------
%Loop for different colors
%--------------------------------------------------------------------------

figure;
for i=1:numcolors
    fpd = meansflatcorr(:,i);
    if min(fpd)<0
        fpd = fpd - min(fpd);  %Crap code occurs when backmode is too high
    end
    logfp = log10(fpd);
    
rfpd = meansflatcorr(:,numcolors);

fpn = fpd./rfpd;
cv = std(fpd)/mean(fpd);
meanfp = nearest(mean(fpd));

subplot(2,2,usedplotpos{i}) % defines were plot goes on figure
numbins = 40;
[bincounts, binpositions] = hist(fpd,numbins);
normcounts = bincounts./sum(bincounts)*100; % Sets y axis as percentage of cells
bar(binpositions, normcounts,1);

dogammafit=0;

switch usedcolors{i}
    case 'YFP'
        dogammafit=1;
    case 'CFP'
        dogammafit=1;
end
%--------------------------------------------------------------------------
%Gammafit 
%--------------------------------------------------------------------------
if dogammafit == 1
    
    %Begin Gammafit Code
[parm,parmci] = gamfit(fpd);
muhat = parm(1);
sigmahat = parm(2);
hold on
binwidth = binpositions(2) - binpositions(1);
histarea = binwidth*sum(normcounts);
x = linspace(binpositions(1),binpositions(end),100);
y = gampdf(x,muhat,sigmahat);
plot(x,histarea*y,'r','LineWidth',2)
%End Gammafit Code
end

% Setting title correct

input=input(2:end-2);
r=strfind(input,'_');
for j=1:length(r)
    input(r(j))=' ';
end


%--------------------------------------------------------------------------
%Making Graphs
%--------------------------------------------------------------------------

set(gca,'FontSize',16)

if i==1
title(input,'Fontsize',16)
end

xlabel(usedcolors{i},'Fontsize',16);
ylabel('% Total cells', 'Fontsize',16);
a = axis;
% a(1,2) = 2000;
axis([0 a(1,2) a(1,3) a(1,4)]);
text(a(1,2)-0.30*a(1,2), a(1,4)-0.1*a(1,4), ['CV ',...
    sprintf('%6.2f',cv)],'Fontsize',14);
text(a(1,2)-0.30*a(1,2), a(1,4)-0.2*a(1,4), ['Mean ',...
    num2str(meanfp)],'Fontsize',14);

text(a(1,2)-0.30*a(1,2), a(1,4)-0.3*a(1,4), ['Numcell ',...
    num2str(length(fpd))],'Fontsize',14); %Checks number of cells
        
end

print('-djpeg',[input2,'-400','.jpg']);









