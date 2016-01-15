function [L,maskThr]= segfluorToSnapScript(ph3,beadmode)
% function [L, phsub, L9, rect, L8, L7, L6, L5, L4, L3, mesub]= segphase(ph3, p)
% 
% segments phase image
%   1. finds edge
%   2. finds mask
%   3. fills cells using imcomplement
%   4. fills edge using cells found from incomplement
%   5. breaks up cell clumps (cutcurvbpts)
%   6. cuts up individual cells (cutcurv)
%   7. cuts up any long cells remaining (breakcellfluor)
%   8. cuts up any kinky cells remaining (breakcellfluor, dekinker)
%   9. removes any small `cells' left
%  10. relabels image
%
%   L:      segmented image (full size)
%   L9:     segmented image (smaller size)
%   phsub:  phase images (smaller size)
%   rect:   transformation required to reconstruct full size image from smaller ones

%Internalize p structure to this functino:
if beadmode
    p.edge_lapofgauss_sigma = 2;
    p.minCellArea = 50;
    p.minCellLengthConservative = 20;
    p.minCellLength = 30;
    p.maxCellWidth = 10;
    p.maxThreshCut = 0.3;
    p.maxThreshCut2 = 0.2;
    p.maxThresh = 0.25;
    p.minThresh = 0.25;
    p.imNumber1 = 2;
    p.imNumber2 = 1;
    p.radius = 5;
    p.angThresh = 2.7;
    p.mThr = 1.9;
    p.ker = 5;
    p.segmentationPhaseSlice = 1;
    p.numphaseslices = 1;
    p.minNumEdgePixels = 250;  % JCR: value used for nature methods paper
    maskThr = 10; %just a made up number, hack
    p.AreaThresh = 1;
    p.FirstThreshFudge = 0.2;
    p.MajorAxisLength = 16;
    p.Area = 0.5*25*pi;
    p.minCellLengthConservative= 5;
else
    p.edge_lapofgauss_sigma = 2;
    p.minCellArea = 50;
    p.minCellLengthConservative = 20;
    p.minCellLength = 30;
    p.maxCellWidth = 10;
    p.maxThreshCut = 0.3;
    p.maxThreshCut2 = 0.2;
    p.maxThresh = 0.25;
    p.minThresh = 0.25;
    p.imNumber1 = 2;
    p.imNumber2 = 1;
    p.radius = 5;
    p.angThresh = 2.7;
    p.mThr = 1.9;
    p.ker = 5;
    p.segmentationPhaseSlice = 1;
    p.numphaseslices = 1;
    p.minNumEdgePixels = 250;  % JCR: value used for nature methods paper
    maskThr = 10; %just a made up number, hack
    p.AreaThresh = 10;
    p.FirstThreshFudge = 0.17;
    p.MajorAxisLength = 50;
    p.Area = 50;
    p.minCellLengthConservative= 12;
end

emptyflag = 0;

if p.imNumber2>size(ph3,3)
    p.imNumber2 = 1;
end
if p.imNumber1>size(ph3,3)
    p.imNumber1 = size(ph3,3);
end

maxsizeofdirt=20000;

% FILTER AND EDGE PHASE
% median filter and rescale phase images
for i= 1:size(ph3,3),
    ph3(:,:,i)= medfilt2(ph3(:,:,i),[3 3]);
    
    x= double(ph3(:,:,i));
    s= sort(x(:));
    small= s(25);
%     big= s(end-25);
signal = s-(small+40);
signal(signal<=0) = [];
if isempty(signal)
    ph3size = size(ph3);
    signal = zeros(ph3size(1:2));
    emptyflag = 1;
else
    emptyflag = 0;
end
signal = signal+(small+40);
big = signal(round(length(signal)*0.93))
    rescaled=(x - small)/(big - small);
    rescaled(rescaled<0)= 0;
    ph3(:,:,i)= uint16(10000*rescaled);
end;
me= edge(ph3(:,:,p.segmentationPhaseSlice),'log',0,p.edge_lapofgauss_sigma);
col_center=[size(me,1),size(me,2)]/2;

%%
% %===========================================================================
% %This code does not seem to contribute anything to the series of L
% %imagess so I have commented it out, OP Jun 6, 2015
% % FIND MASK (smallest region containing cells)
% 
% % MASK: first iteration
% %   Run box of size boxsize over the edge image 
% %   and count the number of white pixels in the box.
% %   Cells have less than minonpixels pixels.
% 
% imout= zeros(size(ph3(:,:,1)));
% replaceboxsize= 5;  % boxsize is shrunk to replaceboxsize on output image 
% boxsize= 40;
% 
% % was minnopixels, now p.minNumEdgePixels in schnitzcell parameter structure
% %minNumEdgePixels= 300; % original for e.coli
% minNumEdgePixels= 299;%150; OP changed back % note change (from 300 to 250 during bacillus tests)
% % JCR: significance of 300 is that's the value of minNumEdgePixels for e.coli
% % JCR: bacillus is smaller, until 2005-08-18 it was 250 (e.g. for methods 
% %      paper work); since then it's 215
% 
% for x= 1:replaceboxsize:size(me,1) - boxsize
%     for y= 1:replaceboxsize:size(me,2) - boxsize
%         xmin= x + boxsize/2;
%         xmax= x + boxsize/2 + replaceboxsize - 1;
%         ymin= y + boxsize/2;
%         ymax= y + boxsize/2 + replaceboxsize - 1;
%         
%         subregion= me(x:x + boxsize, y:y + boxsize);
%         checkiscell= sum(subregion(:)) < p.minNumEdgePixels;
%         imout(xmin:xmax, ymin:ymax)= ones(replaceboxsize)*checkiscell;
%     end
% end
% % dilate output image to make sure ends of cells included
% imout = uint8(imdilate(imout,strel('disk',replaceboxsize)));
% 
% % tidy up imout so that only one highlighted region remains
% % imout = imopen(imout,strel('disk',25)); % new
% % this waxs replaced by the following 5 lines, to make sure something remains!
% imout_temp = imopen(imout,strel('disk',25)); % new
% if max(imout_temp(:))
%     imout = imout_temp;
% end
% clear imout_temp;
% 
% if p.minNumEdgePixels<300
%     imdilate(imout,strel('disk',10));
% end
% imoutlabel= bwlabel(imout);
% r= regionprops(imoutlabel, 'area', 'centroid');
% 
% % cells are in small groups
% % keep the clump closest to the centre of the image
% cen= [r.Centroid];
% xcen= cen(1:2:end);
% ycen= cen(2:2:end);
% 
% % find distance from centre of image
% distances= sqrt((xcen-col_center(2)).^2 + (ycen-col_center(1)).^2);
% [ds, di]= sort(distances);
% % keep closest highlighted region only
% imout2= zeros(size(imout));
% imout2(imoutlabel == di(1))= 1;
% 
% % keep additional cell clumps if the number of colonies is greater than one (Nitzan)
% for i= 2:length(di);
%     imout2(imoutlabel == di(i))= 1;
% end;
% % 
% % EXTRACT SMALLER SUBREGIONS
% % extra= 5; % <- extra number of pixels on either side of highlighted region
% % [fx,fy]= find(imout2);
% % xmin= max(min(fx) - extra, 1);
% % xmax= min(max(fx) + extra, size(imout2,1));
% % ymin= max(min(fy) - extra, 1);
% % ymax= min(max(fy) + extra, size(imout2,2));
% % extract subregions
% % rect= [xmin ymin xmax ymax];
% % phsub= ph3(xmin:xmax, ymin:ymax,  :) ;
% % imout2sub= imout2(xmin:xmax, ymin:ymax);
% % mesub= me(xmin:xmax, ymin:ymax);
% imout2sub= imout2;
% %===========================================================================
%%

rect = [1 1 size(ph3,1) size(ph3,2)];
phsub = ph3(:,:,1);
mesub=me;

% MASK: second iteration
% use a threshold to find where the cells are
imt = phsub;
imt = imt - maxmax(imt)*p.FirstThreshFudge;%2000*1.5;%OP fudge factor
e = edge(imt,'log',0);
% e = me;
f = imfill(e,'holes');
L = bwlabel(f);


r = regionprops(L,'Area');
flittle = find([r.Area]<p.Area);
%fbigs = find([r.Area]>1200);
L(flittle)=[];
%L(fbigs)=[];

L = bwlabel(f);
goodones = setdiff(1:max(max(L)), flittle);
%goodones = setdiff(goodones, fbigs);
bw2 = ismember(L, goodones);
L2 = bwlabel(bw2);
L2= renumberimage(L2);
r2 = regionprops(L,'Solidity');
fsolid = find([r2.Solidity]>0.2);
bw3 = ismember(L2,fsolid);
L3 = bwlabel(bw3);
imshowlabel(L3);
r= regionprops(L3,'majoraxislength');
flittles= find([r.MajorAxisLength] < p.minCellLengthConservative);
disp(['Removing small cells(', num2str(length(flittles)),').']);
for i= 1:length(flittles),
    % delete cell
   L3(L3 == flittles(i))= 0;
end;

L3= renumberimage(L3);
if max(L3(:))>0
    use_mask_cutting = 3;
else
    use_mask_cutting = use_mask_cutting + 1;
end

L4 = L3;
% CUTTING CELL CLUMPS
% finds branchpoints in cell, i.e. points where several cells are grouped
% together in a clump, by morphologically thinning the cell-segment to a
% single-pixel thin line (thinnng with parameter inf). Such branchpoints
% are cut, disentangling the clumps.
%L5= rmsinglepointconnections(L4);
% disp(['Cutting cells(',num2str(max2(L5)),').']);
% for i= 1:max2(L5)
%     %disp([' cutting cell number ',num2str(i)]);
%     Lcell= (L5 == i);
%       p.maxThreshCut = 0.3;
%       p.maxCellWidth = 20;  
%     cutcell= cutcurvbpts(Lcell, p.maxThreshCut, p.maxCellWidth);
%     
%     cellnos= unique(cutcell);  
%     L5(Lcell)= 0;
%     label= max2(L5);
%     for j= 2:length(cellnos)
%         L5(find(cutcell == cellnos(j)))= label+j;
%     end;
% end;
%L5= renumberimage(L5);
L5 = L3;
% Removing single pixels and small speckles
L6= L5;
r= regionprops(L6, 'area');
fpts= find([r.Area]<p.AreaThresh);
for i= 1:length(fpts),
    L6(L6 == fpts(i))= 0;
end;
L6= renumberimage(L6);

% CUTTING INDIVIDUAL CELLS at the narrow waist: makes sure septation events
% are properly identified and that cells are identified seperately. Cuts 
% cells at points where both sides of the cell are sufficiently concave.
%   maxthresh:      smaller maxtresh the more points are cut
%   maxcellwidth:   two points must be closer than this for cut to be accepted
%   mincelllength:  cut ignored if it creates `cell' smaller than this
%COmmented out 8/2/08
% r= regionprops(L6, 'solidity');
% fkinks= find(([r.Solidity] < 0.90));
% disp(['Cutting individual cells(',num2str(length(fkinks)),').']);
% for i= 1:length(fkinks)
%     %disp([' ',num2str(i),':cutting cell number ',num2str(fkinks(i))]);
%     Lcell= (L6 == fkinks(i));
%     cutcell= cutcurv(Lcell, 0.1, 30, 10);
%     
%     cellnos= unique(cutcell);  
%     L6(Lcell)= 0;
%     label= max2(L6);
%     for j= 2:length(cellnos)
%         L6(find(cutcell == cellnos(j)))= label+j;
%     end;
% end;
L6= renumberimage(L6);

% BREAKING CELLS
% break up big cells
% Goes along the THINned cell and looks for places where there is a change 
% in the phase value - i.e. where there could be a space between cells.
% uses PHSUB(:,:,p.imNumber1) (default p.imNumber1=2)
L7= L6;
r= regionprops(L6, 'majoraxislength');
fbiggies= find(([r.MajorAxisLength]>p.MajorAxisLength));
disp(['Breaking up big cells(', num2str(length(fbiggies)),').']);
for i = 1:length(fbiggies),
    if beadmode
        L7(L6 == fbiggies(i))= 0;
    else
        %disp([' ',num2str(i),': checking cell number ',num2str(fbiggies(i))]);
        Lcell= +(L6 == fbiggies(i)); % + converts logical to double
        Lcell(Lcell == 1)= fbiggies(i);
        cutcell= breakcellfluor(Lcell, imt(:,:,1), ...
                           p.maxThresh, p.minThresh, p.minCellLength);

        L7(L6 == fbiggies(i))= 0;
        % place cutcell
        label= max2(L7);
        for j = 1:max2(cutcell),
            L7(cutcell==j)= label+j;
        end;
    end
end;
L7= renumberimage(L7);


% DEKINKING CELLS with low solidity
% Tries again to look for phase extrema along the thin, using another phase
% slice (if available). Then uses "dekinker" to look for sharp angles along
% the bacterial spine (thin) - if those are found, the cell is cut at the
% corner.
% uses PHSUB(:,:,p.imNumber2) (default p.imNumber2=1)
if beadmode
    L8 = L7;
else
    L8= L7;
    r= regionprops(L7, 'solidity');
    fkinkies= find(([r.Solidity] < 0.85));
    disp(['Breaking kinky cells(',num2str(length(fkinkies)),').']);
    for i = 1:length(fkinkies)
        %disp([' ',num2str(i),': dekinking cell number ',num2str(fkinkies(i))]);
        Lcell= +(L7 == fkinkies(i)); % + converts logical to double
        Lcell(Lcell == 1)= fkinkies(i);

        cutcell= breakcellfluor(Lcell, imt(:,:,1), ...
                       p.maxThresh, p.minThresh, p.minCellLength);
        if max2(cutcell) == 1
            cutcell= dekinker(Lcell, p.radius, p.minCellLength, p.angThresh);
        end;

        L8(L7 == fkinkies(i))= 0;
        % place cutcell
        for j = 1:max2(cutcell),
            L8(cutcell==j)= max2(L8)+j;
        end;
    end;
    L8= renumberimage(L8);
end

% FINAL IMAGE
disp('Almost done.');
r= regionprops(L8, 'area');
flittles= find([r.Area] < p.minCellArea);
for i= 1:length(flittles),
    % delete cell
    L8(L8 == flittles(i))= 0;
end;
L8= renumberimage(L8);
L9=L8;


% 
% % FILL IN EDGE
% % reconstruct image from edge and L3
% disp(['Filling cells from edge(',num2str(max2(L8)),').']);
% L8= zeros(size(L3));
% % mesub = mesub | (be & bigmask);    
% for i= 1:max2(L9)
%     %disp([' filling cell number ',num2str(i)]);
%     fx=[];strelsize=4; % Required to get rid of inter-cell garbage.
%     while length(fx)==0 & strelsize>0
%         strelsize=strelsize-1;
%         L9i = imerode(L9==i,strel('disk',strelsize));
%         [fx, fy]= find(L9i);
%     end
%     tempim1= imfill(mesub, [fx, fy]);
%     % The previous procedure could miss cells if they are paired and
%     % connected by a closed-off edge (i.e. do not form one cavity):
%     % [fx, fy]= find(L3 == i);
%     % pt= round(length(fx)/2);
%     % tempim1= imfill(mesub, [fx(pt), fy(pt)]);
%     tempim2= tempim1 & ~mesub;
%     L9(tempim2 > 0)= i;
% end;
% L9(bigmask==0)=0;
% outskirts=bigmask & ~imdilate(mask, strel('disk',3));
% outlist=unique(L4(outskirts));
% if length(outlist)>1;outlist(outlist==0)=[];end
% JCR: commented out next line as an experiment, don't like what it's doing- 
% it can remove an entire colony
% %for olc=outlist';L4(L4==olc)=0;end
% L9= bwlabel(L9>0,4);
%L4= renumberimage(imclearborder(L4));


if beadmode
    if emptyflag
        L = renumberimage(zeros(ph3size(1:2)));
    else
        L = L8;
    end
else
    %this is great - and not at all ad hoc
    %L9 = imerode(L9,strel('disk',2));
    L9= carefuldilate(+L9, strel('diamond',1), 1);

    % removing small cells again...
    r= regionprops(L9,'majoraxislength');
    flittles= find([r.MajorAxisLength] < p.minCellLengthConservative);
    disp(['Removing small cells(', num2str(length(flittles)),').']);
    for i= 1:length(flittles),
        % delete cell
        L9(L9 == flittles(i))= 0;
    end;
    L9= renumberimage(L9);

    %commented by JCWL to avoid cropping areas in schnitz edit 8/02/08
    % extra=5;
    % LNfull=zeros(size(ph3(:,:,1)));
    % LNfull(rect(1):rect(3), rect(2):rect(4))=L9;
    % [fx,fy]= find(LNfull);
    % xmin= max(min(fx) - extra, 1);
    % xmax= min(max(fx) + extra, size(LNfull,1));
    % ymin= max(min(fy) - extra, 1);
    % ymax= min(max(fy) + extra, size(LNfull,2));
    % rect= [xmin ymin xmax ymax];
    % phsub= ph3(xmin:xmax, ymin:ymax,  :) ;
    % L9=LNfull(xmin:xmax, ymin:ymax);

    L= zeros(size(ph3(:,:,1)));
    L = L9;
    if max2(L)==0 % i.e. if there are no cells here...
        disp('Àoh no, there are no cells on this frame...?');
        %keyboard;
    end;
    disp(['Total cells in this frame: ', num2str(max(max(L))),'.']);
end