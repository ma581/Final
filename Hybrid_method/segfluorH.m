function [phsub, L9, rect]= segfluorH(ph3,PhaseIm, p)



% Modified version of SEGFLUORC to incorporate the hybrid method for segmentation.
% Example :
%
% L = segfluorC(Im_phase,Im_RFP,p);
%
% 1. Filter (segfluorC) - to output image col_center for Mask
% 2. Mask (segfluorC) - to output rect
% 3. Hybrid method of segmentation - calls HYBRIDSEG(PhaseIm,ph3) to output L9 the labeled cells
% - Manoj - 








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

if p.imNumber2>size(ph3,3)
    p.imNumber2 = 1;
end
if p.imNumber1>size(ph3,3)
    p.imNumber1 = size(ph3,3);
end

maxsizeofdirt=20000;

%% FILTER AND EDGE PHASE
% median filter and rescale phase images
for i= 1:size(ph3,3),
    ph3(:,:,i)= medfilt2(ph3(:,:,i),[3 3]);
    
    x= double(ph3(:,:,i));
    s= sort(x(:));
    small= s(25);
%     big= s(end-25);
signal = s-(small+40);
signal(signal<=0) = [];
signal = signal+(small+40);
big = signal(round(length(signal)*0.93));
    rescaled=(x - small)/(big - small);
    rescaled(rescaled<0)= 0;
    ph3(:,:,i)= uint16(10000*rescaled);
end;
me= edge(ph3(:,:,p.segmentationPhaseSlice),'log',0,p.edge_lapofgauss_sigma);
col_center=[size(me,1),size(me,2)]/2;

%% FIND MASK (smallest region containing cells)

% MASK: first iteration
%   Run box of size boxsize over the edge image 
%   and count the number of white pixels in the box.
%   Cells have less than minonpixels pixels.

imout= zeros(size(ph3(:,:,1)));
replaceboxsize= 5;  % boxsize is shrunk to replaceboxsize on output image 
boxsize= 40;

% was minnopixels, now p.minNumEdgePixels in schnitzcell parameter structure
%minNumEdgePixels= 300; % original for e.coli
minNumEdgePixels= 150; % note change (from 300 to 250 during bacillus tests)
% JCR: significance of 300 is that's the value of minNumEdgePixels for e.coli
% JCR: bacillus is smaller, until 2005-08-18 it was 250 (e.g. for methods 
%      paper work); since then it's 215

for x= 1:replaceboxsize:size(me,1) - boxsize
    for y= 1:replaceboxsize:size(me,2) - boxsize
        xmin= x + boxsize/2;
        xmax= x + boxsize/2 + replaceboxsize - 1;
        ymin= y + boxsize/2;
        ymax= y + boxsize/2 + replaceboxsize - 1;
        
        subregion= me(x:x + boxsize, y:y + boxsize);
        checkiscell= sum(subregion(:)) < p.minNumEdgePixels;
        imout(xmin:xmax, ymin:ymax)= ones(replaceboxsize)*checkiscell;
    end
end
% dilate output image to make sure ends of cells included
imout = uint8(imdilate(imout,strel('disk',replaceboxsize)));

% tidy up imout so that only one highlighted region remains
% imout = imopen(imout,strel('disk',25)); % new
% this waxs replaced by the following 5 lines, to make sure something remains!
imout_temp = imopen(imout,strel('disk',25)); % new
if max(imout_temp(:))
    imout = imout_temp;
end
clear imout_temp;

if p.minNumEdgePixels<300
    imout = imdilate(imout,strel('disk',10)); % Bruno (seems logical that this is what was intended...)
end
imoutlabel= bwlabel(imout);
r= regionprops(imoutlabel, 'area', 'centroid');

    % cells are in small groups
    % keep the clump closest to the centre of the image
    cen= [r.Centroid];
    xcen= cen(1:2:end);
    ycen= cen(2:2:end);
    
    % find distance from centre of image
    distances= sqrt((xcen-col_center(2)).^2 + (ycen-col_center(1)).^2);
    [ds, di]= sort(distances);
    % keep closest highlighted region only
    imout2= zeros(size(imout));
    imout2(imoutlabel == di(1))= 1;
    
    % keep additional cell clumps if the number of colonies is greater than one (Nitzan)
    for i= 2:length(di);
        imout2(imoutlabel == di(i))= 1;
    end;



rect = [1 1 size(ph3,1) size(ph3,2)];
phsub = ph3(:,:,1);
imout2sub= imout2;
mesub=me;



%% HYBRID METHOD
L7 = hybridseg(PhaseIm,ph3);

L9 = L7;
rgb3 = label2rgb(L9,'hsv','k','shuffle'); imshow(rgb3);






L= zeros(size(ph3(:,:,1)));
L = L9;
if max2(L)==0 % i.e. if there are no cells here...
    disp('Àoh no, there are no cells on this frame...?');
    %keyboard;
end;
disp(['Total cells in this frame: ', num2str(max(max(L))),'.']);