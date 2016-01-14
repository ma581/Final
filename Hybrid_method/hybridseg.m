function segmented_cells = hybridseg(PhaseIm,RFPIm)

% Carries out the hybrid method of segmentation developed. 
% Example :
%
% L = hybridseg(PhaseIm,RFPIm);
%
% Where L is the labelled region of cells
%       PhaseIm = The phase contrast image
%       RFPIm   = The RFP channel image
%
% Hybrid method of segmentation carries out :
%
% 1. Normalising and combining images
% 2. Gamma contrast correction on combined image
% 3. Phase image gamma correction and adaptive threshold
% 4. Phase image breaking big cells
% 5. Multiple edge detection and fill
% 6. Schnitzcells edge detection
% 7. Combining multiple edges+phase thresold + fill
% 8. Morphology for splitting cells
% 9. Cutting cells at narrow waist
% 10. Breaking big cells
% - Manoj - 



%% 1. Normalising and combining images
phase = double(PhaseIm); phase = (phase - minmin(phase))./(maxmax(phase)-minmin(phase));
tritc = double(RFPIm); tritc = (tritc - minmin(tritc))./(maxmax(tritc)-minmin(tritc));
combined = ((1-phase).*tritc);

% Boundary
thresh = im2bw(tritc,graythresh(tritc));%imshow(thresh); 
thresh = imdilate(thresh,strel('disk',1)); %Used to keep only region of interest
x = combined;%imshow(x,[]);

%% 2. Gamma contrast correction on combined image
x_gammacontrast = imadjust(x,[0 1],[0 1],1.4); %imshow(x_gammacontrast);
% IMADJUST(I,[LOW_IN; HIGH_IN],[LOW_OUT; HIGH_OUT],GAMMA)
% Adjusting LOW_IN enables contrast adjustment of only suitable pixels if
% you identify a bimodal distribution in the image histogram

% GAMMA > 1 : mapping is weighted toward lower (darker) output values
% GAMMA < 1 : mapping is weighted toward igher (brighter) output values
% GAMMA = 1 : mapping is linear


%% 3. Phase image gamma correction and adaptive threshold
% Wiener filter on Phase
x_filt = (wiener2((1-phase),[3 3])); %imshow(x_filt); 
% [3 3] region is the standard size 

% Gamma contrast
phase_gammacontrast = imadjust(x_filt,[0 1],[0 1],1.1); %imshow(phase_gammacontrast); %previously from 0.6

% Adaptive Threshold
thresholdimage = adaptivethreshold(phase_gammacontrast,20,0).*thresh;%imshow(thresholdimage);
% ADAPTIVETHRESHOLD(I,WINDOW_SIZE,THRESHOLD) Changing the windowsize
% affects the threshold quality the most


% Filling the holes a little
thresholdimage = ~bwareaopen(~thresholdimage,10,4);%imshow(thresholdimage);
thresholdimage = ~bwareaopen(~thresholdimage,30,8);%imshow(thresholdimage);
thresholdimage = imopen(thresholdimage,strel('disk',2));%imshow(thresholdimage);

%% 4.  Phase image breaking big cells
% % break up big cells
% % Goes along the THINned cell and looks for places where there is a change
% % in the phase value - i.e. where there could be a space between cells.

L7 = breakingcells(thresholdimage, phase,0.02,0.004,30);
% Output = breakingcells(I,PHASE,maxthresh, minthresh, minCellLength)
% maxthresh & minthresh:    threshold values (smaller they are the more cuts are included)
% mincelllength:    cuts which create cells smaller than this will be ignored

% thresholdimage = renumberimage(L7);

thresholdimage = (L7);
%imshowlabel(thresholdimage);

% Little erosion
thresholdimage = imerode(thresholdimage,strel('disk',1)); %imshowlabel(thresholdimage);
thresholdimage = im2bw(thresholdimage);    %imshowlabel(thresholdimage);


%% 5. Multiple edge detection and fill
e = imtophat(x_gammacontrast,strel('disk',100));%imshow(e); % Some cleaning and smoothing
e = edge(e,'log',0,0.9).*thresh;%imshow(e); % Edge detection
e = bwmorph(e,'clean'); %imshow(e); % Some cleaning of noisy pixels
f1 = imfill(e,'holes');%imshow(f1);% Filling cells

% 6. Schnitzcells edge detection
e1 = imtophat(RFPIm,strel('disk',100));%imshow(e1); % Some cleaning and smoothing
e1 = edge(e1,'log',0,0.9).*thresh;%imshow(e1); % Edge detection
e1 = bwmorph(e1,'clean'); %imshow(e1); % Some cleaning of noisy pixels

% Edge before gamma correction
e111 = imtophat(x,strel('disk',100));%imshow(e111); % Some cleaning and smoothing
e111 = edge(e111,'log',0,0.9).*thresh;%imshow(e111);
e111 = bwmorph(e111,'clean'); %imshow(e111);

% Different edge levels to fill the shape
e2 = imopen(x_gammacontrast,strel('disk',2));%imshow(e);

e2 = edge(e2,'log',0,0.01).*thresh+...
    edge(e2,'log',0,0.1).*thresh+...
    edge(e2,'log',0,0.2).*thresh+...
    edge(e2,'log',0,0.3).*thresh+...
    edge(e2,'log',0,0.001).*thresh+...
    edge(e2,'log',0,0.05).*thresh;
e2 = bwmorph(e2,'clean'); %imshow(e2);

% Combined edges
e3 = e+e1+e111+e2; %imshow(e3);
e4 = e3;

% Cleaning edges
e4 = bwmorph(e4,'fill');%imshow(e4);
e5 = bwmorph(e4,'hbreak');%imshow(e5);


%% 7. Combining multiple edges (e5) + phase thresold (thresholdimage) + fill (f1)

f = f1 + e5 + thresholdimage; %imshow(f); % Adding all the stages

% Clean up and fill up
f = ~bwareaopen(~f,20,4);%imshow(f); 
f = ~bwareaopen(~f,22,8);%imshow(f);
f = bwareaopen(f,80,4);%imshow(f);

f = (imerode(f,strel('disk',1)));%imshow(f); % Little erosion

x_out = f;


%% Preliminary result

bw = bwlabel(x_out,8);
L = bwlabel(x_out,8); %For output of this function...
rgb2 = label2rgb(bw,'hsv','k','shuffle'); %imshow(rgb2);
%     imwrite(rgb2,['Preliminary',num2str(i),'.jpg']);

%% 8. Morphology for splitting

h = ~x_out;
for m = 1:5
    h = (~bwmorph(h,'diag'));       %imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle')); % Opening connections
end
h = bwareaopen(h,10,4);             %imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle')); % Cleaning noisy pixels
h = ~bwmorph(~h,'bridge');          %imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle')); % Opening connections
h = bwmorph(h,'hbreak');            %imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle')); % Opening connections
h = imopen(h,strel('disk',4));      %imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle')); % Opening
bw = bwlabel(h,4);
L4 = bw;




 %% 9. CUTTING INDIVIDUAL CELLS at the narrow waist: makes sure septation events
% % are properly identified and that cells are identified seperately. Cuts
% % cells at points where both sides of the cell are sufficiently concave.
% %   maxthresh:      smaller maxtresh the more points are cut
% %   maxcellwidth:   two points must be closer than this for cut to be accepted
% %   mincelllength:  cut ignored if it creates `cell' smaller than this

% L6= renumberimage(L6);

L6 = cuttingcells(L4, 0.4, 30, 10);
% Output  = cuttingcells(I, maxthresh, maxcellwidth, mincelllength)
%   maxthresh:      smaller maxtresh the more points are cut
%   maxcellwidth:   two points must be closer than this for cut to be accepted
%   mincelllength:  cut ignored if it creates `cell' smaller than this



 %% 10. BREAKING CELLS
% % break up big cells
% % Goes along the THINned cell and looks for places where there is a change
% % in the phase value - i.e. where there could be a space between cells.
% % uses PHSUB(:,:,p.imNumber1) (default p.imNumber1=2)
L7= L6;
L7 = breakingcells(L7, phase,0.2,0.04,30);
% Output = breakingcells(I,PHASE,maxthresh, minthresh, minCellLength)
% maxthresh & minthresh:    threshold values (smaller they are the more cuts are included)
% mincelllength:    cuts which create cells smaller than this will be ignored



%% Filling each cell
% This can be used to go cell by cell and close tiny holes or fill up
% incomplete cell shapes

%     L8 = zeros(size(L7));
%     for k = 1:max(max(L7))
%         L2 = zeros(size(L7));
%         L3 = zeros(size(L7));
%         L2(L7==k) = L7(L7==k);
% %         L22 = imclose(L2,strel('disk',1));imshow(L22);
%         L22 = L2;
%         %         L22 = imerode(L22,strel('disk',1));imshow(L22);
% %         L22 = imdilate(L22,strel('disk',1));imshow(L22);
%         %         L22 = imdilate(L2,strel('disk',10));imshow(L22);
%         % %         L33 = imerode(L22,strel('disk',10));imshow(L33);
%         %         L55 = bwmorph(L22,'erode',8);imshow(L55);
%         %         L55 = L55.*k;
%         %         L4(L55==k) = L55(L55==k);imshow(L4);
%         L8(L22==k) = L22(L22==k);imshow(L8);
%     end
% %     L8 = imfill(L8,'holes');imshow(L8);


%% The output

segmented_cells = L7;
rgb3 = label2rgb(segmented_cells,'hsv','k','shuffle'); %imshow(rgb3);


end