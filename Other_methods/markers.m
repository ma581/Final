% This method aims to watershed using both RFP edge detection and the
% composite image x0 (T0 and P0). 

% It runs by:
% 1. Finding edges  using RFP
% 2. Finding edges using Phase
% 3. Finding thick edge using RFP and Phase
% 4. Finding Foreground and Background markers
% 5. Imposing markers on the gradient image
% 6. Watershed

%% 1. Edge using RFP
e = edge(T0,'log',0,0.7);imshow(e);
thresh = im2bw(T0,graythresh(T0));imshow(thresh); thresh = imdilate(thresh,strel('disk',1)); %Used to keep only region of interest
e = e.*thresh;imshow(e);
e = bwmorph(e,'clean');imshow(e);
e = bwareaopen(e,10,8);imshow(e);
% e = bwmorph(e,'thick');imshow(e);

e_out = e;

%% 2. Edge using Phase
e2 = edge(P0,'log',0,0.5);imshow(e2);
threshP = imdilate(thresh,strel('disk',4));imshow(threshP);
e2 = e2.*threshP;imshow(e2);
e2 = bwmorph(e2,'clean');imshow(e2);
e2 = bwareaopen(e2,10,8);imshow(e2);
e2 = bwmorph(e2,'thick');imshow(e2);

e2_out = e2;

%% 3. Thick edge using x0 (T0 and P0)
% x_filt = uint16(wiener2(x0,[50 50])); imshow(x_filt);
% x_adapthist = adapthisteq(x_filt,'NumTiles',[8 8],'Distribution','rayleigh','Alpha',0.9);%imshow(M);[x,counts]=imhist(M);figure,stem(x,counts);
 x_gammacontrast = imadjust(x0,[],[],0.2); imshow(x_gammacontrast);
 x_close = imclose(x_gammacontrast,strel('disk',3));imshow(x_close);
 x_thresh = im2bw(x_gammacontrast,graythresh(x_gammacontrast));imshow(x_thresh);

% x_out = double(x0);
% x_out = x_out./maxmax(x_out);imshow(x_out*5);
% x_out = bwmorph(x_out,'thick');imshow(x_out);

x_out = x_thresh;

% surf(double(x_out),'EdgeColor','none');colormap('jet');

%% Hybrid image for sharper edges

% % Hybrid edge using x0 and T0
% edge_weight = 0.5;
% % h = (1-edge_weight)*x_out + edge_weight*double(e_out);imshow(h);
% h1 = x_out.*e_out;imshow(h1);
% % surf(double(h1),'EdgeColor','none');colormap('jet');
% % h_thresh1 = im2bw(h,graythresh(h1));imshow(h_thresh1);
% % h_close = imclose(h_thresh,strel('disk',1));imshow(h_close);
% % h_close1 = bwmorph(h1,'thick');imshow(h_close1);
% h1 = bwmorph(h1,'clean');imshow(h1);
% h1 = bwareaopen(h1,10,8);imshow(h1);
% % Hybrid edge using P0 and T0
% h1 = bwmorph(h1,'skel',Inf);
% h_out = h1;
% 
% % h_close = imclose(h1,strel('disk',1));imshow(h_close);



%% 4. Background and Foreground


A = edge(T0,'log',0,0.9);imshow(A);
thresh = im2bw(T0,graythresh(T0));imshow(thresh); thresh = imdilate(thresh,strel('disk',1));

A = A.*thresh;imshow(A);
A = imfill(A,'holes');imshow(A);
A = bwareaopen(A, 20,4);imshow(A);
A = imcomplement(A);imshow(A);
A = bwmorph(A,'remove');imshow(A);

bgm = A.*x_out;imshow(bgm);
bgm =  bwareaopen(bgm, 10,8);imshow(bgm);
bgm = bwmorph(bgm,'clean');imshow(bgm);

fgm = imcomplement(bgm);imshow(fgm);
fgm = imerode(fgm,strel('disk',5));imshow(fgm);
fgm = fgm.*thresh;imshow(fgm);
fgm = bwareaopen(fgm, 20,8);imshow(fgm);

%% 5. Using Marker watershed
% Using the hybrid edge image as the bgm
% Trying to find better fgm

% Gradient
I = T0;   
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

% Watershed
gradmag2 = imimposemin(gradmag, fgm | bgm);
L = watershed(gradmag2);
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');imshow(Lrgb);
imshowlabel(bwlabel(L));
i











% rgb1 = label2rgb(bwlabel(L),'jet', 'w', 'shuffle'); imshow(rgb1);
% rgb2 = label2rgb(bwlabel(L),'jet', 'w', 'shuffle'); imshow(rgb2);
% rgb3 = label2rgb(bwlabel(L),'jet', 'w', 'shuffle'); imshow(rgb3);
% rgb4 = label2rgb(bwlabel(L),'jet', 'w', 'shuffle'); imshow(rgb4);
% 
% figure,subplot(2,2,1);imshow(rgb1);title('Movie 1 frame 140');
% subplot(2,2,2);imshow(rgb2);title('Movie 2 frame 140');
% subplot(2,2,3);imshow(rgb3);title('Movie 3 frame 161');
% subplot(2,2,4);imshow(rgb4);title('Movie 4 frame 161');

 %% Randoms
% 
% % Background
% 
% bw =  imadjust(P0,stretchlim(P0),[0 1]).*uint16(threshP);imshow(bw);
% bw = imclose(bw,strel('disk',5)); imshow(bw);
% bw = imregionalmax(bw,4);imshow(bw);
% bw = im2bw(bw,graythresh(bw));imshow(bw);
% 
% 
% 
% bw = imadjust(x0,stretchlim(x0),[0 1]).*uint16(thresh);imshow(bw);
% bw = imclose(bw,strel('disk',6)); imshow(bw);
% bw = im2bw(bw,graythresh(bw)*.5);imshow(bw);
% bw = bwmorph(bw,'close');imshow(bw);
% % bw = imerode(bw,strel('disk',1));imshow(bw);
% % bw = bwperim(imclose(edge(bw,'log',0,2).*logical(thresh),strel('disk',4)));imshow(bw);
% % bw = bwmorph(bw,'skel',Inf);imshow(bw);
% % bw = imopen(bw,strel('disk',1));imshow(bw);
% bw  =(bwmorph(bw,'thin'));imshow(bw);
% 
% 
% h_out = imdilate(h_close1,strel('disk',2));imshow(h_out);
% h_out = bwmorph(h_out,'thick');imshow(h_out);
% imshowpair(h_out,e_out);
% edgeimage = h_out.*e_out;imshow(edgeimage);
% imshow(imfill(edgeimage,'holes'));

 %% Foreground
% % x_gammacontrast = imadjust(x_adapthist,[],[],1); imshow(x_gammacontrast);
% % fgm = im2bw(imcomplement(x_gammacontrast),graythresh(imcomplement(x_gammacontrast))*1.02);imshow(fgm);
% % fgm = imerode(fgm,strel('disk',5));imshow(fgm);
% % fgm = fgm.*thresh;imshow(fgm);
% 
% A = edge(T0,'log',0,0.9);imshow(A);
% thresh = im2bw(T0,graythresh(T0));imshow(thresh); thresh = imdilate(thresh,strel('disk',1));
% 
% A = A.*thresh;imshow(A);
% A = imfill(A,'holes');imshow(A);
% A = bwareaopen(A, 20,4);imshow(A);
% 
% % r = regionprops(A,'Centroid');
% % A1 = zeros(size(A));
% % 
% %     for j=1:length([r.Centroid])/2
% %         x = round(r(j).Centroid(:,1));
% %         y = round(r(j).Centroid(:,2));
% %         A1(y,x) = 1;
% %     end
% 
% 
%  A = imerode(A,strel('disk',1));imshow(A);
% A = bwmorph(A,'thin');imshow(A);
% A = bwmorph(A,'clean');imshow(A);
% A = bwmorph(A,'shrink');imshow(A);
% A = bwmorph(A,'shrink');imshow(A);
% A = bwmorph(A,'clean');imshow(A);
