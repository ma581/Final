% The GAMMA CORRECTION AND THRESHOLD method
% This is currently written to run from my laptop, accessing the four
% movies Movie 1,2,3 and 4 from my pendrive. 
%
% Will need to edit the script to access files from another PC. 
% To run, call the script GAMMAMETHOD

% This method :
% 1. Wiener Filters
% 2. Gamma correction
% 3. Thresholding
% 4. Morphology and splitting



Laptop = input('Laptop? [Yes = 1 No = 0] ');

% Bruno = input('Bruno? ');
set = input('Which Image set? ');

if set == 1
    segRange = 1:147;
    
elseif set == 2
    segRange = 1:170;
    
elseif set == 3
    segRange = 1:170;
elseif set == 4
    segRange = 1:198;
end

% segRange = input('What is segRange? ');
Bruno = 1;

p.minCellArea = 100;
p.maxThresh = 0.050;
p.minThresh = 0.050 ;
p.minCellLength = 20;


%%
for i =segRange
    
    n = 1000+i;
    ns = num2str(n);
    ns(1) = [];
    
    % % On Manoj's Laptop [Bruno files]
    %     P0 = imread(['C:\Users\Manoj\Dropbox\SainsburyLab\RawImages\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-11_1-p-',ns,'.TIF']);
    %     T0 = imread(['C:\Users\Manoj\Dropbox\SainsburyLab\RawImages\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-11_1-t-',ns,'.TIF']);
    %
    if Bruno == 1
        
        if set == 1
            
            if Laptop == 1
                P0 = imread(['E:\1\SubFramedMovies\PSB7-11_1-p-',ns,'.TIF']);
                T0 = imread(['E:\1\SubFramedMovies\PSB7-11_1-t-',ns,'.TIF']);
                
            else
                %         % On Chris' PC  [Bruno files]
                P0 = imread(['C:\Movies\SubFramedMovies\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-11_1-p-',ns,'.TIF']);
                T0 = imread(['C:\Movies\SubFramedMovies\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-11_1-t-',ns,'.TIF']);
            end
            
        elseif set == 2
            
            if Laptop == 1
                P0 = imread(['E:\2\2015-03-11-Bruno-Cyano2\PSB7-16_1-p-',ns,'.TIF']);
                T0 = imread(['E:\2\2015-03-11-Bruno-Cyano2\PSB7-16_1-t-',ns,'.TIF']);
            else
                
                %             % On Chris' PC  [Bruno files]
                P0 = imread(['Z:\Manoj\AnalysedMoviesManoj\RawImages_SubFramed\2015-03-11-Bruno-Cyano2\PSB7-16_1-p-',ns,'.TIF']);
                T0 = imread(['Z:\Manoj\AnalysedMoviesManoj\RawImages_SubFramed\2015-03-11-Bruno-Cyano2\PSB7-16_1-t-',ns,'.TIF']);
            end
            
        elseif set == 3
            
            if Laptop == 1
                P0 = imread(['E:\3\2015-03-18-Bruno-Cyano2\PSB-01_1-p-',ns,'.TIF']);
                T0 = imread(['E:\3\2015-03-18-Bruno-Cyano2\PSB-01_1-t-',ns,'.TIF']);
            else
                
                % On Chris' PC  [Bruno files]
                P0 = imread(['Z:\Manoj\AnalysedMoviesManoj\RawImages_SubFramed\2015-03-18-Bruno-Cyano2\PSB-01_1-p-',ns,'.TIF']);
                T0 = imread(['Z:\Manoj\AnalysedMoviesManoj\RawImages_SubFramed\2015-03-18-Bruno-Cyano2\PSB-01_1-t-',ns,'.TIF']);
            end
            
        elseif set == 4
            
            if Laptop == 1
                P0 = imread(['E:\4\PSB7-13_1-p-',ns,'.TIF']);
                T0 = imread(['E:\4\PSB7-13_1-t-',ns,'.TIF']);
            else
                
                P0 = imread(['Y:\Bruno\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-13_1-p-',ns,'.TIF']);
                T0 = imread(['Y:\Bruno\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-13_1-t-',ns,'.TIF']);
            end
        end
        
    else
        % On Chris' PC [Chris files]
        P0 = imread(['C:\Movies\SubFramedMovies\2014_03_24\bacillus-01_1-p-',ns,'.TIF']);
        T0 = imread(['C:\Movies\SubFramedMovies\2014_03_24\bacillus-01_1-t-',ns,'.TIF']);
    end
    
    
    % Previous technique
    p0 = double(P0);
    t0 = double(T0);
    
    % New technique
    %    p0 = p0*medallp/median(p0(:));
    P0 = uint16(p0); %converts to integer
    
    %   t0 = t0*medallp/median(t0(:));
    T0 = uint16(t0); %converts to integer
    
    % Combining images
    phase = double(P0);phase = (phase - minmin(phase))./(maxmax(phase)-minmin(phase));
    tritc = double(T0);tritc = (tritc - minmin(tritc))./(maxmax(tritc)-minmin(tritc));
    combined = ((1-phase).*tritc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Manoj Edit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Wiener filter on RFP
    x_filt = (wiener2((1-phase),[3 3])); imshow(x_filt);
    thresh = im2bw(tritc,graythresh(tritc));imshow(thresh); thresh = imdilate(thresh,strel('disk',4)); %Used to keep only region of interest
    
    %Adaptive histogram correctiond
    %     x_adapthist = adapthisteq(x_filt,'NumTiles',[8 8],'Distribution','exponential','Alpha',0.9);%imshow(x_adapthist);
    % [x,counts]=imhist(x_adapthist);stem(x,counts);
    
    
  %% Image 
    x = (1-phase);imshow(x,[]);
    
    
    %% Gamma contrast
    x_gammacontrast = imadjust(x,[0.6 1],[0 1],5); imshow(x_gammacontrast);
%     x_gammacontrast = x;
    
    %% Threshold
    e = im2bw(x_gammacontrast,graythresh(x_gammacontrast));imshow(e);
    e = e.*thresh;imshow(e);
    
    %     f = imfill(e,'holes');imshow(f);
%     f = bwmorph(e,'clean'); imshow(f);
        f = bwareaopen(e,10,4);imshow(f);
        
        
  e = e+bwmorph(e,'majority');imshow(e);
    e = bwmorph(e,'hbreak');imshow(e);

    
    f = imfill(e,8,'holes');imshow(f);
    f = bwmorph(f,'clean'); imshow(f);
    f = bwmorph(f,'hbreak');imshow(f);
    
    f = f+bwmorph(f,'majority');imshow(f);
    %         f = bwmorph(f,'hbreak');imshow(f);
    
    f = bwareaopen(f,80,4);imshow(f);
    
    f = bwmorph(f,'thin',3);imshow(f);
    f = imfill(f,8,'holes');imshow(f);
    f = bwmorph(f,'thick',3);imshow(f);
   
    x_out = f;imshow(x_out);
    
    
    %% Preliminary result
    bw = bwlabel(x_out);
    rgb2 = label2rgb(bw,'jet',[.5 .5 .5]); imshow(rgb2);
    %     rgb3 = label2rgb(x_out,'jet', 'w', 'shuffle');imshow(rgb3);
%     imwrite(rgb2,['Prel_thresh',num2str(i),'.jpg']);
    i
    
    %% A little Schnitz!
    h = x_out;
%     h = (bwmorph(x_out,'open',Inf));imshow(label2rgb(bwlabel(h),'jet',[.5 .5 .5]));
%     h = (bwmorph(h,'erode',2));imshow(label2rgb(bwlabel(h),'jet',[.5 .5 .5]));
    h = imopen(h,strel('disk',2));;imshow(label2rgb(bwlabel(h),'jet',[.5 .5 .5]));
    h = bwmorph(h,'hbreak');imshow(label2rgb(bwlabel(h),'jet',[.5 .5 .5]));
%         h = (bwmorph(h,'dilate',2));imshow(label2rgb(bwlabel(h),'jet',[.5 .5 .5]));
    
        h = imfill(h,8,'holes');imshow(h);
    %      h = (~bwmorph(~h,'diag'));imshowlabel(bwlabel(h));
    bw = bwlabel(h);
    rgb3 = label2rgb(bw,'jet',[.5 .5 .5]); imshow(rgb3);
%     imwrite(rgb3,['Schnitz_thresh',num2str(i),'.jpg']);
     
    L4 = bw;
    
    %% CUTTING INDIVIDUAL CELLS at the narrow waist: makes sure septation events
    % are properly identified and that cells are identified seperately. Cuts
    % cells at points where both sides of the cell are sufficiently concave.
    %   maxthresh:      smaller maxtresh the more points are cut
    %   maxcellwidth:   two points must be closer than this for cut to be accepted
    %   mincelllength:  cut ignored if it creates `cell' smaller than this
    %COmmented out 8/2/08
    L6 = L4;
    r= regionprops(L6, 'solidity');
    fkinks1= find(([r.Solidity] > 0.75 ));
    fkinks2= find(([r.Solidity] < 0.90 ));
    fkinks=intersect(fkinks1,fkinks2);
    disp(['Cutting individual cells(',num2str(length(fkinks)),').']);
    for z= 1:length(fkinks)
        %disp([' ',num2str(i),':cutting cell number ',num2str(fkinks(i))]);
        Lcell= (L6 == fkinks(z));
        cutcell= cutcurv(Lcell, 0.4, 30, 10);
        
        cellnos= unique(cutcell);
        L6(Lcell)= 0;
        label= max2(L6);
        for j= 2:length(cellnos)
            L6(find(cutcell == cellnos(j)))= label+j;
        end;
    end;
    L6= renumberimage(L6);
    imshowlabel(L6);
    
    %% BREAKING CELLS
    % break up big cells
    % Goes along the THINned cell and looks for places where there is a change
    % in the phase value - i.e. where there could be a space between cells.
    % uses PHSUB(:,:,p.imNumber1) (default p.imNumber1=2)
    L7= L6;
    
    
    p.maxThresh = 0.050;
    p.minThresh = 0.050 ; %0.050   % changed from 0.05 5/2/08
    p.minCellLength = 30;
    
    
    r= regionprops(L6, 'majoraxislength');
    fbiggies= find(([r.MajorAxisLength]>50));
    disp(['Breaking up big cells(', num2str(length(fbiggies)),').']);
    imshowlabel(L7);
    for z = 1:length(fbiggies),
        %disp([' ',num2str(i),': checking cell number ',num2str(fbiggies(i))]);
        Lcell= +(L6 == fbiggies(z)); % + converts logical to double
        Lcell(Lcell == 1)= fbiggies(z);
        cutcell= breakcellfluor(Lcell, imcomplement(combined), ...
            0.2, 0.005, p.minCellLength);
        
        L7(L6 == fbiggies(z))= 0;
        % place cutcell
        label= max2(L7);
        for j = 1:max2(cutcell),
            L7(cutcell==j)= label+j;
        end;
    end;
    L7= renumberimage(L7);
    imshowlabel(L7);

    
    
end

