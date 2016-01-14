%
% Laptop = input('Laptop? [Yes = 1 No = 0] ');

% Bruno = input('Bruno? ');
% set = input('Which Image set? ');

if set == 1
    segRange = 1:147;
    
elseif set == 2
    segRange = 1:170;
    
elseif set == 3
    segRange = 1:170;
elseif set == 4
    segRange = 1:198;
elseif set == 5
    segRange = 1:184;
end

Bruno = 1;

p.minCellArea = 100;
p.maxThresh = 0.050;
p.minThresh = 0.050 ;
p.minCellLength = 20;


%%
for i = 60:170%segRange
    
    n = 1000+i;
    ns = num2str(n);
    ns(1) = [];
    
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
                
                P0 = imread(['Z:\Bruno\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-13_1-p-',ns,'.TIF']);
                T0 = imread(['Z:\Bruno\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-13_1-t-',ns,'.TIF']);
            end
            
        elseif set == 5
            
            if Laptop == 1
                P0 = imread(['E:\Arijit\SubFramedMovies\Sigma-02_1-p-',ns,'.TIF']);
                T0 = imread(['E:\Arijit\SubFramedMovies\Sigma-02_1-t-',ns,'.TIF']);
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
    %     x_filt = (wiener2(combined,[3 3])); imshow(x_filt);
    thresh = im2bw(tritc,graythresh(tritc));imshow(thresh); thresh = imdilate(thresh,strel('disk',1)); %Used to keep only region of interest
    
    
    %% Image Normalising
    x = combined;imshow(x,[]);
    
    
    %% Gamma contrast on Combined image
    x_gammacontrast = imadjust(x,[0 1],[0 1],1.4); imshow(x_gammacontrast);
    
    
    %% Phase Image
    % Wiener filter on Phase
    x_filt = (wiener2((1-phase),[3 3])); imshow(x_filt);
    
    % Gamma contrast
    phase_gammacontrast = imadjust(x_filt,[0 1],[0 1],1.1); imshow(phase_gammacontrast); %previously from 0.6
    
    %Threshold
    %     thresholdimage = im2bw(phase_gammacontrast,graythresh(phase_gammacontrast));imshow(thresholdimage);
    %     thresholdimage = thresholdimage.*thresh;imshow(thresholdimage);
    %
    thresholdimage = adaptivethreshold(phase_gammacontrast,20,0.0000).*thresh;imshow(thresholdimage);
    
    % Filling the holes a little
    thresholdimage = ~bwareaopen(~thresholdimage,10,4);imshow(thresholdimage);
    thresholdimage = ~bwareaopen(~thresholdimage,30,8);imshow(thresholdimage);
    
    
    thresholdimage = imopen(thresholdimage,strel('disk',2));imshow(thresholdimage);
    
    
    
    % BREAKING CELLS
    % break up big cells
    % Goes along the THINned cell and looks for places where there is a change
    % in the phase value - i.e. where there could be a space between cells.
    % uses PHSUB(:,:,p.imNumber1) (default p.imNumber1=2)
    L6 = bwlabel(thresholdimage);
    L7= bwlabel(thresholdimage);
    
    
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
        %         cutcell= breakcellfluor(Lcell, imcomplement(combined), ...
        %             0.2, 0.005, p.minCellLength);
        cutcell= breakcellfluor(Lcell, imcomplement(phase), ...
            0.02, 0.004, p.minCellLength);
        
        L7(L6 == fbiggies(z))= 0;
        % place cutcell
        label= max2(L7);
        for j = 1:max2(cutcell),
            L7(cutcell==j)= label+j;
        end;
    end;
    thresholdimage= renumberimage(L7);
    imshowlabel(thresholdimage);
    
    
    
    
    
    
    thresholdimage = imerode(thresholdimage,strel('disk',1)); imshowlabel(thresholdimage);
    thresholdimage = im2bw(thresholdimage);    imshowlabel(thresholdimage);
    
    
    %% Edge and fill
    e = imtophat(x_gammacontrast,strel('disk',100));imshow(e);
    e = edge(e,'log',0,0.9).*thresh;imshow(e);
    e = bwmorph(e,'clean'); imshow(e);
    f1 = imfill(e,'holes');imshow(f1);
    
    % Schnitzcells edge
    e1 = imtophat(T0,strel('disk',100));imshow(e1);
    e1 = edge(e1,'log',0,0.9).*thresh;imshow(e1);
    e1 = bwmorph(e1,'clean'); imshow(e1);
    
    % Edge before gamma correction
    e111 = imtophat(x,strel('disk',100));imshow(e111);
    e111 = edge(e111,'log',0,0.9).*thresh;imshow(e111);
    e111 = bwmorph(e111,'clean'); imshow(e111);
    
    % Different edge levels to fill the shape
    e2 = imopen(x_gammacontrast,strel('disk',2));imshow(e);
    
    e2 = edge(e2,'log',0,0.01).*thresh+...
        edge(e2,'log',0,0.1).*thresh+...
        edge(e2,'log',0,0.2).*thresh+...
        edge(e2,'log',0,0.3).*thresh+...
        edge(e2,'log',0,0.001).*thresh+...
        edge(e2,'log',0,0.05).*thresh;
    e2 = bwmorph(e2,'clean'); imshow(e2);
    
    % Combined edges
    e3 = e+e1+e111+e2; imshow(e3);
    e4 = e3;
    % Cleaning edges
    e4 = bwmorph(e4,'fill');imshow(e4);
    e5 = bwmorph(e4,'hbreak');imshow(e5);
    
    
    % Superposition of imfill, edges and phase threshold
    f = f1 + e5 + thresholdimage; imshow(f);
    %     f = adaptivethreshold(5*(imerode(1-phase,strel('disk',3)))+f1+e5,80,0).*thresh;imshow(f);
    
    %         f = adaptivethreshold(1.5*(thresholdimage)+f1+e5,80,0).*thresh;imshow(f);
    
    % Clean up and fill up
    f = ~bwareaopen(~f,20,4);imshow(f);
    f = ~bwareaopen(~f,22,8);imshow(f);
    f = bwareaopen(f,80,4);imshow(f);
    
    f = (imerode(f,strel('disk',1)));imshow(f);
    %             f = ~bwareaopen(~f,200,8);imshow(f);
    
    %     f = imfill(f,8,'holes');imshow(f);
    x_out = f;
    
    
    
    %% Preliminary result
    bw = bwlabel(x_out,8);
    rgb2 = label2rgb(bw,'hsv','k','shuffle'); imshow(rgb2);
    %     imwrite(rgb2,['Preliminary',num2str(i),'.jpg']);
    i
    
    %% A little Schnitz!
    h = ~x_out;
    for m = 1:5
        h = (~bwmorph(h,'diag'));imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    end
    
    %     h = bwmorph(h,'majority');imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    h = bwareaopen(h,10,4);imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    
    h = ~bwmorph(~h,'bridge');imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    %     h = imfill(h,8,'holes');imshow(h);
    
    h = bwmorph(h,'hbreak');imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    %      h = (~bwmorph(~h,'diag'));imshowlabel(bwlabel(h));
    h = imopen(h,strel('disk',4));imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    %     h = bwmorph(h,'open',Inf);imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    %     h = ~bwmorph(~h,'bridge',Inf);imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    %         %     h = bwmorph(h,'erode',1);imshow(label2rgb(bwlabel(h),'hsv','k','shuffle'));
    %     h = ~bwmorph(~h,'bridge',Inf);imshow(label2rgb(bwlabel(h,4),'hsv','k','shuffle'));
    %
    bw = bwlabel(h,4);
    L4 = bw;
    
    
    
    %% Schnitz!
    
    %     h = bwmorph(f,'hbreak');imshowlabel(bwlabel(h));
    %     % h = bwmorph(f,'remove');imshowlabel(bwlabel(h));
    %     h = bwmorph(h,'open',1);imshowlabel(bwlabel(h));
    %     h = bwmorph(h,'erode',1);imshowlabel(bwlabel(h));
    %     h = bwmorph(h,'open',1);imshowlabel(bwlabel(h));
    %     h = bwmorph(h,'thicken');imshowlabel(bwlabel(h));
    
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
        %         cutcell= breakcellfluor(Lcell, imcomplement(combined), ...
        %             0.2, 0.005, p.minCellLength);
        cutcell= breakcellfluor(Lcell, imcomplement(phase), ...
            0.2, 0.04, p.minCellLength);
        
        L7(L6 == fbiggies(z))= 0;
        % place cutcell
        label= max2(L7);
        for j = 1:max2(cutcell),
            L7(cutcell==j)= label+j;
        end;
    end;
    L7= renumberimage(L7);
    imshowlabel(L7);
    
    
    %% Filling each cell
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
    
    L8 = L7;
    rgb3 = label2rgb(L8,'hsv','k','shuffle'); imshow(rgb3);
    %     imwrite(rgb3,['Final',num2str(i),'.jpg']);
    ns
    
end

