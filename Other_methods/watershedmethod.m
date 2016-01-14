% The WATERSHED method
% This is currently written to run from my laptop, accessing the four
% movies Movie 1,2,3 and 4 from my pendrive. 
%
% Will need to edit the script to access files from another PC. 
% To run, call the script WATERSHEDMEHTOD

% This method :
% 1. Runs through the images and calculates the median
% 2. Normalises by the median
% 3. Standard deviation filter
% 4. Combine the RFP and Phase image by multiplication
% 5. Adaptive wiener filter
% 6. Gamma contrast correcton
% 7. Global threshold
% 8. Morphology
% 9. Removing small regions
% 10. Watershed
% 11. Improving the shape


% Bruno = input('Bruno? ');
Bruno = 1; %To run on Bruno's image files

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


p.minCellArea = 100;
p.maxThresh = 0.050;
p.minThresh = 0.050 ;
p.minCellLength = 20;


%%
medallp = 0;
medallt = 0;
iteration = 0;

for i = segRange
    
    n = 1000+i;
    ns = num2str(n);
    ns(1) = [];
    
    % On Manoj's Laptop [Bruno files]
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
    
    p0 = double(P0);
    t0 = double(T0);
    
    medip = median(p0(:)); %Median in image i
    medit = median(t0(:));  %Median in image i
    
    medallp = ((i-1)*medallp + medip)/i; %Median over all images
    medallt = ((i-1)*medallt + medit)/i; %Median over all images
end


%%
for i =segRange
    
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
                
                P0 = imread(['Y:\Bruno\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-13_1-p-',ns,'.TIF']);
                T0 = imread(['Y:\Bruno\2015-03-04-Bruno-Cyano2\SubFramedMovies\PSB7-13_1-t-',ns,'.TIF']);
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
    
    % 2. Normalise by median
    p0 = double(P0);
    t0 = double(T0);
    
%    p0 = p0*medallp/median(p0(:));
    P0 = uint16(p0); %converts to integer
    
 %   t0 = t0*medallp/median(t0(:));
    T0 = uint16(t0); %converts to integer
    
    % 3. Standard deviation filter
    g0 = stdfilt(P0,ones(3));g1 = stdfilt(T0,ones(3));
    g0 = g0/maxmax(g0)*255;g1 = g1/maxmax(g1)*255;
    
    % 4. Multiplying images
    x0 = uint16(g1).*uint16(g0);
    
    %Using only phase
    %     x0 = uint16(g0);
    %           x0 = P0;
    %
    %Using only RFP
    %     x0 = uint16(g1);
    %           x0 = T0;
    
    
    %% Manoj Edit %----------------------------------------------
    
    % 5. Wiener filter
    x_filt = uint16(wiener2(x0,[3 3])); imshow(x_filt);
    %     x_dilate = imdilate(x_filt,strel('disk',1));imshow(x_dilate);
    %     x_close = imerode(x_dilate,strel('disk',2));imshow(x_close);
    
    % 6. Adaptive histogram correction
    x_adapthist = adapthisteq(x_filt,'NumTiles',[8 8],'Distribution','exponential','Alpha',1);%imshow(M);
    % [x,counts]=imhist(x_adapthist);stem(x,counts);
    
    % 7. Gamma contrast correction
    x_gammacontrast = imadjust(x_adapthist,[],[],0.8); imshow(x_gammacontrast);
%     x_gammacontrast = imadjust(x_adapthist,stretchlim(x_adapthist),[0 0.7]);imshow(x_gammacontrast);

       
    % 8. Automatic thresholding
    x_thresh = im2bw(x_gammacontrast,graythresh(x_gammacontrast));imshow(x_thresh);
    %     J1 = imopen(J,strel('disk',1));imshow(J1);
    
    % 9. Closing small gaps (Morphological)
    x_close =  imclose(x_thresh,strel('disk',3));
    
    % Preliminary result
    rgb2 = label2rgb(bwlabel(~x_close),'jet',[.5 .5 .5]); imshow(rgb2);
    %     rgb2 = label2rgb(watershed(J1),'jet',[.5 .5 .5]); imshow(rgb2);
    %     imwrite(rgb1,['gamma',num2str(i),'.jpg']);
    
    
    y0(:,:,i) = x_close; %Storing results in a matrix

    x_inverse = ~x_close; % Image inverse
    
    %% Remove small black regions before watershed (Schnitzcells)
    r = regionprops(x_inverse,'Area');
    flittle = find([r.Area]<50);
    L(flittle)=[];
    L = bwlabel(x_inverse);
    imshow(L);
    goodones = setdiff(1:max(max(L)), flittle);
    bw2 = ismember(L, goodones);
    L2 = bwlabel(bw2);
    
    % C = imclose(~L2,strel('disk',2));imshow(C);
    
    %Removing pixel noise
    L3 = bwareaopen(~L2,500,4);imshow(L3); 
    % C = medfilt2(C,[2 2]);imshow(C);
    
    
    rgb2 = label2rgb(bwlabel(~L3),'jet', 'w', 'shuffle'); imshow(rgb2);
%     imwrite(rgb2,['gamma',num2str(i),'.jpg']);
    
    
    %% 10. Watershed
    W1 = watershed(L3);
    

    
    %% 11. Making watershed cells nicer by eroding each cell individually
    L4 = zeros(size(W1));
    for k = 1:max(max(W1))
        L2 = zeros(size(W1));
        L3 = zeros(size(W1));
        L2(W1==k) = W1(W1==k);
        %         L22 = imdilate(L2,strel('disk',1));
        L33 = imerode(L2,strel('disk',2));
        L4(L33==k) = L33(L33==k);
    end
    
    
    
    %% Saving Watershed output
    rgb3 = label2rgb(L4,'jet', 'w', 'shuffle');
    %     imshowpair(rgb,P0*50);
    imshow(rgb3);
%     imwrite(rgb3,['watershed_',num2str(i),'.jpg']);
    %        a = imfuse(rgb,P0*50);
    %        imwrite(a,['merge_',num2str(i),'.tif']);
    i
    %        pause
    
end

