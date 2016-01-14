% The MARKER WATERSHED method
% This is currently written to run from my laptop, accessing the four
% movies Movie 1,2,3 and 4 from my pendrive. 
%
% Will need to edit the script to access files from another PC. 
% To run, call the script MARKER_WATERSHED

% This method :
% 1. Calls markers


% Bruno = input('Bruno? ');
set = input('Which Image set? ');
Laptop = 1;
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
medallp = 0;
medallt = 0;
iteration = 0;



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
    
    p0 = double(P0);
    t0 = double(T0);
    
%    p0 = p0*medallp/median(p0(:));
    P0 = uint16(p0); %converts to integer
    
 %   t0 = t0*medallp/median(t0(:));
    T0 = uint16(t0); %converts to integer
    
    g0 = stdfilt(P0,ones(3));g1 = stdfilt(T0,ones(3));
    g0 = g0/maxmax(g0)*255;g1 = g1/maxmax(g1)*255;
    x0 = uint16(g1).*uint16(g0);
    
    %Using only phase
    %     x0 = uint16(g0);
    %           x0 = P0;
    %
    %Using only RFP
    %     x0 = uint16(g1);
    %           x0 = T0;
    
    
    %% Manoj Edit %----------------------------------------------

    
    % Calling sub-script to run the marker watershed method
    markers;
    rgb2 = label2rgb(L,'jet',[.5 .5 .5]); imshow(rgb2);
    imwrite(rgb2,['Markerwatershed',num2str(i),'.jpg']);

    
end

