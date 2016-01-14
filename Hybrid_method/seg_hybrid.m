
% This can run the hybrid method outside Schnitzcells.
% It is purely for testing out on my Laptop, reading in the image files
% from my pendrive and then calling HYBRIDSEG.M


% Which  Movie set
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
    P0 = uint16(p0); %converts to integer
    T0 = uint16(t0); %converts to integer
    
    L9 =  hybridseg(P0,T0);
    
    rgb3 = label2rgb(L9,'hsv','k','shuffle'); imshow(rgb3);
    %     imwrite(rgb3,['Final',num2str(i),'.jpg']);
    ns
end

