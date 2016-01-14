% Used to produce the distribution of fluorescence in the field of view
% Example :
%
% Intensityheatmap
% When user input is requested 'What is the filename? '
% Key in eg Sbeads_001-1
%
% Produces a heat map with the following steps:

    % 1. Loading data from segmented beads (output file of REDSEGBEADS eg
    % Sbeads_001-1.mat)
    % 2. Calculating (i)mean intensities 
    %                (ii) number of beads 
    %                (iii) standard dev 
    %    within a box (box_height x box_width)
    % 3. Normalising the mean intensity
    % 4. Producing surf plots (viewed in 2D)
    % 5. Optional scaling of calculated mean intensities matrix to full image size matrix for future use
% - Manoj -

%% Loading data
clear;
filename = input('What is the filename? ','s');
    
%File names eg
    % Sqim_2_001-1


load(filename);
clear filename exptime ;
AA = who;
A =[];

%Putting data into A
for i = 1:length(AA)-1
    A = [A , eval(sprintf([AA{i}])) ];
end

var = eval(sprintf([AA{1}]));
[Image_height, Image_width ] = size(var); %Final image dimensions


 %% Copying all stats we need
 var = eval(sprintf([AA{length(AA)}]));
 rfp_means = zeros(1,length(var));
 centroids_y = zeros(1, length(var));
 centroids_x = zeros(1, length(var));
 
 for i = 1:length(var)
     rfp_means(i) = var(i).mean(:,2); % Change 1 = YFP and 2 = RFP
     centroids_y(i) = var(i).Centroid(:,2);
     centroids_x(i) = var(i).Centroid(:,1);
 end
 
 
 %% Calculating mean intensities within sections
 x_bins = 24;
 y_bins = 20;
 box_height = Image_height/(y_bins);   %Heat map sections
 box_width = Image_width/(x_bins); 
 mean_intensities = zeros(Image_height/box_height, Image_width/box_width); % For Means, Number, StDev
 counter_beads = zeros(Image_height/box_height, Image_width/box_width);
 std_beads = zeros(Image_height/box_height, Image_width/box_width);
 no_of_beads = length(var);
 
 xmin = 0; 
 xmax = box_width;
 ymin = 0;
 ymax = box_height;

 for j = 1:Image_width/box_width %Loop through across x
     for k = 1:Image_height/box_height %Loop through all y

         relevant_intensity = [];
         counter = 0;

         for m = 1:no_of_beads % Loop through all beads

           
            if(centroids_y(m) < ymax && centroids_y(m) >= ymin  ) %within box
                if(centroids_x(m) < xmax && centroids_x(m) >= xmin) %within box
                    relevant_intensity = [relevant_intensity,rfp_means(m)];
                    counter = counter +1;
                end
            end
            
         end

         % Calculating mean intensity in box
         
         if counter == 0 
                mean_intensities(k,j) = 0;
                counter_beads(k,j) = 0;
                std_beads(k,j) = 0;
            else
                mean_intensities(k,j) = mean(relevant_intensity); %within box mean intensity
                counter_beads(k,j) = counter; % number of beads per box
                std_beads(k,j) = std(relevant_intensity);
         end

         ymin = ymin + box_height;
         ymax = ymax + box_height;
         
%          figure(1); %Shows box 
%          axis([0 Image_width 0 Image_height]);
%          plot(xmin,ymin,'.',xmax,ymax,'.',xmin,ymax,'.',xmax,ymin,'.');
%          hold on;
     end
     
     ymin = 0;
     ymax = box_height;
     
     xmin = xmin + box_width;
     xmax = xmax + box_width;
 end

 
 %% Normalising
 mean_intensities = mean_intensities./max(max(mean_intensities));
 
  
 %% Producing the final map

fig = figure; colormap('jet');

subplot(2,2,1);
surf(mean_intensities);
view(2);
colorbar;
title('Mean Intensity');
xlabel('Bin number');
ylabel('Bin number');
axis([0,x_bins,0,y_bins]);


subplot(2,2,2);
surf(counter_beads);
view(2);
colorbar;
title('Number of beads');
xlabel('Bin number');
ylabel('Bin number');
axis([0,x_bins,0,y_bins]);

subplot(2,2,3);
surf(std_beads);
view(2);
colorbar;
title('Standard deviation');
xlabel('Bin number');
ylabel('Bin number');
axis([0,x_bins,0,y_bins]);

subplot(2,2,4);
plot(centroids_x,centroids_y,'o','MarkerSize',2,'Color',[0,255,127 ]./255);
title('Beads');xlabel('X (pixels)');ylabel('Y (pixels)');
axis([0,Image_width,0,Image_height]);

print(fig,'Heat-map','-dpng');
 % Binning outliers in the heat map to produce a better colour range for
 % the heat-map
% set(gca,'clim',[minval maxval]);


%% Full scale map
% This resizes the bins to the full size of all pixels in the image. This
% can be used to normalise the fluorescence to account for non-uniform
% fluorescence in the field of view. 

% fullsize = reshapematrix(mean_intensities,counter_beads,std_beads,box_height,box_width);
% figure;
% h = surf(fullsize(:,:,1));
% view(2);
% colorbar;
% title('Mean Intensity (Full scale)')           
% set(h,'LineStyle','none')

