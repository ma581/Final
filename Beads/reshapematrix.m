%% Reshaping matrix

function fullsize = reshapematrix(MatA,MatB,MatC,box_height,box_width)

     fullsize = zeros(1040,1392,1); %Full image size

     m = 1;
     n = 1;
     j = 1;
     k = 1;
     
     
     while j <= size(fullsize,2) %Loop through across x
         while k <= size(fullsize,1) %Loop through all y
             
                fullsize(k: k + box_height-1,j: j + box_width-1,1) = repmat(MatA(m,n),box_height,box_width);
                fullsize(k: k + box_height-1,j: j + box_width-1,2) = repmat(MatB(m,n),box_height,box_width);
                fullsize(k: k + box_height-1,j: j + box_width-1,3) = repmat(MatC(m,n),box_height,box_width);
                k = k + box_height;
                m = m + 1;
         end
         n = n + 1;
         
         m = 1;
         k = 1;
         j = j + box_width;
         
     end
         
% figure;
% h = surf(fullsize(:,:,1));
% view(2);
% colorbar;
% title('Mean Intensity (Full scale)')           
% set(h,'LineStyle','none')


end

% reshapematrix(mean_intensities,counter_beads,std_beads,box_height,box_width)
