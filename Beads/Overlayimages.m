
% This script overlays images to display beads identified and the good beads kept
% after processing


%% %----------------------------------------------------------------------
 % ALL BEADS
%----------------------------------------------------------------------
disp('Warning! This function clears the workspace!')

% Loading data
clear;
filename = input('What is the filename for ALL BEADS? ','s');
load(filename);
clear filename exptime Sbeads_001;
A = who; %Identifying all the processed images
I = []; %Dimensions of image

for i = 1:length(A)-1
    I = [I , eval(sprintf([A{i}])) ];
end
%I = [Lqbeads_00101, Lqbeads_00102, Lqbeads_00103, Lqbeads_00104, Lqbeads_00105, Lqbeads_00106, Lqbeads_00107, Lqbeads_00108, Lqbeads_00109, Lqbeads_00110, Lqbeads_00111, Lqbeads_00112, Lqbeads_00113, Lqbeads_00114, Lqbeads_00115, Lqbeads_00116, Lqbeads_00117, Lqbeads_00118, Lqbeads_00119, Lqbeads_00120];

% All coloured white
im_rgb1(:,:,1) = double(logical(I));
im_rgb1(:,:,2) = double(logical(I));
im_rgb1(:,:,3) = double(logical(I));
figure,imshow(im_rgb1);
 

%% %----------------------------------------------------------------------
 % Good beads
%----------------------------------------------------------------------

% Loading Data
clearvars -except im_rgb1
filename = input('What is the filename for GOOD BEADS? ','s');
load(filename);
clear filename exptime Sbeads_001;
% J = [Lqqim_2_00101, Lqqim_2_00102, Lqqim_2_00103, Lqqim_2_00104, Lqqim_2_00105, Lqqim_2_00106, Lqqim_2_00107, Lqqim_2_00108, Lqqim_2_00109, Lqqim_2_00110, Lqqim_2_00111, Lqqim_2_00112, Lqqim_2_00113, Lqqim_2_00114, Lqqim_2_00115, Lqqim_2_00116, Lqqim_2_00118, Lqqim_2_00119, Lqqim_2_00120, Lqqim_2_00121, Lqqim_2_00122, Lqqim_2_00123, Lqqim_2_00124, Lqqim_2_00125, Lqqim_2_00126, Lqqim_2_00127, Lqqim_2_00128, Lqqim_2_00129];
% J = [Lqbeads_00101, Lqbeads_00102, Lqbeads_00103, Lqbeads_00104, Lqbeads_00105, Lqbeads_00106, Lqbeads_00107, Lqbeads_00108, Lqbeads_00109, Lqbeads_00110, Lqbeads_00111, Lqbeads_00112, Lqbeads_00113, Lqbeads_00114, Lqbeads_00115, Lqbeads_00116, Lqbeads_00117, Lqbeads_00118, Lqbeads_00119, Lqbeads_00120];
A = who; %Identifying all the processed images
J =[]; 

for i = 1:(length(A)-2)
    J = [J , eval(sprintf([A{i}])) ];
end

% All coloured green
im_rgb2(:,:,1) = zeros(size(J));
im_rgb2(:,:,2) = double(logical(J));
im_rgb2(:,:,3) = zeros(size(J));
    figure,imshow(im_rgb2);

%% %----------------------------------------------------------------------
 % Merge
%----------------------------------------------------------------------

C = imfuse(im_rgb1,im_rgb2,'blend');
j = 1;

% D = zeros(1040,1392,3,length(A)-2);
for i = 1:length(A)-2
    D(:,:,:,i) = C(:,j:j+1391,:);
    j = j + 1392;
    figure,imshow(D(:,:,:,i));
    imwrite(D(:,:,:,i),['im_',num2str(i),'.tif']);
end


% %% %----------------------------------------------------------------------
%  % Small beads
% %----------------------------------------------------------------------
% 
% load('Sqim_2_001-1 - Smallbeads.mat');
%  J = [Lqqim_2_00101, Lqqim_2_00102, Lqqim_2_00103, Lqqim_2_00104, Lqqim_2_00105, Lqqim_2_00106, zeros(size(Lqqim_2_00106)),Lqqim_2_00108, Lqqim_2_00109, Lqqim_2_00110, Lqqim_2_00111, Lqqim_2_00112, Lqqim_2_00113, Lqqim_2_00114, zeros(size(Lqqim_2_00106)),Lqqim_2_00116, zeros(size(Lqqim_2_00106)),zeros(size(Lqqim_2_00106)),Lqqim_2_00120, Lqqim_2_00121, Lqqim_2_00122, Lqqim_2_00123, Lqqim_2_00124, Lqqim_2_00125, Lqqim_2_00126, Lqqim_2_00127, Lqqim_2_00128, Lqqim_2_00129];
% 
%  % All coloured red
% im_rgb3(:,:,1) = double(logical(J));
% im_rgb3(:,:,2) = zeros(size(J));
% im_rgb3(:,:,3) = zeros(size(J));
% figure,imshow(im_rgb3);
% 
% E = imfuse(C,im_rgb3,'blend');
% j = 1;
% for i = 1:28
%     F(:,:,:,i) = E(:,j:j+1391,:);
%     j = j + 1392;
% %     figure,imshow(D(:,:,:,i));
%     imwrite(F(:,:,:,i),['im_',num2str(i),'.tif']);
% end


%% %----------------------------------------------------------------------
 % Low intensity beads
%----------------------------------------------------------------------

% load('Sqim_2_001-1 - Lowintensity.mat');
%  J = [zeros(size(Lqqim_2_00104)),zeros(size(Lqqim_2_00104)),zeros(size(Lqqim_2_00104)),...
%      Lqqim_2_00104, zeros(size(Lqqim_2_00104)), Lqqim_2_00106, Lqqim_2_00107, zeros(size(Lqqim_2_00104)),...
%      Lqqim_2_00109, Lqqim_2_00110, Lqqim_2_00111, Lqqim_2_00112, Lqqim_2_00113, Lqqim_2_00114, Lqqim_2_00115,zeros(size(Lqqim_2_00104)), ...
%      zeros(size(Lqqim_2_00104)),zeros(size(Lqqim_2_00104)),zeros(size(Lqqim_2_00104)),...
%      Lqqim_2_00121, Lqqim_2_00122, Lqqim_2_00123, Lqqim_2_00124, Lqqim_2_00125, Lqqim_2_00126, Lqqim_2_00127, Lqqim_2_00128, zeros(size(Lqqim_2_00104))];
%  % All coloured red
% im_rgb3(:,:,1) = double(logical(J));
% im_rgb3(:,:,2) = zeros(size(J));
% im_rgb3(:,:,3) = zeros(size(J));
% % figure,imshow(im_rgb3);
% 
% E = imfuse(C,im_rgb3,'blend');
% j = 1;
% for i = 1:28
%     F(:,:,:,i) = E(:,j:j+1391,:);
%     j = j + 1392;
% %     figure,imshow(F(:,:,:,i));
%     imwrite(F(:,:,:,i),['im_',num2str(i),'.tif']);
% end


%% Aside stuff

% srcFiles = dir('C:\Users\Manoj\Dropbox\SainsburyLab\Beads\output\*.tif');  % the folder in which ur images exists
% for i = 1 : length(srcFiles)
%     filename = strcat('C:\Users\Manoj\Dropbox\SainsburyLab\Beads\output\',srcFiles(i).name);
%     I = imread(filename);
%     figure, imshow(I);
% end


% All Beads
% input_data = load('Sqim_2_001-1 -Allbeads.mat');
% srcFiles = dir('C:\Users\Manoj\Dropbox\SainsburyLab\Beads\output\*.tif');  % the folder in which ur images exists
% for i = 1 : length(srcFiles)
%     filename = strcat('C:\Users\Manoj\Dropbox\SainsburyLab\Beads\output\',srcFiles(i).name);
%     I = [Lqqim_2_00101, Lqqim_2_00102, Lqqim_2_00103, Lqqim_2_00104, Lqqim_2_00105, Lqqim_2_00106, Lqqim_2_00107, Lqqim_2_00108, Lqqim_2_00109, Lqqim_2_00110, Lqqim_2_00111, Lqqim_2_00112, Lqqim_2_00113, Lqqim_2_00114, Lqqim_2_00115, Lqqim_2_00116, Lqqim_2_00118, Lqqim_2_00119, Lqqim_2_00120, Lqqim_2_00121, Lqqim_2_00122, Lqqim_2_00123, Lqqim_2_00124, Lqqim_2_00125, Lqqim_2_00126, Lqqim_2_00127, Lqqim_2_00128, Lqqim_2_00129];
% 
%     % All coloured white
%     im_rgb1(:,:,1) = logical(I);
%     im_rgb1(:,:,2) = logical(I);
%     im_rgb1(:,:,3) = logical(I);
% end
