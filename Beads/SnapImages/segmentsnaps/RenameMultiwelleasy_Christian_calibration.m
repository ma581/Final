%First the global parameters for the renaming are defined (numWells, numSubwells, numSubwells)

numWells = 1; %input('Enter number of wells visited: ');
numSubwells = 5; %input('Enter number of snapshots within each well: ');
outputFile =  'RFP30_YFP_100'; %input('Enter new file name, surrounded by single quotes: ');

% Then the paramters for each for loop are defined. And each for Loop
% renames the files in such a way that it is obvious from which pad the
% image is.

inputFile = 'awell-y-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'y'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end
  
inputFile = 'awell-t-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'r'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end

inputFile = 'awell-p-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'p'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end


%inputFile = 'awell-c-'; %input('Enter original file header, surrounded by single quotes: ');

%channel = 'c'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

%for i = 1 : numWells,
 %   for j = 1 : numSubwells,
  %      disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
   %     im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
    %    imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
   % end
%end

%inputFile = 'awell-g-'; 
%channel = 'g'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

%for i = 1 : numWells,
 %   for j = 1 : numSubwells,
  %      disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
   %     im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
    %    imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    %end
%end
%     for j = 1 : numFrames
%         im = imread([inputFile num2str(j) '.tif']);
%         imwrite(im,[outputFile str3(j) '.tif'],'tiff');
%         j
%     end
